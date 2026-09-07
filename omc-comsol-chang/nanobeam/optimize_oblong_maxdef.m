function [best_oblong, best_maxdef, fitnessBest, dsBest] = optimize_oblong_maxdef()
%OPTIMIZE_OBLONG_MAXDEF  Optimize the nanobeam cavity defect taper.
%
%   Searches two parameters, where each objective evaluation is a full COMSOL
%   mechanical + optical + PML solve via RunNanobeamFEM:
%       x(1) = defectAspectRatio  (mapped to P.oblong, see below)
%       x(2) = maxdef             (maximum defect ratio, dimensionless)
%
%   This is the optimizer counterpart of sweep_oblong_maxdef.m: instead of
%   evaluating a fixed grid, it lets an optimizer choose where to spend the
%   evaluation budget.  See optimization_README.md for the full guide.
%
%   OPTIMIZER (OPT.optimizer)
%     'neldermead' - fminsearch simplex.  Base MATLAB, works today.
%     'bayesopt'   - Gaussian-process Bayesian optimization.  RESERVED: the
%                    branch raises a specific error explaining that the
%                    Statistics and Machine Learning Toolbox is licensed but
%                    not installed.  Design: bayesopt_oblong_maxdef_PROPOSAL.md
%
%   FITNESS (OPT.fitnessFcn, default fitness_coopProduct)
%
%       fitness = g_OM^2 * Q_opt * exp(-((lambda_opt - targetLambda)/lambdaTol)^2)
%
%   subject to a BINARY mechanical-Q acceptance gate:
%
%       Q_mech <  OPT.QmechMin  ->  fitness = 0   (design rejected outright)
%
%   The wavelength factor is a GAUSSIAN in the detuning, so lambdaTol is the
%   1/e half-width.  See fitness_coopProduct for the full definition and for
%   the alternative objectives kept alongside it.
%
%   PARAMETER MAPPING (identical to the sweep)
%       mirrorAspectRatio = P.hy / P.hx;
%       r = defectAspectRatio / mirrorAspectRatio;
%       Pe.oblong = log(r) / (2*log(1 - maxdef));
%
%   OBJECTIVE SCALING
%   Optimizers minimize, and the raw fitness spans many orders of magnitude,
%   so the objective returned is -log10(fitness).  That keeps optimset's
%   absolute TolFun meaningful (it becomes an effectively relative test on the
%   fitness) instead of unreachable.  Nothing can abort a run: every
%   evaluation is wrapped in try/catch, and the infeasible outcomes are
%   penalized in increasing order of badness --
%       gated (Q_mech too low)        ->  1 * OPT.Jpenalty
%       unusable (no mode / error)    ->  2 * OPT.Jpenalty
%       out of bounds                 ->  OPT.Jpenalty * (2 + distance)
%   The out-of-bounds penalty grows with how far outside the box the point
%   lies, so the simplex always has a direction back into the feasible set.
%
%   STOP AND RESUME (OptimStore.m)
%   Every evaluation is written to a SQLite database the moment it completes,
%   so an interrupted study never loses a COMSOL solve.  Because fminsearch is
%   deterministic given x0, a resumed session retraces its earlier path
%   exactly -- and those repeats are served from the database instead of being
%   re-solved, turning the replay into seconds rather than hours.
%     OPT.runId  : names the study.  KEEP it to continue a study, CHANGE it to
%                  start a fresh one.
%     OPT.resume : 1 = serve stored points, 0 = re-solve (still records).
%     OPT.dbPath : the database, shared across studies.
%   Ctrl-C is safe: an onCleanup guard marks the run 'interrupted', and the
%   completed evaluations are already persisted.
%
%   Returns:
%     best_oblong  - oblong value at the optimum (mapped from defectAspectRatio)
%     best_maxdef  - maxdef value at the optimum
%     fitnessBest  - fitness (positive, un-scaled) at the optimum
%     dsBest       - ds struct from the final re-run of the best point
%
%   Writes into OPT.rootLoc (which is keyed on OPT.runId, not the date, so a
%   resumed study keeps all of its output together):
%     eval_NNNN_obX_mdY/           per-evaluation FEM output and mode profiles
%     BEST/                        final re-run of the best design, plots on
%     optimize_log.csv             one row per evaluation, appended
%     optimize_convergence.png/fig fitness trace + simplex trajectory
%   and into OPT.dbPath: the SQLite run database, plus a .jsonl mirror.
%
%   See also OPTIMSTORE, TEST_OPTIMSTORE_RESUME, SWEEP_OBLONG_MAXDEF,
%   RUNNANOBEAMFEM.
%
%   Based on sweep_oblong_maxdef.m; adapted to use an optimizer.

close all;

%% ===================== TUNABLE KNOBS =====================

% --- nominal mirror hole dimensions (m) ---
% Defined here rather than only in the P struct below because the starting
% point and the objective's oblong mapping BOTH depend on the mirror aspect
% ratio hy/hx.  If the two are written out separately they silently drift
% apart, and x0 then no longer corresponds to oblong_0 at all.  P.hx / P.hy
% are assigned from these.
hx_0 = 343e-9;      % nominal mirror hole height
hy_0 = 617e-9;      % nominal mirror hole width

% --- starting point (same as sweep centre) ---
% hy ~ (1-maxdef)^(1+oblong), hx ~ (1-maxdef)^(1-oblong).
% oblong=0 -> equal scaling; oblong=1 -> hx constant, only hy shrinks.
oblong_0 = 1.15;
maxdef_0 = 0.16;
defectAspectRatio_0 = (hy_0/hx_0) * (1-maxdef_0)^(2*oblong_0);  % hy/hx * (1-maxdef)^(2*oblong)

% --- defectAspectRatio box constraints (dimensionless) ---
OPT.defectAspectRatio_min = defectAspectRatio_0*0.5;
OPT.defectAspectRatio_max = defectAspectRatio_0*1.5;

% --- maxdef box constraints (dimensionless, 0–<1) ---
% Fraction by which the mirror lattice constant is reduced at cavity centre.
% Lower bound is an ABSOLUTE 0.05, not a fraction of maxdef_0, so the box no
% longer moves if the start point is retuned.  Note the box is now asymmetric
% about x0 (0.05 .. 0.33 around 0.22), which is fine for both optimizers --
% neither assumes a centred box -- but it widens the INDUCED oblong range a
% long way, because oblong = log(r)/(2*log(1-maxdef)) diverges as maxdef -> 0.
% See the feasibility scan printed at startup.
OPT.maxdef_min = 0.05;
OPT.maxdef_max = maxdef_0*1.5;

% --- target mechanical frequency (used for solver target) ---
OPT.targetFreq = 7e9;  % Hz

% --- optical-wavelength fitness target and tolerance ---
% The wavelength factor is GAUSSIAN: exp(-((lambda - target)/lambdaTol)^2),
% so lambdaTol is the 1/e half-width, NOT an exponential decay length.
% A lambdaTol miss costs a factor 1/e; a 2*lambdaTol miss costs 1/e^4.
OPT.targetLambda = 1550;   % nm  — target optical resonance wavelength
OPT.lambdaTol    = 100;    % nm  — Gaussian 1/e half-width for wavelength mismatch

% --- mechanical-Q acceptance gate ---
% Q_mech enters the fitness as a hard pass/fail rather than as a smooth
% factor: any design whose radiative Q_mech falls below this scores exactly 0
% and is rejected.  Requires P.solveMechPML = 1 (otherwise Q_mech is NaN and
% the point is treated as unusable rather than gated).
OPT.QmechMin = 1e7;

% --- objective penalty scale ---
% Objective values for feasible designs are -log10(fitness), i.e. roughly
% -10 to -18, so any penalty comfortably above zero dominates them.  Keeping
% the penalty O(1e3) rather than O(1e12) avoids reintroducing the enormous
% dynamic range that the log scaling exists to remove.
OPT.Jpenalty = 1e3;

% --- optimizer selection ---
% 'bayesopt'   : Gaussian-process Bayesian optimization.  Spends the
%                evaluation budget where the surrogate says it is most
%                informative, models the Q_mech gate as a constraint rather
%                than a penalty, and handles the flat gated region that
%                stalls a simplex.  Needs the Statistics and Machine
%                Learning Toolbox.  Settings under OPT.bo below.
% 'neldermead' : fminsearch simplex (base MATLAB, no toolbox).  Kept as a
%                switchable alternative -- both write the same physics to
%                the same run database, so a study started with one can be
%                resumed or seeded with the other.
OPT.optimizer = 'bayesopt';

% --- run persistence / stop-and-resume (see OptimStore.m) -----------------
% Every evaluation is written to SQLite as it completes, so an interrupted
% study never loses a COMSOL solve.  Re-running with the SAME OPT.runId and
% OPT.resume = 1 replays the stored points for free, then carries on.
%   OPT.runId  : names the study.  Change it to start a fresh search, keep
%                it to continue one.  Deliberately NOT date-stamped -- a
%                date-stamped id would silently start a new study tomorrow.
%   OPT.resume : 1 = serve already-evaluated points from the database,
%                0 = re-solve everything (the database is still written).
%   OPT.storeBackend : which storage layer to use.  '' auto-detects; force
%                one of 'dbtoolbox' | 'python' | 'cli' | 'jsonl' to pin it.
%                An unknown name is an error, never a silent fall-through.
%                OptimStore.availableBackends() lists what this machine has,
%                and OptimStore.selfTestAll() exercises each of them.
OPT.useStore     = 1;

% PINNED to 'jsonl' deliberately, not left on auto-detect.  The remote
% workstation has the Microsoft Store build of Python, which MATLAB cannot
% launch at all: the ACLs on %ProgramFiles%\WindowsApps deny process creation,
% so every probe costs a failed spawn plus an OS-level "Access is denied"
% message.  Pinning skips detection entirely.
%
% Nothing is given up by this.  'jsonl' needs only a writable directory, and
% test_OptimStore_resume runs its whole suite against EVERY backend and
% asserts they produce identical numbers -- runs, resume, best, listRuns and
% GP seeding all behave the same.  What is unavailable is ad-hoc SQL querying
% of the results.
%
% Set back to '' once SQLite is reachable everywhere this runs (a python.org
% interpreter via pyenv, or a sqlite3 executable on the PATH).
OPT.storeBackend = 'jsonl';
OPT.runId    = 'oblongMaxdef_trial1';
OPT.dbPath   = fullfile('.', 'test', '1D_OMC_hole', 'optim_runs.sqlite3');
OPT.resume   = 1;

% --- root folder for per-evaluation subfolders ---
% Keyed on OPT.runId, deliberately NOT on the date.  A study resumed on a
% later date must write its new eval folders, CSV and convergence figure into
% the SAME place as the session that started it -- a date stamp here scatters
% one logical study across several dated directories, while the database
% (OPT.dbPath) keeps treating it as a single run.
OPT.rootLoc = [fullfile('.','test','1D_OMC_hole',['optimize_',OPT.runId]),filesep];

% --- FITNESS FUNCTION (user-configurable) ---------------------------------
% Cooperativity-like figure of merit with a binary Q_mech gate:
%   Q_mech >= OPT.QmechMin:
%       fitness = g_OM^2 * Q_opt * exp(-((lambda_opt - targetLambda)/lambdaTol)^2)
%   Q_mech <  OPT.QmechMin:
%       fitness = 0
%   Q_mech unavailable (PML off / mech solve failed):
%       fitness = NaN
% Alternative objectives (fitness_QxQxLambda, fitness_maxGOM,
% fitness_maxLSiV) are defined at the bottom of this file -- swap one in here.
OPT.fitnessFcn = @(ds) fitness_coopProduct(ds, OPT.targetLambda, OPT.lambdaTol, OPT.QmechMin);

% The SAME figure of merit with the Q_mech gate disabled.  Nelder-Mead uses
% the gated score (a rejected design is simply worthless to it), but a GP
% surrogate needs the ungated value: a gated 0 is the same number everywhere
% below the threshold and tells the model nothing about WHERE the boundary
% is.  Under 'bayesopt' the gate becomes a coupled constraint instead, and
% this is the objective.  Both are stored, so the two optimizers can share a
% run database and resume from each other.
OPT.fitnessFcnRaw = @(ds) fitness_coopProduct(ds, OPT.targetLambda, OPT.lambdaTol, -Inf);

% --- bayesopt settings (used only when OPT.optimizer = 'bayesopt') --------
% maxEvals is the TOTAL size of the study, and points seeded from the store
% COUNT AGAINST IT -- bayesopt's MaxObjectiveEvaluations includes InitialX.
% So resuming a study that already holds maxEvals points does nothing: to
% extend it, raise maxEvals, exactly as resume() would.  Every point not
% served from the store is a full COMSOL solve, so time one first.
OPT.bo.maxEvals       = 60;
OPT.bo.numSeedPoints  = 8;
OPT.bo.acquisition    = 'expected-improvement-plus';
% FALSE on purpose: adaptive meshing (P.mAdjMesh / P.oAdjMesh) makes a
% repeated evaluation of the same geometry differ slightly.  Declaring the
% objective deterministic would force the GP to interpolate that jitter
% exactly and produce a badly conditioned surrogate.
OPT.bo.deterministic  = false;
OPT.bo.plotFcn        = {@plotMinObjective, @plotObjectiveModel};

%% ===================== P STRUCT DEFAULTS =====================
% Geometry matches test_nanobeamRectFEM_withPML.m (fabricated OMC device).
P.xsect     = 'rect';     % beam cross sectional shape - 'tri' or 'rect' or 'isoFit'
P.beamMat   = 'diamond';
P.celltype  = 'hole';
P.anisoMat  = 1;
P.rxtal     = 0;
P.rxtalInFilename = 1;

P.a   = 650e-9;     % nominal lattice constant (m)
P.w   = 800e-9;     % beam width (m)
P.theta = 45;       % etch angle (degrees; no effect for isoFit/rect)
P.th  = 250e-9;     % beam thickness (m)
P.hx  = hx_0;       % nominal mirror hole height (m) - set in TUNABLE KNOBS
P.hy  = hy_0;       % nominal mirror hole width  (m) - set in TUNABLE KNOBS

P.nholes    = 18;   % number of holes in half-beam
P.ndef      = 8;    % number of holes in half-defect region
P.maxdef    = maxdef_0;  % start value (22 %) — overwritten every evaluation
P.oblong    = oblong_0;  % start value — overwritten every evaluation
P.taperFunc = 'cubic';
P.holeatctr = 1;    % hole at cavity centre (matches fabricated device)

P.useManualDefect = 0;

P.wvgmir          = 0;
P.wgmTaper.func   = 'cubic';
P.wgmTaper.endtype = 'custom';
P.wgmTaper.a_end   = 650e-9;
P.wgmTaper.hx_end  = 343e-9;
P.wgmTaper.hy_end  = 300e-9;

P.asymCav = 0;      % disable asymmetric cavity for speed during optimization

P.lambda = 1550e-9;
P.nbeam  = 2.386;

P.stdDev    = [0, 0];
P.stdDevPos = 0;
P.asym      = 0;

% --- physics toggles (mechanical + optics + PML) ---
P.solveMech    = 1;
P.solveOpt     = 1;
P.calcG        = 1 * (P.solveMech && P.solveOpt);
P.calcS        = 0;
P.solveMechPML = 1;     % mechanical PML for radiative Q estimation

% Filter out non-localised modes; RunNanobeamFEM records warning when none pass.
P.extractLocMechModes = 1;

% --- plot / save: mode profiles saved to each per-eval subfolder ---
P.plotgeom   = 0;
P.storeMPH   = 0;
P.plotMech   = 1 * P.solveMech;  % displacement/strain plot per eval
P.plotOpt    = 1 * P.solveOpt;
P.plotStrCpl = 1 * P.calcS;
P.saveplots  = 1;   % write .png/.fig into the per-eval evLoc subfolder
P.savedat    = 1;

% --- mechanical solver ---
P.mevenx  = 1;
P.meveny  = 1;
P.mevenz  = 1;
P.freq    = OPT.targetFreq;
P.mneigs  = 20;
P.mMesh   = 1;     % fine mesh (matches withPML test script)
P.mAdjMesh = 1;
P.max_dof = 5e6;

% --- optical solver (used only if solveOpt = 1) ---
P.oevenx  = (-1)^P.holeatctr;  % -1 for fundamental mode with hole at centre
P.oeveny  = -1;
P.oevenz  = 1;
P.oneigs  = 1;
P.oMesh   = 3;
P.oAdjMesh = 1;
P.airrad  = 2 * P.lambda + P.w / 2;

P.g0min = 80e3;

% SiV strain coupling axes
P.zSiV  = {[1 1 1]};
P.xSlc  = 0;
P.ySlc  = 0;
P.zSlc  = 0;
P.LStats.xmin = [];
P.LStats.xmax = [];
P.LStats.ymin = 0;
P.LStats.ymax = 60e-9;
P.LStats.zmin = P.th / 2 - 80e-9;
P.LStats.zmax = P.th / 2;

% PML settings match test_nanobeamRectFEM_withPML.m
P.PMLmeshDiv = 20;
P.PMLLen     = 10e-6;
P.PMLstr     = 0.008;

%% ===================== OPEN CSV LOG =====================
% APPEND, never truncate.  OPT.rootLoc is keyed on OPT.runId rather than
% the date, so a resumed session reopens the SAME file as the session that
% started the study -- opening it 'w' would silently erase the history of
% every earlier session.  The header is written only when the file is being
% created, so a resumed log stays a single valid CSV.
if ~exist(OPT.rootLoc, 'dir'); mkdir(OPT.rootLoc); end
logFile = [OPT.rootLoc, 'optimize_log.csv'];
newLog  = (exist(logFile, 'file') ~= 2);
fid = fopen(logFile, 'a');
if newLog
    fprintf(fid, ['eval,oblong,maxdef,wM_GHz,gOM_kHz,gMB_kHz,gPE_kHz,' ...
                  'LSiV_MHz,Q_mech,Q_opt,lambda_opt_nm,fitness,status\n']);
else
    % Mark the seam so a multi-session log can be read back unambiguously.
    fprintf(fid, '# resumed %s\n', ...
        char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
end
fclose(fid);

% Eval counter (shared with the objective via nested scope).
evalCount = 0;

% Track the best (positive) fitness seen so far for the final re-run.
bestSoFar.fitness = -Inf;
bestSoFar.oblong  = NaN;
bestSoFar.maxdef  = NaN;
bestSoFar.defectAspectRatio = NaN;

% Full search history, appended by the objective and plotted after the run.
history.eval    = [];
history.dAR     = [];
history.maxdef  = [];
history.oblong  = [];
history.fitness = [];
history.Qmech   = [];
history.J       = [];

% Number of evaluations served from the database instead of from COMSOL.
nCacheHits = 0;

% Number of candidates rejected by the lithography pre-filter (no solve).
nInfeasible = 0;

%% ===================== RUN STORE (SQLite stop / resume) ==================
store = [];
if OPT.useStore
    store = OptimStore(OPT.dbPath, OPT.storeBackend);
    store.startRun(OPT.runId, 'optimize_oblong_maxdef', OPT.optimizer, OPT);

    nPrior = store.evalCount(OPT.runId);
    fprintf('Run store    : %s  [backend: %s]\n', store.dbPath, store.backend);
    fprintf('Run id       : %s  (%d evaluation(s) already stored)\n', ...
        OPT.runId, nPrior);

    if nPrior > 0 && OPT.resume
        % Seed the incumbent from what the earlier session already found,
        % so a resumed run can never report a worse design than one that
        % has already been paid for -- even if every eval this session fails.
        [xPrior, fPrior, iPrior] = store.best(OPT.runId);
        if ~isempty(xPrior)
            bestSoFar.fitness = fPrior;
            bestSoFar.defectAspectRatio = xPrior(1);
            bestSoFar.maxdef  = xPrior(2);
            bestSoFar.oblong  = log(xPrior(1) / (P.hy / P.hx)) / ...
                                (2 * log(1 - xPrior(2)));
            fprintf(['Resuming     : best stored fitness %.4g at eval %d ' ...
                     '(dAR=%.4f, maxdef=%.4f)\n'], ...
                fPrior, iPrior, xPrior(1), xPrior(2));
        else
            fprintf('Resuming     : no successful evaluation stored yet\n');
        end
        % Keep new indices clear of the ones already in the database.
        evalCount = nPrior;
    elseif nPrior > 0
        fprintf(['NOT resuming : OPT.resume = 0, so all %d stored point(s) ' ...
                 'will be re-solved\n'], nPrior);
        evalCount = nPrior;
    end

    % Mark the run interrupted if this function exits without reaching the
    % end -- an error, or a Ctrl-C.  The close-out block at the bottom
    % overwrites it with 'complete' on the happy path.
    runGuard = onCleanup(@() store.finishRun(OPT.runId, 'interrupted'));
end

%% ===================== FMINSEARCH SETUP =====================
x0 = [defectAspectRatio_0, maxdef_0];

% TolFun is an ABSOLUTE tolerance on the spread of objective values across the
% simplex.  Because the objective is -log10(fitness) -- an O(10) quantity --
% rather than -fitness, which runs to ~1e18, TolFun here acts as an
% effectively RELATIVE tolerance on the fitness itself: 1e-3 in log10 is about
% a 0.23 % change.  Tightening it much further is pointless: FEM mesh
% adaptation (P.mAdjMesh / P.oAdjMesh) makes the objective slightly
% discontinuous as the geometry changes, so the numerical noise floor sits
% above a 1e-6 tolerance and the simplex would only chase mesh artifacts.
%
% CAVEAT: Nelder-Mead has no noise handling of its own.  If the run converges
% onto a point whose fitness does not reproduce when re-run, that is mesh
% noise, not an optimum -- coarsen the adaptive meshing or seed this script
% from a coarse sweep_oblong_maxdef pass instead of a cold start.
options = optimset('Display', 'iter', ...
                   'TolX',        1e-4, ...
                   'TolFun',      1e-3, ...
                   'MaxIter',     100, ...
                   'MaxFunEvals', 200);

fprintf('\n==================== OPTIMIZATION START ====================\n');
fprintf('Start defectAspectRatio = %.6f  (oblong_0 = %.4f)\n', defectAspectRatio_0, oblong_0);
fprintf('Start maxdef            = %.6f\n', maxdef_0);
fprintf('defectAspectRatio bounds = [%.6f, %.6f]\n', ...
    OPT.defectAspectRatio_min, OPT.defectAspectRatio_max);
fprintf('maxdef bounds            = [%.6f, %.6f]\n', ...
    OPT.maxdef_min, OPT.maxdef_max);
fprintf('Q_mech acceptance gate   = %.3g  (fitness = 0 below this)\n', OPT.QmechMin);
fprintf('Log written -> %s\n\n', logFile);

%% ===================== FABRICATION FEASIBILITY SCAN ==================
% isFabricable reproduces CreateNanobeamGeom's lithography guards as pure
% arithmetic, so an unfabricable candidate costs microseconds instead of a
% multi-minute solve that throws.  Scanning the box up front says how much
% of the search space is actually reachable -- a mostly-infeasible box is
% worth knowing about BEFORE spending a day of COMSOL on it.
nScan = 41;
[scanDAR, scanMD] = meshgrid( ...
    linspace(OPT.defectAspectRatio_min, OPT.defectAspectRatio_max, nScan), ...
    linspace(OPT.maxdef_min, OPT.maxdef_max, nScan));
[scanOK, scanMargin] = isFabricable(scanDAR, scanMD, P);
fracFab = mean(scanOK(:));
fprintf('Fabricable   : %.1f%% of the search box (tightest margin %.1f nm)\n', ...
    100*fracFab, min(scanMargin(:))*1e9);
if fracFab == 0
    error('optimize_oblong_maxdef:boxInfeasible', ...
        ['No point on a %dx%d scan of the bounds is fabricable, so the ' ...
         'study was NOT started. Widen the bounds, or relax ' ...
         'P.minHoleGap / P.minSidewallGap / P.minFeature.'], nScan, nScan);
elseif fracFab < 0.25
    warning('optimize_oblong_maxdef:boxMostlyInfeasible', ...
        ['Only %.1f%% of the search box is fabricable. The optimizer can ' ...
         'work here, but seed points are drawn from the whole box, so ' ...
         'most will be rejected. Consider tightening the bounds around ' ...
         'the feasible region.'], 100*fracFab);
end

% Wrap the objective in a try/catch so the optimizer never crashes.
objFcn = @(x) safeObjective(x);

switch lower(OPT.optimizer)
    case 'neldermead'
        [xBest, JBest] = fminsearch(objFcn, x0, options);

    case 'bayesopt'
        % Gaussian-process Bayesian optimization.  The check below is
        % deliberately specific: a bare 'Unrecognized function ''bayesopt'''
        % reads like a licence problem when it is usually a missing install,
        % and sends you debugging the wrong thing.
        if exist('bayesopt', 'file') == 0 && exist('bayesopt', 'builtin') == 0
            error('optimize_oblong_maxdef:noBayesopt', ...
                ['bayesopt is not available in this MATLAB install.\n' ...
                 '  licence permits it (1=yes) : %d\n' ...
                 '  installed on disk (1=yes)  : %d\n' ...
                 'Licence 1 with install 0 means: install the Statistics ' ...
                 'and Machine Learning Toolbox from the Add-On Explorer. ' ...
                 'No re-licensing is needed.'], ...
                license('test','Statistics_Toolbox'), ...
                isfolder(fullfile(matlabroot,'toolbox','stats')));
        end

        % The box is the variable range, so bayesopt never proposes an
        % out-of-bounds point and the sloped barrier is simply unused here.
        optVars = [ ...
            optimizableVariable('dAR', ...
                [OPT.defectAspectRatio_min, OPT.defectAspectRatio_max], ...
                'Type','real'), ...
            optimizableVariable('maxdef', ...
                [OPT.maxdef_min, OPT.maxdef_max], 'Type','real')];

        % Seed the GP with whatever the store already holds for this run, so
        % resuming a Bayesian study -- or picking one up from a Nelder-Mead
        % study that used the same runId -- costs no COMSOL time.
        [initX, initObj, initCon] = seedBayesFromStore(store, OPT);
        nSeeded  = height(initX);
        seedLeft = max(0, OPT.bo.numSeedPoints - nSeeded);
        if nSeeded > 0
            fprintf(['Seeding GP with %d stored point(s); %d fresh seed ' ...
                     'point(s) to go.\n'], nSeeded, seedLeft);
        end
        if nSeeded >= OPT.bo.maxEvals
            % Silently doing nothing is the worst outcome here: it looks
            % like a finished study that simply made no progress.
            warning('optimize_oblong_maxdef:budgetExhausted', ...
                ['The store already holds %d point(s) for run "%s", ' ...
                 'which meets or exceeds OPT.bo.maxEvals = %d. bayesopt ' ...
                 'counts seeded points against its budget, so NO new ' ...
                 'evaluations will be made. Raise OPT.bo.maxEvals to ' ...
                 'extend the study.'], ...
                nSeeded, OPT.runId, OPT.bo.maxEvals);
        end

        % XConstraintFcn is evaluated on CANDIDATES, before the objective is
        % called at all, so an unfabricable design is never proposed and
        % never costs a solve.  This is the cheapest possible rejection --
        % strictly better than the coupled constraints, which require an
        % evaluation to learn from.  isFabricable is vectorized because
        % bayesopt passes a table of many candidate rows at once.
        xConFcn = @(t) isFabricable(t.dAR, t.maxdef, P);

        boArgs = {'MaxObjectiveEvaluations',  OPT.bo.maxEvals, ...
                  'XConstraintFcn',           xConFcn, ...
                  'NumSeedPoints',            max(1, seedLeft), ...
                  'AcquisitionFunctionName',  OPT.bo.acquisition, ...
                  'IsObjectiveDeterministic', OPT.bo.deterministic, ...
                  'NumCoupledConstraints',    2, ...
                  'AreCoupledConstraintsDeterministic', [false false], ...
                  'PlotFcn',                  OPT.bo.plotFcn, ...
                  'Verbose',                  1};
        if nSeeded > 0
            boArgs = [boArgs, {'InitialX', initX, ...
                               'InitialObjective', initObj, ...
                               'InitialConstraintViolations', initCon}];
        end

        boResults = bayesopt(@bayesObjective, optVars, boArgs{:});

        % bayesopt has several notions of "best" (observed vs. estimated,
        % feasible vs. not).  Rather than pick one, take its feasible best as
        % the reported point and let the shared resolve-best block below
        % reconcile it against bestSoFar, which both optimizers maintain
        % identically from the GATED fitness.  JBest = NaN routes it there.
        try
            tBest = bestPoint(boResults);
            xBest = [tBest.dAR, tBest.maxdef];
        catch
            % No feasible point found: every design was gated or failed.
            warning('optimize_oblong_maxdef:noFeasibleBO', ...
                ['bayesopt found no point satisfying the constraints. ' ...
                 'Falling back to the best evaluated point.']);
            xBest = [defectAspectRatio_0, maxdef_0];
        end
        JBest = NaN;

        % Keep the study object with the run's other artifacts.
        try
            save([OPT.rootLoc, 'bayesopt_results.mat'], 'boResults');
        catch ME
            warning('optimize_oblong_maxdef:boSaveFailed', ...
                'Could not save the bayesopt results object: %s', ME.message);
        end

    otherwise
        error('optimize_oblong_maxdef:unknownOptimizer', ...
            ['Unknown OPT.optimizer "%s". Valid values are ''neldermead'' ' ...
             'and ''bayesopt''.'], OPT.optimizer);
end

%% ===================== RESOLVE BEST POINT =====================
best_defectAspectRatio = xBest(1);
best_maxdef            = xBest(2);

% Recompute the oblong mapping at the optimum.
mirrorAspectRatio = P.hy / P.hx;
r_best            = best_defectAspectRatio / mirrorAspectRatio;
best_oblong       = log(r_best) / (2*log(1 - best_maxdef));

% Undo the log scaling to report the maximized fitness on its natural scale.
% A J at or above the penalty floor means fminsearch finished on an
% out-of-bounds, gated, or failed vertex, which carries no fitness at all.
% The bayesopt branch sets JBest = NaN deliberately, which lands in the same
% place: the incumbent is then resolved from bestSoFar below, so both
% optimizers report the best point by exactly the same rule.
if isnan(JBest) || JBest >= OPT.Jpenalty
    fitnessBest = NaN;
else
    fitnessBest = 10^(-JBest);
end

% If fminsearch returned a penalized/failed point but a valid point was seen
% earlier, prefer the best valid point that was actually evaluated.
if (~isfinite(fitnessBest) || fitnessBest <= 0) && isfinite(bestSoFar.fitness)
    fprintf(['fminsearch endpoint was penalized/invalid; ' ...
             'falling back to best valid eval seen.\n']);
    best_defectAspectRatio = bestSoFar.defectAspectRatio;
    best_maxdef            = bestSoFar.maxdef;
    best_oblong            = bestSoFar.oblong;
    fitnessBest            = bestSoFar.fitness;
end

fprintf('\n==================== OPTIMIZATION COMPLETE ====================\n');
fprintf('Best defectAspectRatio = %.6f\n', best_defectAspectRatio);
fprintf('Best oblong            = %.4f\n', best_oblong);
fprintf('Best maxdef            = %.4f\n', best_maxdef);
fprintf('Best fitness           = %.4g\n', fitnessBest);
fprintf('Total evaluations      = %d\n', evalCount);
fprintf('Rejected unfabricable  = %d (no solve attempted)\n', nInfeasible);
if OPT.useStore
    fprintf('Served from store      = %d (COMSOL solves skipped)\n', ...
        nCacheHits);
    fprintf('Run store              -> %s [%s]\n', ...
        store.dbPath, store.backend);
end
fprintf('Log written            -> %s\n', logFile);

%% ===================== CONVERGENCE PLOTS =====================
% Nelder-Mead produces no landscape map the way the grid sweep does, so the
% only way to judge whether the run actually converged -- as opposed to
% exhausting MaxFunEvals, or stalling on a flat region of the Q_mech gate --
% is to look at the trace of evaluated points.
if isempty(history.eval)
    warning('optimize_oblong_maxdef:noEvals', ...
        'No in-bounds evaluations were made; skipping convergence plots.');
else
    % Gated (fit == 0) and failed (fit == NaN) points cannot go on a log axis.
    validFit = history.fitness;
    validFit(~isfinite(validFit) | validFit <= 0) = NaN;

    runningBest = nan(size(validFit));
    bestYet     = -Inf;
    for k = 1:numel(validFit)
        if isfinite(validFit(k)) && validFit(k) > bestYet
            bestYet = validFit(k);
        end
        if isfinite(bestYet); runningBest(k) = bestYet; end
    end

    figConv = figure('Name', [OPT.optimizer, ' convergence'], 'Color', 'w', ...
                     'Position', [100, 100, 1000, 420]);
    tiledlayout(figConv, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % --- panel 1: fitness vs evaluation number ---
    nexttile;
    semilogy(history.eval, validFit, 'o', 'MarkerSize', 5, ...
        'MarkerFaceColor', [0.60 0.75 0.95], 'MarkerEdgeColor', [0.20 0.40 0.70]);
    hold on;
    semilogy(history.eval, runningBest, '-', 'LineWidth', 1.8, ...
        'Color', [0.85 0.20 0.20]);
    rejected = history.eval(~(isfinite(history.fitness) & history.fitness > 0));
    if ~isempty(rejected) && any(isfinite(validFit))
        % Park the rejected evals on the bottom of the axis so a run that is
        % mostly gated is immediately obvious rather than showing as gaps.
        yFloor = min(validFit(isfinite(validFit))) * ones(size(rejected));
        semilogy(rejected, yFloor, 'x', 'Color', [0.5 0.5 0.5], 'MarkerSize', 7);
        legend({'evaluated', 'running best', 'gated / failed'}, ...
            'Location', 'southeast');
    else
        legend({'evaluated', 'running best'}, 'Location', 'southeast');
    end
    grid on;
    xlabel('evaluation number');
    ylabel('fitness');
    title('Fitness vs. evaluation');

    % --- panel 2: where the search went, inside the parameter box ---
    % Nelder-Mead walks a connected path, so joining the points is
    % meaningful.  bayesopt jumps wherever the acquisition function points,
    % so a connecting line would imply a trajectory that does not exist.
    isSimplex = strcmpi(OPT.optimizer, 'neldermead');
    nexttile;
    if isSimplex
        plot(history.dAR, history.maxdef, '-', 'Color', [0.7 0.7 0.7]);
        hold on;
    end
    scatter(history.dAR, history.maxdef, 34, history.eval, 'filled');
    hold on;
    cb = colorbar;
    cb.Label.String = 'evaluation number';
    plot(defectAspectRatio_0, maxdef_0, 'ks', 'MarkerSize', 11, 'LineWidth', 1.4);
    plot(best_defectAspectRatio, best_maxdef, 'r*', 'MarkerSize', 15, 'LineWidth', 1.6);
    plot([OPT.defectAspectRatio_min, OPT.defectAspectRatio_max, ...
          OPT.defectAspectRatio_max, OPT.defectAspectRatio_min, ...
          OPT.defectAspectRatio_min], ...
         [OPT.maxdef_min, OPT.maxdef_min, OPT.maxdef_max, ...
          OPT.maxdef_max, OPT.maxdef_min], 'k--', 'LineWidth', 1);
    grid on;
    xlabel('defectAspectRatio');
    ylabel('maxdef');
    if isSimplex
        title('Simplex trajectory');
        legend({'path', 'evals', 'start', 'best', 'bounds'}, 'Location', 'best');
    else
        title('Sampled designs');
        legend({'evals', 'start', 'best', 'bounds'}, 'Location', 'best');
    end

    saveas(figConv,  [OPT.rootLoc, 'optimize_convergence.png']);
    savefig(figConv, [OPT.rootLoc, 'optimize_convergence.fig']);
end

%% ===================== FINAL RE-RUN WITH PLOTS =====================
Pfinal = P;
Pfinal.oblong     = best_oblong;
Pfinal.maxdef     = best_maxdef;
Pfinal.plotgeom   = 1;
Pfinal.plotMech   = 1 * P.solveMech;
Pfinal.plotOpt    = 1 * P.solveOpt;
Pfinal.plotStrCpl = 1 * P.calcS;
Pfinal.saveplots  = 1;
if isfield(Pfinal, 'fileBase'); Pfinal = rmfield(Pfinal, 'fileBase'); end

bestLoc = [OPT.rootLoc, 'BEST', filesep];
fprintf('\nRe-running best design with plots -> %s\n', bestLoc);
[dsBest, ~] = RunNanobeamFEM(Pfinal, bestLoc);

%% ===================== CLOSE OUT THE RUN =====================
if OPT.useStore
    % Disarm the interrupted-marker BEFORE writing 'complete': clearing an
    % onCleanup fires its callback, so the other order would let the guard
    % overwrite the completion status on the way out.
    clear runGuard;
    store.finishRun(OPT.runId, 'complete');
    fprintf(['\nRun "%s" marked complete. To extend it, raise ' ...
             'MaxFunEvals and re-run\nwith the same OPT.runId -- the %d ' ...
             'stored point(s) replay for free.\n'], ...
        OPT.runId, store.evalCount(OPT.runId));
end

%% ===================== NESTED OBJECTIVE =====================
    function J = safeObjective(x)
    % Wrapper around objective() so that nothing can abort the optimization.
    % An uncaught error is scored the same as an unusable evaluation
    % (2 * OPT.Jpenalty) instead of a hard-coded magic number, keeping every
    % penalty level consistent with the log-scaled objective below.
        try
            J = objective(x);
        catch ME
            warning('optimize_oblong_maxdef:objError', ...
                'objective failed at x=[%.4f, %.4f]: %s', x(1), x(2), ME.message);
            J = 2 * OPT.Jpenalty;
        end
    end

    function [obj, con] = bayesObjective(t)
    % bayesopt adapter.  bayesopt passes a one-row TABLE, and wants the
    % objective plus a vector of coupled-constraint values (<= 0 means
    % satisfied).
    %
    % The Q_mech gate is expressed as a CONSTRAINT here rather than folded
    % into the objective.  Feeding the Nelder-Mead penalty ladder (1000,
    % 2000) to a GP that also sees feasible values around -16 would wreck the
    % kernel's length scale: the surrogate would spend its capacity modelling
    % the penalty cliff instead of the physics.  With a constraint, the GP
    % instead learns WHERE the gate boundary is, which is the actual design
    % question.
    %   con(1) : Q_mech gate,  1 - Qmech/QmechMin  (<= 0 when Qmech is high enough)
    %   con(2) : usable solve, -1 on success, +1 when there is no usable mode
    % A NaN objective marks an evaluation bayesopt could not use; it models
    % those with a separate error classifier rather than treating NaN as a
    % value.
        R = evalPoint([t.dAR, t.maxdef]);

        if isfinite(R.fitRaw) && R.fitRaw > 0
            obj = -log10(R.fitRaw);
        else
            obj = NaN;
        end

        if isfinite(R.Qmech)
            con1 = 1 - R.Qmech/OPT.QmechMin;
        else
            con1 = 1;   % no Q_mech at all counts as violating the gate
        end
        con2 = -1;
        if ~(isfinite(R.fitRaw) && R.fitRaw > 0)
            con2 = 1;
        end
        con = [con1, con2];
    end

    function J = objective(x)
    % Nelder-Mead adapter: sloped box barrier, then the penalty ladder.
        defectAspectRatio_i = x(1);
        maxdef_i            = x(2);

        % --- box constraints: sloped barrier rather than a flat plateau ---
        % A constant out-of-bounds penalty gives the simplex no direction to
        % move in once a vertex leaves the box, so it can collapse against a
        % bound and report false convergence.  Scaling the penalty by how far
        % outside the box the point sits -- normalized by each parameter's own
        % range so the two axes are comparable despite their different scales
        % -- restores a downhill direction back into the feasible set.
        % No FEM call is made for an infeasible point.
        dARrange = OPT.defectAspectRatio_max - OPT.defectAspectRatio_min;
        mdRange  = OPT.maxdef_max - OPT.maxdef_min;
        viol = max([(OPT.defectAspectRatio_min - defectAspectRatio_i) / dARrange, ...
                    (defectAspectRatio_i - OPT.defectAspectRatio_max) / dARrange, ...
                    (OPT.maxdef_min - maxdef_i) / mdRange, ...
                    (maxdef_i - OPT.maxdef_max) / mdRange, ...
                    0]);
        if viol > 0
            % Worse than any evaluated point (gated = 1x, unusable = 2x), and
            % strictly increasing with distance outside the box.
            J = OPT.Jpenalty * (2 + viol);
            return;
        end

        R = evalPoint(x);
        J = R.J;
    end

    function R = evalPoint(x)
    % Shared evaluation core for BOTH optimizers: maps x -> (oblong, maxdef),
    % serves the point from the store if it was already solved, otherwise
    % runs RunNanobeamFEM, then logs to the CSV, the store and the plotting
    % history.  Returns every score the callers might want:
    %   R.fit    - GATED fitness (0 below OPT.QmechMin) -- Nelder-Mead
    %   R.fitRaw - UNGATED fitness                      -- bayesopt
    %   R.Qmech  - mechanical radiative Q               -- the BO constraint
    %   R.J      - Nelder-Mead penalty-ladder objective
    % Keeping the expensive part in one place is what lets the two optimizers
    % share a run database and resume from each other's evaluations.
        defectAspectRatio_i = x(1);
        maxdef_i            = x(2);
        R = struct('fit', NaN, 'fitRaw', NaN, 'Qmech', NaN, ...
                   'J', 2*OPT.Jpenalty, 'cached', false, ...
                   'fabricable', true);

        % --- lithography pre-filter: cheapest rejection there is ---------
        % Same guards CreateNanobeamGeom throws on, but as arithmetic.
        % Without this the FEM pipeline starts building geometry and only
        % then errors, so an unfabricable candidate costs a real fraction of
        % a solve.  bayesopt also has this as its XConstraintFcn, so under
        % 'bayesopt' this branch is a backstop rather than the main line.
        [fabOK, fabMargin, fabWhy] = isFabricable( ...
            defectAspectRatio_i, maxdef_i, P);
        if ~fabOK
            nInfeasible = nInfeasible + 1;
            % Graded like the out-of-bounds barrier, and on the same tier:
            % both are points that are never solved at all.  Normalizing by
            % the minimum hole gap keeps the violation O(1).
            viol   = min(10, -fabMargin / 50e-9);
            R.J    = OPT.Jpenalty * (2 + viol);
            R.fabricable = false;
            fprintf(['  eval  --: dAR=%.4f  maxdef=%.4f  ' ...
                     '[UNFABRICABLE: %s, %.1f nm short]\n'], ...
                defectAspectRatio_i, maxdef_i, fabWhy{1}, -fabMargin*1e9);
            return;
        end

        % --- resume fast-forward: was this exact point already solved? ---
        % fminsearch is deterministic, so a resumed run retraces its earlier
        % path exactly before extending it.  Serving those repeats from the
        % database turns the replay into seconds instead of re-solving hours
        % of COMSOL.  Matching is on OptimStore.paramKeyOf, which rounds to
        % 10 significant digits so a value that has round-tripped through
        % SQLite's REAL storage still matches.
        if OPT.useStore && OPT.resume
            hit = store.lookup(OPT.runId, [defectAspectRatio_i, maxdef_i]);
            if ~isempty(hit)
                % Re-derive rather than trusting the stored objective_j: the
                % row may have been written by the OTHER optimizer, whose
                % scoring rule differs.  The physics (fitness, fitnessRaw,
                % Qmech) is what is actually shared.
                R.fit    = hit.fitness;
                R.fitRaw = hit.fitnessRaw;
                R.Qmech  = hit.Qmech;
                R.cached = true;
                if isfinite(R.fit) && R.fit > 0
                    R.J = -log10(R.fit);
                elseif isfinite(R.fit)
                    R.J = OPT.Jpenalty;
                else
                    R.J = 2*OPT.Jpenalty;
                end
                J = R.J;
                nCacheHits = nCacheHits + 1;
                fprintf(['  eval %3d: oblong=%.4f  maxdef=%.4f  ' ...
                         '[CACHED  fitness=%.4g  %s]\n'], ...
                    hit.eval_idx, hit.oblong, maxdef_i, hit.fitness, ...
                    hit.status);

                if isfinite(hit.fitness) && hit.fitness > bestSoFar.fitness
                    bestSoFar.fitness = hit.fitness;
                    bestSoFar.oblong  = hit.oblong;
                    bestSoFar.maxdef  = maxdef_i;
                    bestSoFar.defectAspectRatio = defectAspectRatio_i;
                end

                % Replayed points join the history so the convergence plots
                % show the whole search, not just this session's part of it.
                history.eval(end+1, 1)    = hit.eval_idx;
                history.dAR(end+1, 1)     = defectAspectRatio_i;
                history.maxdef(end+1, 1)  = maxdef_i;
                history.oblong(end+1, 1)  = hit.oblong;
                history.fitness(end+1, 1) = hit.fitness;
                history.Qmech(end+1, 1)   = NaN;
                history.J(end+1, 1)       = J;

                % Keep the counter ahead of every index already stored, so a
                % genuinely new point cannot collide with a stored one on
                % the (run_id, eval_idx) primary key.
                evalCount = max(evalCount, hit.eval_idx);
                return;
            end
        end

        evalCount = evalCount + 1;

        Pe = P;

        % --- map defect aspect ratio to oblong (same as the sweep) ---
        mirrorAspectRatio = P.hy / P.hx;
        r = defectAspectRatio_i / mirrorAspectRatio;

        Pe.maxdef = maxdef_i;
        Pe.oblong = log(r) / (2*log(1 - maxdef_i));

        % Unique folder per evaluation — RunNanobeamFEM never reuses a cache.
        evLoc = sprintf('%seval_%04d_ob%.4f_md%.4f%s', ...
            OPT.rootLoc, evalCount, Pe.oblong, maxdef_i, filesep);

        % Strip any stale fileBase so CreateFileBase regenerates it.
        if isfield(Pe, 'fileBase'); Pe = rmfield(Pe, 'fileBase'); end

        fprintf('  eval %3d: oblong=%.4f  maxdef=%.4f  -> %s\n', ...
            evalCount, Pe.oblong, maxdef_i, evLoc);

        wM_i     = NaN;
        gOM_i    = NaN;
        gMB_i    = NaN;
        gPE_i    = NaN;
        LSiV_i   = NaN;
        Qmech_i  = NaN;
        Qopt_i   = NaN;
        lambda_i = NaN;
        fit      = NaN;
        fitRaw   = NaN;
        status_i = 'failed';

        try
            [ds, ~] = RunNanobeamFEM(Pe, evLoc);

            % Treat a missing or empty ds as a silent failure.
            if isempty(ds) || ~isstruct(ds)
                error('optimize_oblong_maxdef:emptyDs', ...
                    'RunNanobeamFEM returned empty ds');
            end

            % --- check for localised mechanical modes ---
            hasMech = false;
            if isfield(ds, 'mfem') && isfield(ds.mfem, 'freqs') && ...
                    ~isempty(ds.mfem.freqs)
                hasMech = true;
            end
            if ~hasMech
                warning('optimize_oblong_maxdef:noMechModes', ...
                    'eval %d (oblong=%.4f, maxdef=%.4f): no localised mechanical modes found', ...
                    evalCount, Pe.oblong, maxdef_i);
                status_i = 'no_mech_modes';
            else
                status_i = 'ok';
            end

            % --- extract mechanical frequency (prefer coupled-mode result) ---
            if isfield(ds, 'cpl') && isfield(ds.cpl, 'mechFreq') && ...
                    ~isempty(ds.cpl.mechFreq)
                wM_i = max(real(ds.cpl.mechFreq));
            elseif hasMech
                wM_i = max(real(ds.mfem.freqs));
            end

            % --- check for high-Q optical modes (when solveOpt = 1) ---
            if Pe.solveOpt
                hasOpt = isfield(ds, 'ofem') && isfield(ds.ofem, 'freqs') && ...
                         ~isempty(ds.ofem.freqs);
                if ~hasOpt
                    warning('optimize_oblong_maxdef:noOptModes', ...
                        'eval %d (oblong=%.4f, maxdef=%.4f): no high-Q optical modes found', ...
                        evalCount, Pe.oblong, maxdef_i);
                    if strcmp(status_i, 'ok')
                        status_i = 'no_opt_modes';
                    else
                        status_i = [status_i, '+no_opt_modes'];
                    end
                end
            end

            % --- extract gOM and its mechanism breakdown ---
            % gMB (moving boundary) and gPE (photoelastic) are recorded by
            % CalcGOM at the SAME mode pair as the reported maximum, so the
            % three numbers always correspond to one another.  They carry
            % opposite signs and routinely partially cancel, which is why
            % the net alone does not tell you whether a design is robust.
            [gOM_i, gMB_i, gPE_i] = local_getGOM(ds);

            % --- extract LSiV ---
            if isfield(ds, 'cpl') && isfield(ds.cpl, 'SiV') && ...
                    isfield(ds.cpl.SiV, 'LSiVMaxMode') && ...
                    ~isempty(ds.cpl.SiV.LSiVMaxMode)
                LSiV_i = max(real(ds.cpl.SiV.LSiVMaxMode));
            end

            % --- extract mechanical radiative Q (requires solveMechPML) ---
            [Qmech_i, ~] = local_getQmech(ds);

            % --- extract optical Q and wavelength (requires solveOpt) ---
            [Qopt_i, lambda_i] = local_getQopt(ds, OPT.targetLambda);

            % Fitness is NaN when modes are absent.  Both the gated and the
            % ungated score are computed: Nelder-Mead uses the first, the GP
            % surrogate the second (see OPT.fitnessFcnRaw).
            if hasMech
                fit = OPT.fitnessFcn(ds);
                if isempty(fit) || ~isfinite(fit)
                    fit = NaN;
                end
                fitRaw = OPT.fitnessFcnRaw(ds);
                if isempty(fitRaw) || ~isfinite(fitRaw)
                    fitRaw = NaN;
                end
            end

        catch ME
            warning('optimize_oblong_maxdef:evalFailed', ...
                'eval %d (oblong=%.4f, maxdef=%.4f) failed: %s', ...
                evalCount, Pe.oblong, maxdef_i, ME.message);
            status_i = 'failed';
        end

        % --- flag designs rejected by the binary Q_mech gate ---
        % fit == 0 means the point solved cleanly but Q_mech fell below
        % OPT.QmechMin.  That is a very different outcome from 'failed', and
        % conflating the two makes the CSV impossible to interpret afterwards.
        if isfinite(fit) && fit == 0 && strcmp(status_i, 'ok')
            status_i = 'Qmech_gated';
        end

        fprintf(['         wM=%.3f GHz   gOM=%.2f kHz (MB %.2f + PE %.2f)   ' ...
                 'LSiV=%.4f MHz\n         Qm=%.3g   Qo=%.3g   lam=%.1f nm   ' ...
                 'fitness=%.4g   [%s]\n'], ...
            wM_i*1e-9, gOM_i*1e-3, gMB_i*1e-3, gPE_i*1e-3, LSiV_i*1e-6, ...
            Qmech_i, Qopt_i, lambda_i, fit, status_i);

        % Append to CSV immediately so the log survives a mid-run crash.
        fidLog = fopen(logFile, 'a');
        fprintf(fidLog, ['%d,%.6f,%.6f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,', ...
                         '%.6g,%.6g,%.6g,%s\n'], ...
            evalCount, Pe.oblong, maxdef_i, ...
            wM_i*1e-9, gOM_i*1e-3, gMB_i*1e-3, gPE_i*1e-3, LSiV_i*1e-6, ...
            Qmech_i, Qopt_i, lambda_i, fit, status_i);
        fclose(fidLog);

        % --- record the search trajectory for the convergence plots ---
        history.eval(end+1, 1)    = evalCount;
        history.dAR(end+1, 1)     = defectAspectRatio_i;
        history.maxdef(end+1, 1)  = maxdef_i;
        history.oblong(end+1, 1)  = Pe.oblong;
        history.fitness(end+1, 1) = fit;
        history.Qmech(end+1, 1)   = Qmech_i;

        % --- form the negated, log-scaled objective ---
        % Raw fitness spans ~1e10-1e18 (g_OM^2 * Q_opt), and optimset's TolFun
        % is an ABSOLUTE tolerance on the spread of objective values.  With
        % J = -fit, TolFun = 1e-6 can never be met, so the run always burns
        % through MaxFunEvals regardless of convergence.  Minimizing
        % -log10(fitness) puts J on an O(10) scale, which turns TolFun into an
        % effectively relative test on the fitness and lets the run stop early.
        %   fit  > 0   -> usable design
        %   fit == 0   -> solved, but rejected by the binary Q_mech gate
        %   fit  = NaN -> no usable mode / solve failed (not evaluable at all)
        % Both infeasible cases are penalized, the un-evaluable one more
        % heavily, so the simplex prefers drifting back through
        % gated-but-solvable territory rather than through the dead zone.
        if isfinite(fit) && fit > 0
            J = -log10(fit);
            if fit > bestSoFar.fitness
                bestSoFar.fitness = fit;
                bestSoFar.oblong  = Pe.oblong;
                bestSoFar.maxdef  = maxdef_i;
                bestSoFar.defectAspectRatio = defectAspectRatio_i;
            end
        elseif isfinite(fit)
            J = OPT.Jpenalty;        % fit == 0: rejected by the Q_mech gate
        else
            J = 2 * OPT.Jpenalty;    % fit == NaN: no usable result
        end

        % --- persist immediately so an interrupted run resumes from here ---
        % Written after J exists but before returning, so a stored row is
        % always a COMPLETE record: a run killed between the FEM solve and
        % this line simply re-solves that one point next time.
        if OPT.useStore
            try
                store.logEval(OPT.runId, evalCount, struct( ...
                    'x',        [defectAspectRatio_i, maxdef_i], ...
                    'oblong',   Pe.oblong, ...
                    'wM',       wM_i,     'gOM',      gOM_i, ...
                    'gMB',      gMB_i,    'gPE',      gPE_i, ...
                    'LSiV',     LSiV_i,   'Qmech',    Qmech_i, ...
                    'Qopt',     Qopt_i,   'lambdaNm', lambda_i, ...
                    'fitness',  fit,      'fitnessRaw', fitRaw, ...
                    'J',        J,        'status',     status_i, ...
                    'evLoc',    evLoc));
            catch ME
                % Persistence must never kill a study that is already
                % several COMSOL-hours deep.
                warning('optimize_oblong_maxdef:storeFailed', ...
                    'Could not persist eval %d: %s', evalCount, ME.message);
            end
        end

        history.J(end+1, 1) = J;

        R.fit    = fit;
        R.fitRaw = fitRaw;
        R.Qmech  = Qmech_i;
        R.J      = J;
    end

end  % function optimize_oblong_maxdef

%% ===================== LOCAL FITNESS HELPERS =====================
%
% ACTIVE fitness  : fitness_coopProduct  (wired up via OPT.fitnessFcn above).
%
% The other fitness_* functions in this section are NOT dead code left behind
% by accident -- they are ready-made ALTERNATIVE objectives.  To switch, edit
% the single OPT.fitnessFcn line in the TUNABLE KNOBS section:
%
%   OPT.fitnessFcn = @(ds) fitness_QxQxLambda(ds, OPT.targetLambda, OPT.lambdaTol);
%   OPT.fitnessFcn = @(ds) fitness_maxGOM(ds);
%   OPT.fitnessFcn = @(ds) fitness_maxLSiV(ds);
%
% Anything returning NaN / non-finite is treated as an unusable point and
% penalized in the objective; anything returning 0 is treated as a rejected
% (but successfully evaluated) design.

function fit = fitness_coopProduct(ds, targetLambdaNm, tolNm, QmechMin)
%FITNESS_COOPPRODUCT  Cooperativity-like figure of merit.  [ACTIVE]
%
%   Q_mech enters as a BINARY acceptance gate, not as a multiplicative factor:
%
%       Q_mech <  QmechMin  ->  fit = 0    (design rejected outright)
%       Q_mech >= QmechMin  ->  fit = g_OM^2 * Q_opt * L(lambda_opt)
%
%   where L is a GAUSSIAN wavelength-matching factor
%
%       L(lambda) = exp( -((lambda - targetLambda)/tolNm)^2 )
%
%   so tolNm is the 1/e half-width: a tolNm miss costs a factor 1/e, a
%   2*tolNm miss costs 1/e^4.
%
%   Graceful degradation:
%     no Q_mech at all (PML off / mech solve failed) -> NaN (unusable point)
%     no g_OM (P.calcG = 0)                          -> g^2 factor dropped
%     optics off / no optical mode                   -> g_OM^2 alone
%
%   NOTE ON THE BINARY GATE: it makes the objective piecewise-flat -- every
%   design below QmechMin scores exactly 0 and is indistinguishable from every
%   other, so the simplex gets no direction to move in while it sits inside a
%   rejected region.  If the run stalls immediately, either lower
%   OPT.QmechMin or switch to fitness_QxQxLambda, which varies smoothly with
%   Q_mech and therefore always supplies a gradient.
fit = NaN;
[Qmech, ~]        = local_getQmech(ds);
[Qopt,  lambdaNm] = local_getQopt(ds, targetLambdaNm);

% No Q_mech at all: the point could not be evaluated -- that is different from
% the point being bad, so return NaN rather than 0.
if isnan(Qmech); return; end

% --- binary Q_mech acceptance gate ---
if Qmech < QmechMin
    fit = 0;
    return;
end

gOM     = local_getGOM(ds);
gFactor = 1;
if ~isnan(gOM); gFactor = gOM^2; end

if isnan(Qopt) || isnan(lambdaNm)
    % Optics unavailable: keep only the factors we actually have.  Scores from
    % this branch are NOT on the same scale as full-branch scores, so a run
    % that mixes the two cannot be ranked meaningfully.
    fit = gFactor;
else
    fit = gFactor * Qopt * exp(-((lambdaNm - targetLambdaNm) / tolNm)^2);
end
end

function [Qmech, fMech] = local_getQmech(ds)
% Mechanical radiative Q of the best localised cavity mode (solveMechPML=1).
% Selects the localised mode with the highest localization ratio; falls back
% to the highest-Q solved mode if localisation filter was not applied.
Qmech = NaN;  fMech = NaN;
if ~isfield(ds, 'mfem'); return; end
mfem = ds.mfem;
if ~isfield(mfem, 'QAll') || isempty(mfem.QAll); return; end
if isfield(mfem, 'locInd') && ~isempty(mfem.locInd) && ...
        isfield(mfem, 'locRatio') && ~isempty(mfem.locRatio) && ...
        max(mfem.locInd) <= numel(mfem.locRatio)
    [~, k] = max(mfem.locRatio(mfem.locInd));
    idx    = mfem.locInd(k);
else
    [~, idx] = max(mfem.QAll);
end
if idx <= numel(mfem.QAll)
    Qmech = mfem.QAll(idx);
end
if isfield(mfem, 'freqs') && idx <= numel(mfem.freqs)
    fMech = mfem.freqs(idx);
end
end

function [Qopt, lambdaNm] = local_getQopt(ds, targetLambdaNm)
% Optical Q and wavelength (nm) of the high-Q mode nearest targetLambdaNm.
% Uses the high-Q subset ds.ofem.Q / ds.ofem.lambda (Q>1e4); falls back to
% the full ds.ofem.QAll / ds.ofem.lambdaAll set if the subset is empty.
% ds.ofem.lambda is in meters; returned lambdaNm is in nm.
Qopt = NaN;  lambdaNm = NaN;
if ~isfield(ds, 'ofem'); return; end
ofem = ds.ofem;
if isfield(ofem, 'lambda') && ~isempty(ofem.lambda) && ...
        isfield(ofem, 'Q')   && ~isempty(ofem.Q)
    lam = ofem.lambda(:) * 1e9;  Q = ofem.Q(:);
elseif isfield(ofem, 'lambdaAll') && ~isempty(ofem.lambdaAll) && ...
        isfield(ofem, 'QAll')      && ~isempty(ofem.QAll)
    lam = ofem.lambdaAll(:) * 1e9;  Q = ofem.QAll(:);
else
    return;
end
[~, k]   = min(abs(lam - targetLambdaNm));
Qopt     = Q(k);
lambdaNm = lam(k);
end

function [gOM, gMB, gPE] = local_getGOM(ds)
%LOCAL_GETGOM  Net optomechanical coupling and its mechanism breakdown.
%
%   gOM - max |g_OM| in Hz.  NaN if not computed (P.calcG = 0).
%   gMB - moving-boundary part, Hz, SIGNED.
%   gPE - photoelastic part, Hz, SIGNED.
%
%   gMB and gPE come from cpl.gMBmax / cpl.gPEmax, which CalcGOM records at
%   the same (oSol, mSol) pair as the reported maximum, so gMB + gPE is the
%   signed net whose magnitude is gOM.  They are returned SIGNED on purpose:
%   the two mechanisms carry opposite signs and routinely partially cancel,
%   and taking abs() of each would hide exactly that.  A design sitting in a
%   cancellation notch has a small gOM despite large individual parts, and
%   is fragile to fabrication error -- see cpl.breakdown for the full table.
gOM = NaN;  gMB = NaN;  gPE = NaN;
if ~isfield(ds, 'cpl'); return; end
c = ds.cpl;
if isfield(c, 'gMax') && ~isempty(c.gMax)
    gOM = max(abs(c.gMax));
end
if isfield(c, 'gMBmax') && ~isempty(c.gMBmax)
    gMB = c.gMBmax(1);
end
if isfield(c, 'gPEmax') && ~isempty(c.gPEmax)
    gPE = c.gPEmax(1);
end
end

function fit = fitness_QxQxLambda(ds, targetLambdaNm, tolNm)
%FITNESS_QXQXLAMBDA  Q_opt * Q_mech, exponential wavelength gate.  [ALTERNATIVE]
%
%   fitness = Q_opt * Q_mech * exp(-|lambda_opt - targetLambda| / tol)
%
%   Differs from the active fitness_coopProduct in two ways, both deliberate:
%     1. Q_mech is a smooth multiplicative factor, not a binary gate -- use
%        this one if the binary gate flattens the landscape and stalls the
%        simplex, or if you want to see how Q_mech trades against Q_opt.
%     2. The wavelength penalty decays EXPONENTIALLY (Laplacian) in |dlambda|,
%        not as a Gaussian, so it punishes small detunings harder and large
%        detunings more gently than the active fitness.
%     3. g_OM does not appear at all -- use fitness_coopProduct when
%        optomechanical coupling is what you actually care about.
%
%   Graceful degradation:
%     no Q_mech (PML off / absent)    -> NaN (unusable point)
%     optics off / no optical mode    -> Q_mech alone
%
%   To use:  OPT.fitnessFcn = @(ds) fitness_QxQxLambda(ds, OPT.targetLambda, OPT.lambdaTol);
fit = NaN;
[Qmech, ~]        = local_getQmech(ds);
[Qopt,  lambdaNm] = local_getQopt(ds, targetLambdaNm);
if isnan(Qmech); return; end
if isnan(Qopt) || isnan(lambdaNm)
    fit = Qmech;
else
    fit = Qopt * Qmech * exp(-abs(lambdaNm - targetLambdaNm) / tolNm);
end
end

function fit = fitness_maxGOM(ds)
%FITNESS_MAXGOM  Maximize bare optomechanical coupling.  [ALTERNATIVE]
%
%   Returns max |g_OM| in Hz.  NaN if gOM not computed (P.calcG = 0).
%   Ignores both Q's and the optical wavelength entirely -- use it only for a
%   quick "where is the coupling strongest" scan, since nothing stops the
%   optimizer from parking on a design with a useless optical resonance.
%
%   To use:  OPT.fitnessFcn = @(ds) fitness_maxGOM(ds);
fit = NaN;
if isfield(ds, 'cpl') && isfield(ds.cpl, 'gMax') && ~isempty(ds.cpl.gMax)
    fit = max(abs(ds.cpl.gMax));
end
end

function fit = fitness_maxLSiV(ds)
%FITNESS_MAXLSIV  Maximize SiV strain coupling.  [ALTERNATIVE]
%
%   Returns max lambda_SiV in Hz.  NaN if strain coupling was not computed --
%   requires P.calcS = 1, which is currently 0 in the knobs above, so remember
%   to enable it before selecting this objective.
%
%   To use:  P.calcS = 1;  OPT.fitnessFcn = @(ds) fitness_maxLSiV(ds);
fit = NaN;
if isfield(ds, 'cpl') && isfield(ds.cpl, 'SiV') && ...
        isfield(ds.cpl.SiV, 'LSiVMaxMode') && ~isempty(ds.cpl.SiV.LSiVMaxMode)
    fit = max(real(ds.cpl.SiV.LSiVMaxMode));
end
end
