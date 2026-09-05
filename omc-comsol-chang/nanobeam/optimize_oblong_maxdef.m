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
hx_0 = 196e-9;      % nominal mirror hole height
hy_0 = 578e-9;      % nominal mirror hole width

% --- starting point (same as sweep centre) ---
% hy ~ (1-maxdef)^(1+oblong), hx ~ (1-maxdef)^(1-oblong).
% oblong=0 -> equal scaling; oblong=1 -> hx constant, only hy shrinks.
oblong_0 = 1.96;
maxdef_0 = 0.22;
defectAspectRatio_0 = (hy_0/hx_0) * (1-maxdef_0)^(2*oblong_0);  % hy/hx * (1-maxdef)^(2*oblong)

% --- defectAspectRatio box constraints (dimensionless) ---
OPT.defectAspectRatio_min = defectAspectRatio_0*0.5;
OPT.defectAspectRatio_max = defectAspectRatio_0*1.5;

% --- maxdef box constraints (dimensionless, 0–<1) ---
% Fraction by which the mirror lattice constant is reduced at cavity centre.
OPT.maxdef_min = maxdef_0*0.5;
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
% 'neldermead' : fminsearch simplex (base MATLAB, works today).
% 'bayesopt'   : Gaussian-process Bayesian optimization -- reserved; see
%                bayesopt_oblong_maxdef_PROPOSAL.md.  Needs the Statistics
%                and Machine Learning Toolbox, which is licensed but NOT
%                installed on this machine.
OPT.optimizer = 'neldermead';

% --- run persistence / stop-and-resume (see OptimStore.m) -----------------
% Every evaluation is written to SQLite as it completes, so an interrupted
% study never loses a COMSOL solve.  Re-running with the SAME OPT.runId and
% OPT.resume = 1 replays the stored points for free, then carries on.
%   OPT.runId  : names the study.  Change it to start a fresh search, keep
%                it to continue one.  Deliberately NOT date-stamped -- a
%                date-stamped id would silently start a new study tomorrow.
%   OPT.resume : 1 = serve already-evaluated points from the database,
%                0 = re-solve everything (the database is still written).
OPT.useStore = 1;
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

%% ===================== P STRUCT DEFAULTS =====================
% Geometry matches test_nanobeamRectFEM_withPML.m (fabricated OMC device).
P.xsect     = 'isoFit';     % beam cross sectional shape - 'tri' or 'rect' or 'isoFit'
P.beamMat   = 'diamond';
P.celltype  = 'hole';
P.anisoMat  = 1;
P.rxtal     = 0;
P.rxtalInFilename = 1;

P.a   = 529e-9;     % nominal lattice constant (m)
P.w   = 750e-9;     % beam width (m)
P.theta = 45;       % etch angle (degrees; no effect for isoFit/rect)
P.th  = 500e-9;     % beam thickness (m)
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
    fprintf(fid, ['eval,oblong,maxdef,wM_GHz,gOM_kHz,LSiV_MHz,' ...
                  'Q_mech,Q_opt,lambda_opt_nm,fitness,status\n']);
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

%% ===================== RUN STORE (SQLite stop / resume) ==================
store = [];
if OPT.useStore
    store = OptimStore(OPT.dbPath);
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

% Wrap the objective in a try/catch so the optimizer never crashes.
objFcn = @(x) safeObjective(x);

switch lower(OPT.optimizer)
    case 'neldermead'
        [xBest, JBest] = fminsearch(objFcn, x0, options);

    case 'bayesopt'
        % Gaussian-process Bayesian optimization. NOT WIRED UP YET -- the
        % design is in bayesopt_oblong_maxdef_PROPOSAL.md, and
        % ../bayesopt_boomerang.m is the working reference implementation.
        % The check below is deliberately specific: on this machine the
        % Statistics and Machine Learning Toolbox is LICENSED but NOT
        % INSTALLED, which otherwise surfaces as a bare 'Unrecognized
        % function' that reads like a licence problem instead of an
        % install one, and sends you debugging the wrong thing.
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
        error('optimize_oblong_maxdef:bayesoptNotWired', ...
            ['OPT.optimizer = ''bayesopt'' is reserved but not implemented ' ...
             'yet.\nSee bayesopt_oblong_maxdef_PROPOSAL.md for the design; ' ...
             'set OPT.optimizer = ''neldermead'' to run today.']);

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
if JBest >= OPT.Jpenalty
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

    figConv = figure('Name', 'fminsearch convergence', 'Color', 'w', ...
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

    % --- panel 2: simplex trajectory inside the parameter box ---
    nexttile;
    plot(history.dAR, history.maxdef, '-', 'Color', [0.7 0.7 0.7]);
    hold on;
    scatter(history.dAR, history.maxdef, 34, history.eval, 'filled');
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
    title('Simplex trajectory');
    legend({'path', 'evals', 'start', 'best', 'bounds'}, 'Location', 'best');

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

    function J = objective(x)
    % Core objective: maps x -> (oblong, maxdef), runs RunNanobeamFEM in a
    % unique folder, evaluates the fitness, logs to CSV, and returns the
    % negated log10 fitness (see the objective-forming block at the bottom).
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
                J = hit.J;
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
        LSiV_i   = NaN;
        Qmech_i  = NaN;
        Qopt_i   = NaN;
        lambda_i = NaN;
        fit      = NaN;
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

            % --- extract gOM ---
            if isfield(ds, 'cpl') && isfield(ds.cpl, 'gMax') && ...
                    ~isempty(ds.cpl.gMax)
                gOM_i = max(abs(ds.cpl.gMax));
            end

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

            % Fitness is NaN when modes are absent.
            if hasMech
                fit = OPT.fitnessFcn(ds);
                if isempty(fit) || ~isfinite(fit)
                    fit = NaN;
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

        fprintf(['         wM=%.3f GHz   gOM=%.2f kHz   LSiV=%.4f MHz   ' ...
                 'Qm=%.3g   Qo=%.3g   lam=%.1f nm   fitness=%.4g   [%s]\n'], ...
            wM_i*1e-9, gOM_i*1e-3, LSiV_i*1e-6, Qmech_i, Qopt_i, lambda_i, fit, status_i);

        % Append to CSV immediately so the log survives a mid-run crash.
        fidLog = fopen(logFile, 'a');
        fprintf(fidLog, '%d,%.6f,%.6f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%s\n', ...
            evalCount, Pe.oblong, maxdef_i, ...
            wM_i*1e-9, gOM_i*1e-3, LSiV_i*1e-6, Qmech_i, Qopt_i, lambda_i, fit, status_i);
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
                    'LSiV',     LSiV_i,   'Qmech',    Qmech_i, ...
                    'Qopt',     Qopt_i,   'lambdaNm', lambda_i, ...
                    'fitness',  fit,      'J',        J, ...
                    'status',   status_i, 'evLoc',    evLoc));
            catch ME
                % Persistence must never kill a study that is already
                % several COMSOL-hours deep.
                warning('optimize_oblong_maxdef:storeFailed', ...
                    'Could not persist eval %d: %s', evalCount, ME.message);
            end
        end

        history.J(end+1, 1) = J;
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

function gOM = local_getGOM(ds)
% Max |g_OM| in Hz of the cavity mode, or NaN if not computed (P.calcG = 0).
gOM = NaN;
if isfield(ds, 'cpl') && isfield(ds.cpl, 'gMax') && ~isempty(ds.cpl.gMax)
    gOM = max(abs(ds.cpl.gMax));
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
