function [best_oblong, best_maxdef, fitnessBest, dsBest, OPT, stats] = optimize_rib_oblong_maxdef(optOverrides)
%OPTIMIZE_RIB_OBLONG_MAXDEF  Optimize the rib-cavity defect taper.
%
%   Rib counterpart of optimize_oblong_maxdef.m. Searches two parameters of
%   the rib cavity's hole-style taper (P.celltype = 'rib',
%   P.ribTaperScheme = 'hole'); each objective evaluation is a full COMSOL
%   mechanical + optical solve via RunNanobeamFEM:
%       x(1) = oblong   taper exponent P.oblong, searched DIRECTLY
%       x(2) = maxdef   maximum defect ratio, dimensionless
%   (Unlike optimize_oblong_maxdef.m, which searches the centre-cell aspect
%   ratio and maps it to oblong. The ratio is still derived and recorded.)
%
%   [best_oblong, best_maxdef, fitnessBest, dsBest, OPT, stats] = ...
%       optimize_rib_oblong_maxdef()
%   [...] = optimize_rib_oblong_maxdef(optOverrides)
%       optOverrides is merged over the TUNABLE KNOBS below, nested structs
%       (OPT.bo, OPT.nm) field by field, e.g.
%           optimize_rib_oblong_maxdef(struct('dryRun', true))
%           optimize_rib_oblong_maxdef(struct('optimizer', 'neldermead'))
%           optimize_rib_oblong_maxdef(struct('bo', struct('maxEvals', 90)))
%       An unknown field name is an error, never silently ignored.
%
%   See optimization_README.md ("Rib cavity counterpart") for the guide.
%
%   OPTIMIZER (OPT.optimizer)
%     'bayesopt'   - Gaussian-process Bayesian optimization (DEFAULT). Needs
%                    the Statistics and Machine Learning Toolbox.
%     'neldermead' - fminsearch simplex, base MATLAB.
%   Both write the same physics to the same run store, so a study started
%   with one can be resumed or seeded with the other.
%
%   SEARCH SPACE AND THE DERIVED ASPECT RATIO
%   The 'hole' taper scheme (createRibCavityGeom.m, ribLength.m) gives, with
%   sig_k = 1 - maxdef*f((k-1)/ndef):
%       a_k = a*sig_k,  t_k = t*sig_k^(1-oblong),  w_k = w_mirror*sig_k^(1+oblong)
%   x = [oblong, maxdef] goes into P.oblong / P.maxdef unchanged. The
%   centre-cell aspect ratio w/t (sig = 1 - maxdef) follows from it,
%       defectAspectRatio = (w_mirror/t) * (1-maxdef)^(2*oblong)
%   and is what isFabricableRib takes and what the run store is KEYED on
%   (its 'dar' column; 'oblong' is stored beside it). The store is therefore
%   independent of the search space: a study recorded by the earlier version
%   of this script, which searched defectAspectRatio, resumes here unchanged.
%   The buildable region is not a rectangle in (oblong, maxdef): the corners
%   of the box at large maxdef are unbuildable. They are never solved (see
%   FABRICATION PRE-FILTER) and the start-up scan reports their share.
%
%   FITNESS (OPT.fitnessFcn, default fitness_ribCoop)
%
%       fitness = g_OM^2 * Q_opt * exp(-((lambda_opt - targetLambda)/lambdaTol)^2)
%
%   This is the hole optimizer's fitness_coopProduct WITHOUT its binary
%   Q_mech gate, because the rib cavity has no mechanical PML:
%   buildRibCavity errors on P.solveMechPML = 1, and SolveNanobeamFEM only
%   fills ds.mfem.QAll under a PML. Nothing in the default objective
%   therefore rewards mechanical confinement; only SolveNanobeamFEM's
%   localization filter (P.extractLocMechModes) acts on the mechanics.
%   Two opt-in substitutes:
%     OPT.useQmechGate = 1  - gate on a Q_mech ESTIMATE taken from the
%         complex eigenvalue that the low-reflecting beam end
%         (P.mEndBC = 'lowreflecting') produces. A first-order absorbing
%         boundary, not a PML: use it to rank designs, not as a Q value.
%         The estimate is recorded in the Q_mech column either way.
%     OPT.mechFreqWindow = [fLo fHi]  - smooth roll-off of the fitness when
%         the mechanical mode of the best g_OM pair leaves the window (e.g.
%         the mirror's phononic band gap).
%
%   OBJECTIVE SCALING
%   Both optimizers minimize -log10(fitness). Under Nelder-Mead the
%   infeasible outcomes are penalized in increasing order of badness:
%       gated (Q_mech too low, gate on) ->  1 * OPT.Jpenalty
%       unusable (no mode / error)      ->  2 * OPT.Jpenalty
%       unfabricable / out of bounds    ->  OPT.Jpenalty * (2 + distance)
%   Under bayesopt they are coupled constraints plus an XConstraintFcn.
%
%   FABRICATION PRE-FILTER (isFabricableRib.m)
%   The builder's validity rules and the 50 nm lithography floors, as
%   vectorized arithmetic. No unbuildable candidate ever costs a solve.
%
%   STOP AND RESUME (OptimStore.m)
%   Every evaluation is persisted as it completes.
%     OPT.runId  : names the study. KEEP it to continue, CHANGE it to restart.
%     OPT.resume : 1 = serve stored points, 0 = re-solve (still records).
%     OPT.dbPath : the run store (separate from the hole optimizer's, because
%                  the stored 'dAR' means w/t here and hy/hx there).
%   Stored points are matched on (defectAspectRatio, maxdef), so they are
%   reused whichever variable the session that wrote them searched. Stored
%   points outside the current oblong box still seed the Bayesian surrogate.
%   Ctrl-C is safe: an onCleanup guard marks the run 'interrupted'.
%
%   DRY RUN (OPT.dryRun = 1)
%   Replaces RunNanobeamFEM by an analytic stand-in (dryRunDs) and writes to
%   tempdir instead of ./test, so the whole script -- pre-filter, both
%   optimizers, store, resume, plots -- runs in seconds with no COMSOL.
%   NOT physics. See test_optimize_rib_dryrun.m.
%
%   Returns:
%     best_oblong  - oblong at the optimum
%     best_maxdef  - maxdef at the optimum
%     fitnessBest  - fitness (positive, un-scaled) at the optimum
%     dsBest       - ds struct from the final re-run of the best point
%     OPT          - the settings actually used (after overrides)
%     stats        - counters: nEval (highest evaluation index), nSolved,
%                    nCacheHits, nInfeasible; plus stats.bestCpl, the
%                    coupling of the reported best design in Hz:
%                      .gOM  net |g_OM|       .gMB  moving boundary, signed
%                      .gPE  photoelastic, signed (gMB + gPE = +/- gOM)
%                      .cancelRatio  (|gMB| + |gPE|)/gOM, > 1 = they oppose
%                      .rerun        the same four from the BEST re-run
%
%   Writes into OPT.rootLoc (keyed on OPT.runId, not the date):
%     eval_NNNN_obX_mdY/           per-evaluation FEM output and mode profiles
%     BEST/                        final re-run of the best design, plots on
%     optimize_log.csv             one row per evaluation, appended
%     optimize_convergence.png/fig fitness trace + sampled designs
%     bayesopt_results.mat         the BayesianOptimization object (bayesopt)
%   and into OPT.dbPath: the run store, plus a .jsonl mirror.
%
%   See also OPTIMIZE_OBLONG_MAXDEF, ISFABRICABLERIB, OPTIMSTORE,
%   SEEDBAYESFROMSTORE, CREATERIBCAVITYGEOM, RUNNANOBEAMFEM,
%   TEST_NANOBEAMRIBFEM.

if nargin < 1 || isempty(optOverrides)
    optOverrides = struct();
end

close all;

%% ===================== PATH =====================
% The shared rib helpers (ribDefaultParams, validateRibParams, ribLength)
% live in the omc-comsol-chang root. Added at the END of the path so this
% folder's own CreateFileBase.m / LoadMaterialParams.m still win.
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir)
    thisDir = pwd;
end
addpath(fileparts(thisDir), '-end');

%% ===================== TUNABLE KNOBS =====================

% --- mirror-cell rib dimensions (m) ---
% Defined here rather than only in the P struct below because the derived
% centre-cell aspect ratio (the pre-filter's input and the run store's key)
% depends on the mirror aspect ratio w_mirror/t. P.t / P.w_mirror are
% assigned from these.
t_0  = 280e-9;      % rib width along x at the mirror cells
wm_0 = 1300e-9;     % rib tip-to-tip width at the mirror cells

% --- starting point (the design of test_nanobeamRibFEM.m) ---
% w ~ (1-maxdef)^(1+oblong), t ~ (1-maxdef)^(1-oblong).
% oblong = 0 -> t and w scale like a; oblong < -1 -> ribs lengthen and
% narrow toward the centre; oblong > 1 -> ribs shorten and widen.
oblong_0 = -2;
maxdef_0 = 0.2;
mirrorAspectRatio   = wm_0 / t_0;                                        % 4.643
defectAspectRatio_0 = dARofOblong(oblong_0, maxdef_0, mirrorAspectRatio); % 11.335

% --- oblong box constraints (dimensionless taper exponent) ---
% Searched directly, on a LINEAR axis. Both taper directions are inside:
% oblong < 0 (ribs lengthen and narrow toward the centre, as the start
% design does) and oblong > 0 (ribs shorten and widen).
% The buildable oblong interval narrows as maxdef grows (base design):
%     maxdef <= 0.10 : wider than [-8, 5]     maxdef = 0.25 : [-3.35, 1.96]
%     maxdef  = 0.15 : [-6.70, 3.63]          maxdef = 0.30 : [-2.51, 1.54]
%     maxdef  = 0.20 : [-4.61, 2.59]
% bounded below by the tip fillet (t_ctr -> 2*r2) and above by the root
% fillet ((a - t)/2 -> r1 at the centre). So no box is both wide and fully
% buildable. With the default below 91.9 % of a 41x41 scan passes
% isFabricableRib; the rest (the two corners at large maxdef) is excluded
% by XConstraintFcn under bayesopt and by the graded penalty under
% Nelder-Mead, never solved.
% Alternatives: [-2.5, 1.5] is the widest box that is 100 % buildable
% (0.35 nm margin), but it stops 0.5 below the start design; [-3, 3] gives
% 91.7 %, [-5, 2] 89.1 %.
OPT.oblong_min = -4;
OPT.oblong_max = 2.5;

% --- maxdef box constraints (dimensionless, 0-<1) ---
% Absolute, not fractions of maxdef_0, so the box does not move if the start
% point is retuned.
OPT.maxdef_min = 0.05;
OPT.maxdef_max = 0.30;

% --- target mechanical frequency (solver eigen-shift, P.freq) ---
OPT.targetFreq = 7.2e9;    % Hz  (test_nanobeamRibFEM.m)

% --- optical-wavelength fitness target and tolerance ---
% The wavelength factor is GAUSSIAN: exp(-((lambda - target)/lambdaTol)^2),
% so lambdaTol is the 1/e half-width.
OPT.targetLambda = 1550;   % nm
OPT.lambdaTol    = 100;    % nm

% --- mechanical-Q gate: OFF by default (no rib PML, see the header) ---
% 1 = designs whose ESTIMATED Q_mech (low-reflecting beam end) is below
% OPT.QmechMin score 0. Needs P.mEndBC = 'lowreflecting'. OPT.QmechMin is
% the hole optimizer's PML value; recalibrate it from the Q_mech column of
% a first ungated run before switching the gate on.
OPT.useQmechGate = 0;
OPT.QmechMin     = 1e7;

% --- mechanical-frequency window: OFF by default ([]) ---
% [fLo fHi] in Hz multiplies the fitness by exp(-(df/mechFreqTol)^2), where
% df is how far the mechanical mode of the best g_OM pair lies outside the
% window (0 inside it).
OPT.mechFreqWindow = [];
OPT.mechFreqTol    = 0.5e9;    % Hz, Gaussian 1/e half-width of the roll-off

% --- objective penalty scale (Nelder-Mead ladder) ---
OPT.Jpenalty = 1e3;

% --- optimizer selection ---
% 'bayesopt'   : Gaussian-process Bayesian optimization (default).
% 'neldermead' : fminsearch simplex (base MATLAB, no toolbox).
OPT.optimizer = 'bayesopt';

% --- run persistence / stop-and-resume (see OptimStore.m) -----------------
%   OPT.runId  : names the study. Deliberately NOT date-stamped.
%   OPT.resume : 1 = serve already-evaluated points from the store.
%   OPT.storeBackend : '' auto-detects; 'dbtoolbox' | 'python' | 'cli' |
%                'jsonl' pins it. Pinned to 'jsonl' for the same reason as
%                in optimize_oblong_maxdef.m (Microsoft Store Python on the
%                remote workstation cannot be launched by MATLAB).
OPT.useStore     = 1;
OPT.storeBackend = 'jsonl';
OPT.runId        = 'ribOblongMaxdef_trial1';
OPT.dbPath       = fullfile('.', 'test', '1D_OMC_rib', 'optim_runs.sqlite3');
OPT.resume       = 1;

% --- output file names ---
% 1 = each evaluation's files are named rib_ob<oblong>_md<maxdef>...
% instead of CreateFileBase's full rib name (141 characters for the start
% design). With the full name, <rootLoc>/eval_.../<name>_wl_..._gO_....mat
% is over 300 characters from this repo's location, beyond Windows' 260
% character MAX_PATH, and the save fails. The full parameter set is still in
% ds.P inside each .mat, and the start design's CreateFileBase name is
% printed at startup. 0 = use CreateFileBase (needs long paths enabled).
OPT.shortFileBase = 1;

% --- dry run: no COMSOL, analytic stand-in, output under tempdir ---
OPT.dryRun = 0;

% --- bayesopt settings (used only when OPT.optimizer = 'bayesopt') --------
% maxEvals is the TOTAL size of the study, and points seeded from the store
% COUNT AGAINST IT. To extend a finished study, raise maxEvals. Every point
% not served from the store is a full COMSOL solve, so time one first.
OPT.bo.maxEvals       = 60;
OPT.bo.numSeedPoints  = 8;
OPT.bo.acquisition    = 'expected-improvement-plus';
% FALSE on purpose: adaptive meshing (P.mAdjMesh / P.oAdjMesh) makes a
% repeated evaluation of the same geometry differ slightly.
OPT.bo.deterministic  = false;
OPT.bo.plotFcn        = {@plotMinObjective, @plotObjectiveModel};
% 1 = on a cold start the start point above is the first evaluation, so
% every study contains the baseline design. (bayesopt would otherwise only
% draw random seed points.)
OPT.bo.evalStartFirst = 1;

% --- Nelder-Mead settings (used only when OPT.optimizer = 'neldermead') ---
% TolFun is an absolute tolerance on -log10(fitness), i.e. effectively a
% RELATIVE tolerance on the fitness: 1e-3 is about 0.23 %.
OPT.nm.TolX        = 1e-4;
OPT.nm.TolFun      = 1e-3;
OPT.nm.MaxIter     = 100;
OPT.nm.MaxFunEvals = 200;

% --- apply caller overrides (unknown names are an error) ---
OPT = mergeStruct(OPT, optOverrides, 'OPT');

% --- derived settings: do not edit ---
if OPT.dryRun
    % Never touch the real study: a fresh run id and a tempdir store / output
    % root, unless the caller names them explicitly (the resume test does).
    dryRoot = fullfile(tempdir, 'optimize_rib_oblong_maxdef_dryrun');
    boOver  = struct();
    nmOver  = struct();
    if isfield(optOverrides, 'bo'); boOver = optOverrides.bo; end
    if isfield(optOverrides, 'nm'); nmOver = optOverrides.nm; end
    if ~isfield(optOverrides, 'runId')
        OPT.runId = [OPT.runId, '_dryrun_', ...
            char(datetime('now', 'Format', 'yyyyMMddHHmmssSSS'))];
    end
    if ~isfield(optOverrides, 'dbPath')
        OPT.dbPath = fullfile(dryRoot, 'optim_runs.sqlite3');
    end
    if ~isfield(boOver, 'plotFcn');     OPT.bo.plotFcn = []; end
    if ~isfield(boOver, 'maxEvals');    OPT.bo.maxEvals = min(OPT.bo.maxEvals, 15); end
    if ~isfield(nmOver, 'MaxFunEvals'); OPT.nm.MaxFunEvals = min(OPT.nm.MaxFunEvals, 40); end
    rootBase = dryRoot;
else
    rootBase = fullfile('.', 'test', '1D_OMC_rib');
end

% Root folder for per-evaluation subfolders. Keyed on OPT.runId, NOT on the
% date, so a study resumed later keeps all of its output together.
OPT.rootLoc = [fullfile(rootBase, ['optimize_', OPT.runId]), filesep];

% --- FITNESS FUNCTION (user-configurable) ---------------------------------
%   fitness = g_OM^2 * Q_opt * exp(-((lambda_opt - targetLambda)/lambdaTol)^2)
%             [* mechanical-frequency window factor, if OPT.mechFreqWindow]
%   OPT.useQmechGate = 1 and estimated Q_mech < OPT.QmechMin  ->  fitness = 0
%   no usable optical mode / g_OM (or, gate on, no Q_mech)    ->  fitness = NaN
% Alternative: fitness_maxGOM, at the bottom of this file.
if OPT.useQmechGate
    QmechGate = OPT.QmechMin;
else
    QmechGate = -Inf;
end
OPT.fitnessFcn = @(ds) fitness_ribCoop(ds, OPT.targetLambda, OPT.lambdaTol, ...
    QmechGate, OPT.mechFreqWindow, OPT.mechFreqTol);

% The SAME figure of merit with the Q_mech gate disabled. Nelder-Mead uses
% the gated score; the GP surrogate needs the ungated value, and the gate
% becomes a coupled constraint instead. With OPT.useQmechGate = 0 the two
% are identical. Both are stored, so the optimizers can share a run store.
OPT.fitnessFcnRaw = @(ds) fitness_ribCoop(ds, OPT.targetLambda, OPT.lambdaTol, ...
    -Inf, OPT.mechFreqWindow, OPT.mechFreqTol);

%% ===================== P STRUCT DEFAULTS =====================
% The settings of test_nanobeamRibFEM.m. Field meanings: rib_README.md.
P.xsect    = 'rect';
P.celltype = 'rib';
P.beamMat  = 'diamond';
P.anisoMat = 1;
P.rxtal    = 0;
P.rxtalInFilename = 1;

P.a  = 600e-9;      % mirror lattice constant (cell width)
P.t  = t_0;         % rib width along x - set in TUNABLE KNOBS
P.s  = 100e-9;      % spine width along y
P.th = 250e-9;      % slab thickness
P.r1 = 40e-9;       % concave rib/spine (root) fillet
P.r2 = 40e-9;       % convex rib-tip fillet
P.L_ref = 400e-9;   % L(a) rule; unused while P.w_mirror is set
P.a_ref = 240e-9;
P.dLda  = -0.4998;
P.w_mirror = wm_0;  % rib tip-to-tip width at the mirror cells - TUNABLE KNOBS
P.d = 0;            % sidewall angle (rad); d ~= 0 needs r1 = r2 = 0

P.ncellHalf = 14;
P.ribatctr  = 1;    % rib at the cavity centre
P.ndef      = 5;
P.maxdef    = maxdef_0;     % start value - overwritten every evaluation
P.taperFunc = 'quadratic';
P.taperConvention = 'centerspacing';
P.ribTaperScheme  = 'hole';     % required: this is what P.oblong acts on
P.oblong    = oblong_0;     % start value - overwritten every evaluation

P.lambda = OPT.targetLambda*1e-9;   % optical eigen-search target (m)
P.nbeam  = 2.386;

P.asymCav = 0;
P.asym    = 0;

% --- physics toggles ---
P.solveMech    = 1;
P.solveOpt     = 1;
P.calcG        = 1 * (P.solveMech && P.solveOpt);
P.calcS        = 0;
P.solveMechPML = 0;     % not supported by buildRibCavity (buildRibCavity.m:60-62)
P.extractLocMechModes = 1;

% --- plot / save: mode profiles saved to each per-eval subfolder ---
P.plotgeom   = 0;
P.storeMPH   = 0;
P.plotMech   = 1 * P.solveMech;
P.plotOpt    = 1 * P.solveOpt;
P.plotStrCpl = 1 * P.calcS;

% --- mechanical solver ---
P.mevenx  = 1;
P.meveny  = 1;
P.mevenz  = 1;
P.freq    = OPT.targetFreq;
P.mneigs  = 20;
P.mMesh   = 1;
P.mAdjMesh = 1;
P.mEndBC  = 'lowreflecting';    % beam-end BC: 'lowreflecting' or 'fixed'
P.max_dof = 5e6;

% --- optical solver ---
P.oevenx  = 1;      % forced to 1 by SetupNanobeamFEM anyway when ribatctr = 1
P.oeveny  = -1;     % TE-like
P.oevenz  = 1;
P.oneigs  = 1;
P.oMesh   = 3;
P.oAdjMesh = 1;
P.airrad  = 3.1e-6;
P.airZOffset = 0;   % an axis offset inside the slab breaks the mesh

P.g0min = 80e3;

% SiV strain coupling axes (read by RunNanobeamFEM even with calcS = 0)
P.zSiV  = {[1 1 1]};
P.xSlc  = 0;
P.ySlc  = 0;
P.zSlc  = 0;

%% ===================== PRE-FLIGHT (pure MATLAB) =====================
% Configuration mistakes should stop the study before the first solve, not
% turn every evaluation into a caught error.
if ~any(strcmpi(OPT.optimizer, {'bayesopt', 'neldermead'}))
    error('optimize_rib_oblong_maxdef:unknownOptimizer', ...
        ['Unknown OPT.optimizer "%s". Valid values are ''neldermead'' ' ...
         'and ''bayesopt''.'], OPT.optimizer);
end
if ~(isfield(P, 'ribTaperScheme') && strcmp(P.ribTaperScheme, 'hole'))
    error('optimize_rib_oblong_maxdef:needsHoleScheme', ...
        'P.ribTaperScheme must be ''hole'': P.oblong has no effect in the ''rib'' scheme.');
end
if P.solveMechPML
    error('optimize_rib_oblong_maxdef:noRibPML', ...
        ['P.solveMechPML = 1, but buildRibCavity does not build a mechanical ' ...
         'PML (buildRibCavity.m:60-62). Set it to 0. Once rib PML support ' ...
         'exists, remove this check: local_getQmechRib already prefers ' ...
         'ds.mfem.QAll.']);
end
if ~(isempty(OPT.mechFreqWindow) || (isnumeric(OPT.mechFreqWindow) && ...
        numel(OPT.mechFreqWindow) == 2 && OPT.mechFreqWindow(1) < OPT.mechFreqWindow(2)))
    error('optimize_rib_oblong_maxdef:badFreqWindow', ...
        'OPT.mechFreqWindow must be [] or [fLo fHi] in Hz with fLo < fHi.');
end
if OPT.useQmechGate
    if ~strcmp(P.mEndBC, 'lowreflecting')
        error('optimize_rib_oblong_maxdef:noQmech', ...
            ['OPT.useQmechGate = 1 needs a mechanical loss channel, and with ' ...
             'P.mEndBC = ''%s'' and no PML the eigenvalues are lossless. Use ' ...
             'P.mEndBC = ''lowreflecting'' or set OPT.useQmechGate = 0.'], P.mEndBC);
    end
    warning('optimize_rib_oblong_maxdef:QmechEstimate', ...
        ['OPT.useQmechGate = 1: Q_mech is an ESTIMATE from the low-reflecting ' ...
         'beam end (a first-order absorbing boundary), not a PML result. ' ...
         'OPT.QmechMin = %.3g should be calibrated against an ungated run.'], ...
        OPT.QmechMin);
end

% The start design must build, and the knobs must describe the same mirror
% cell as P does (they drift apart if P.w_mirror / P.t are edited directly).
try
    Pchk = createRibCavityGeom(P);
catch ME
    error('optimize_rib_oblong_maxdef:startInvalid', ...
        'The start design (oblong_0 = %g, maxdef_0 = %g) does not build: %s', ...
        oblong_0, maxdef_0, ME.message);
end
Pdef  = ribDefaultParams(P);
wmChk = 2*ribLength(Pdef, Pdef.a);
if abs(wmChk/Pdef.t - mirrorAspectRatio) > 1e-9*mirrorAspectRatio
    error('optimize_rib_oblong_maxdef:mappingDrift', ...
        ['P gives a mirror aspect ratio w/t = %.6f, but the knobs wm_0/t_0 ' ...
         'give %.6f. Set the mirror cell through t_0 and wm_0 only.'], ...
        wmChk/Pdef.t, mirrorAspectRatio);
end
[x0ok, ~, x0why] = isFabricableRib(defectAspectRatio_0, maxdef_0, P);
if ~(isscalar(OPT.oblong_min) && isscalar(OPT.oblong_max) && ...
        isfinite(OPT.oblong_min) && isfinite(OPT.oblong_max) && ...
        OPT.oblong_min < OPT.oblong_max)
    error('optimize_rib_oblong_maxdef:badOblongBounds', ...
        'OPT.oblong_min and OPT.oblong_max must be finite scalars with min < max.');
end
x0inBox = oblong_0 >= OPT.oblong_min && oblong_0 <= OPT.oblong_max && ...
          maxdef_0 >= OPT.maxdef_min && maxdef_0 <= OPT.maxdef_max;
if ~x0ok
    warning('optimize_rib_oblong_maxdef:startUnfabricable', ...
        'The start design fails the fabrication pre-filter (%s).', x0why{1});
end

% Output path length (Windows MAX_PATH; see OPT.shortFileBase).
PfullName = CreateFileBase(Pchk);
if OPT.shortFileBase
    baseLen = numel(shortFileBase('', -99.9999, 0.9999));
else
    baseLen = numel(PfullName.fileBase);
end
if strncmp(OPT.rootLoc, '.', 1)
    absRoot = fullfile(pwd, OPT.rootLoc);
else
    absRoot = OPT.rootLoc;
end
% eval folder + file base + RunNanobeamFEM's '_wl_..nm_wM_..GHz_gO_..kHz.mat'
worstPathLen = numel(absRoot) + numel('eval_0000_ob-00.0000_md0.0000') + 1 + baseLen + 45;
if ispc && worstPathLen > 259
    warning('optimize_rib_oblong_maxdef:pathTooLong', ...
        ['Output paths reach about %d characters, beyond Windows'' 260 ' ...
         'character limit, so saving results will fail unless long paths ' ...
         'are enabled. Set OPT.shortFileBase = 1, shorten OPT.runId, or ' ...
         'run from a shorter folder.'], worstPathLen);
end

%% ===================== OPEN CSV LOG =====================
% APPEND, never truncate: OPT.rootLoc is keyed on OPT.runId, so a resumed
% session reopens the SAME file. The first 13 columns are the hole
% optimizer's; the rib-specific ones follow.
if ~exist(OPT.rootLoc, 'dir'); mkdir(OPT.rootLoc); end
logFile = [OPT.rootLoc, 'optimize_log.csv'];
newLog  = (exist(logFile, 'file') ~= 2);
fid = fopen(logFile, 'a');
if newLog
    fprintf(fid, ['eval,oblong,maxdef,wM_GHz,gOM_kHz,gMB_kHz,gPE_kHz,' ...
                  'LSiV_MHz,Q_mech,Q_opt,lambda_opt_nm,fitness,status,' ...
                  'dAR,t_ctr_nm,w_ctr_nm\n']);
else
    % Mark the seam so a multi-session log can be read back unambiguously.
    fprintf(fid, '# resumed %s\n', ...
        char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
end
fclose(fid);

% Counters and trackers shared with the nested objective.
evalCount   = 0;    % highest evaluation index in use
nSolved     = 0;    % evaluations actually computed this session
nCacheHits  = 0;    % evaluations served from the store
nInfeasible = 0;    % candidates rejected by the pre-filter (no solve)

% Best (positive) fitness seen so far, for the final re-run.
bestSoFar.fitness = -Inf;
bestSoFar.oblong  = NaN;
bestSoFar.maxdef  = NaN;
bestSoFar.defectAspectRatio = NaN;
% Its coupling in Hz: net, and the signed moving-boundary and photoelastic
% parts of the same mode pair (see local_getGOM).
bestSoFar.gOM = NaN;
bestSoFar.gMB = NaN;
bestSoFar.gPE = NaN;

% Full search history, appended by the objective and plotted after the run.
history.eval    = [];
history.dAR     = [];
history.maxdef  = [];
history.oblong  = [];
history.fitness = [];
history.Qmech   = [];
history.J       = [];

%% ===================== RUN STORE (stop / resume) ==================
store = [];
if OPT.useStore
    store = OptimStore(OPT.dbPath, OPT.storeBackend);
    % jsonencode cannot represent a cell array of function handles, so the
    % plot-function list is stringified for the configuration snapshot.
    cfgSnap = OPT;
    % For the record only: which variables this session searched. Resume
    % does not depend on it, because the store is keyed on (dAR, maxdef).
    cfgSnap.searchSpace = 'oblong';
    if iscell(cfgSnap.bo.plotFcn)
        cfgSnap.bo.plotFcn = cellfun(@func2str, cfgSnap.bo.plotFcn, ...
            'UniformOutput', false);
    end
    store.startRun(OPT.runId, 'optimize_rib_oblong_maxdef', OPT.optimizer, cfgSnap);

    nPrior = store.evalCount(OPT.runId);
    fprintf('Run store    : %s  [backend: %s]\n', store.dbPath, store.backend);
    fprintf('Run id       : %s  (%d evaluation(s) already stored)\n', ...
        OPT.runId, nPrior);

    if nPrior > 0 && OPT.resume
        % Seed the incumbent from what the earlier session already found,
        % so a resumed run can never report a worse design than one that
        % has already been paid for.
        % The store is keyed on x = [defectAspectRatio, maxdef].
        [xPrior, fPrior, iPrior] = store.best(OPT.runId);
        if ~isempty(xPrior)
            bestSoFar.fitness = fPrior;
            bestSoFar.defectAspectRatio = xPrior(1);
            bestSoFar.maxdef  = xPrior(2);
            bestSoFar.oblong  = oblongOfDAR(xPrior(1), xPrior(2), mirrorAspectRatio);
            hitPrior = store.lookup(OPT.runId, xPrior);
            if ~isempty(hitPrior)
                bestSoFar.gOM = hitPrior.gOM;
                bestSoFar.gMB = hitPrior.gMB;
                bestSoFar.gPE = hitPrior.gPE;
            end
            fprintf(['Resuming     : best stored fitness %.4g at eval %d ' ...
                     '(oblong=%.4f, maxdef=%.4f, dAR=%.4f)\n' ...
                     '               gOM = %s\n'], ...
                fPrior, iPrior, bestSoFar.oblong, xPrior(2), xPrior(1), ...
                cplText(cplSummary(bestSoFar.gOM, bestSoFar.gMB, bestSoFar.gPE)));
            if bestSoFar.oblong < OPT.oblong_min - 1e-9 || ...
                    bestSoFar.oblong > OPT.oblong_max + 1e-9 || ...
                    xPrior(2) < OPT.maxdef_min || xPrior(2) > OPT.maxdef_max
                fprintf(['               NOTE: that design lies OUTSIDE the ' ...
                         'current search box. It stays the incumbent\n' ...
                         '               until a better one is found; widen ' ...
                         'the bounds to search around it.\n']);
            end
        else
            fprintf('Resuming     : no successful evaluation stored yet\n');
        end
        % Keep new indices clear of the ones already in the store.
        evalCount = nPrior;
    elseif nPrior > 0
        fprintf(['NOT resuming : OPT.resume = 0, so all %d stored point(s) ' ...
                 'will be re-solved\n'], nPrior);
        evalCount = nPrior;
    end

    % Mark the run interrupted if this function exits without reaching the
    % end -- an error, or a Ctrl-C. The close-out block overwrites it.
    runGuard = onCleanup(@() store.finishRun(OPT.runId, 'interrupted'));
end

%% ===================== START-UP REPORT =====================
x0 = [oblong_0, maxdef_0];

fprintf('\n==================== OPTIMIZATION START ====================\n');
if OPT.dryRun
    fprintf('*** DRY RUN: analytic stand-in for RunNanobeamFEM, NOT physics ***\n');
end
fprintf('Optimizer                = %s\n', OPT.optimizer);
fprintf('Mirror cell              : t = %.0f nm, w = %.0f nm (w/t = %.4f)\n', ...
    t_0*1e9, wm_0*1e9, mirrorAspectRatio);
fprintf('Search variables         = oblong, maxdef\n');
fprintf('Start oblong             = %.6f  (centre-cell w/t = %.4f)\n', oblong_0, defectAspectRatio_0);
fprintf('Start maxdef             = %.6f\n', maxdef_0);
fprintf('oblong bounds            = [%.6f, %.6f]\n', ...
    OPT.oblong_min, OPT.oblong_max);
fprintf('maxdef bounds            = [%.6f, %.6f]\n', ...
    OPT.maxdef_min, OPT.maxdef_max);
if OPT.useQmechGate
    fprintf('Q_mech gate              = %.3g on the low-reflecting-end ESTIMATE\n', OPT.QmechMin);
else
    fprintf('Q_mech gate              = off (no rib PML; Q_mech estimate is logged only)\n');
end
if isempty(OPT.mechFreqWindow)
    fprintf('Mech. frequency window   = off\n');
else
    fprintf('Mech. frequency window   = [%.2f, %.2f] GHz (roll-off %.2f GHz)\n', ...
        OPT.mechFreqWindow(1)*1e-9, OPT.mechFreqWindow(2)*1e-9, OPT.mechFreqTol*1e-9);
end
fprintf('Start design file name   : %s\n', PfullName.fileBase);
fprintf('Log written -> %s\n\n', logFile);

%% ===================== FABRICATION FEASIBILITY SCAN ==================
% isFabricableRib reproduces the builder's validity rules and the 50 nm
% lithography floors as arithmetic. Scanning the box up front says how much
% of the search space is reachable before any COMSOL time is spent. The grid
% is uniform in (oblong, maxdef), which is also how bayesopt draws its seed
% points, so the fraction is the share of random draws that are admissible.
% Part of the box being unbuildable is EXPECTED in these variables (see the
% oblong bounds above); the report then names the rules that exclude it.
nScan = 41;
[scanOB, scanMD] = meshgrid(linspace(OPT.oblong_min, OPT.oblong_max, nScan), ...
                            linspace(OPT.maxdef_min, OPT.maxdef_max, nScan));
[scanOK, scanMargin, scanWhy] = isFabricableRib( ...
    dARofOblong(scanOB, scanMD, mirrorAspectRatio), scanMD, P);
fracFab = mean(scanOK(:));
if all(scanOK(:))
    fprintf('Fabricable   : 100.0%% of the search box (tightest margin %.1f nm)\n', ...
        min(scanMargin(:))*1e9);
else
    scanRules = unique(scanWhy(~scanOK));
    scanText  = cell(size(scanRules));
    for iRule = 1:numel(scanRules)
        scanText{iRule} = sprintf('%s %.1f%%', scanRules{iRule}, ...
            100*nnz(strcmp(scanWhy(:), scanRules{iRule}))/numel(scanOK));
    end
    fprintf('Fabricable   : %.1f%% of the search box (excluded: %s)\n', ...
        100*fracFab, strjoin(scanText(:)', ', '));
end
if fracFab == 0
    error('optimize_rib_oblong_maxdef:boxInfeasible', ...
        ['No point on a %dx%d scan of the bounds is fabricable, so the ' ...
         'study was NOT started. Change the bounds, or relax ' ...
         'P.minRibFeature / P.minRibGap / P.minRibThickness / P.ribFabStrict.'], ...
        nScan, nScan);
elseif fracFab < 0.25
    warning('optimize_rib_oblong_maxdef:boxMostlyInfeasible', ...
        ['Only %.1f%% of the search box is fabricable. The optimizer can ' ...
         'work here, but consider tightening the bounds around the ' ...
         'feasible region.'], 100*fracFab);
end

%% ===================== OPTIMIZER DISPATCH =====================
switch lower(OPT.optimizer)
    case 'neldermead'
        % CAVEAT: Nelder-Mead has no noise handling of its own. If the run
        % converges onto a point whose fitness does not reproduce when
        % re-run, that is mesh noise, not an optimum.
        options = optimset('Display', 'iter', ...
                           'TolX',        OPT.nm.TolX, ...
                           'TolFun',      OPT.nm.TolFun, ...
                           'MaxIter',     OPT.nm.MaxIter, ...
                           'MaxFunEvals', OPT.nm.MaxFunEvals);
        % Wrapped in try/catch so the optimizer never crashes.
        [xBest, JBest] = fminsearch(@safeObjective, x0, options);

    case 'bayesopt'
        % A bare 'Unrecognized function ''bayesopt''' reads like a licence
        % problem when it is usually a missing install.
        if exist('bayesopt', 'file') == 0 && exist('bayesopt', 'builtin') == 0
            error('optimize_rib_oblong_maxdef:noBayesopt', ...
                ['bayesopt is not available in this MATLAB install.\n' ...
                 '  licence permits it (1=yes) : %d\n' ...
                 '  installed on disk (1=yes)  : %d\n' ...
                 'Licence 1 with install 0 means: install the Statistics ' ...
                 'and Machine Learning Toolbox from the Add-On Explorer. ' ...
                 'Or set OPT.optimizer = ''neldermead''.'], ...
                license('test', 'Statistics_Toolbox'), ...
                isfolder(fullfile(matlabroot, 'toolbox', 'stats')));
        end

        % The box is the variable range, so bayesopt never proposes an
        % out-of-bounds point. oblong is on a linear axis (no transform).
        optVars = [ ...
            optimizableVariable('oblong', ...
                [OPT.oblong_min, OPT.oblong_max], 'Type', 'real'), ...
            optimizableVariable('maxdef', ...
                [OPT.maxdef_min, OPT.maxdef_max], 'Type', 'real')];

        % Coupled constraints (<= 0 means satisfied):
        %   gate on  : [Q_mech gate, usable solve]
        %   gate off : [usable solve]
        % seedBayesFromStore always returns both columns; without the gate
        % only the second is meaningful (its first treats a missing Q_mech
        % as a violation), so that one is dropped.
        [initX, initObj, initCon] = seedBayesFromStore(store, OPT);
        % seedBayesFromStore (shared with the hole optimizer) returns the
        % store's own coordinates, (dAR, maxdef). Map them to the search
        % variables; the mapping is exact, so rows written by a session that
        % searched dAR seed this one just as well. A value that comes back a
        % rounding error outside the box is snapped onto its edge.
        seedOb = oblongOfDAR(initX.dAR, initX.maxdef, mirrorAspectRatio);
        seedOb(abs(seedOb - OPT.oblong_min) < 1e-9) = OPT.oblong_min;
        seedOb(abs(seedOb - OPT.oblong_max) < 1e-9) = OPT.oblong_max;
        nSeedOut = nnz(seedOb < OPT.oblong_min | seedOb > OPT.oblong_max);
        initX = table(seedOb, initX.maxdef, 'VariableNames', {'oblong', 'maxdef'});
        if OPT.useQmechGate
            numCon = 2;
        else
            numCon  = 1;
            initCon = initCon(:, 2);
        end
        nSeeded  = height(initX);
        seedLeft = max(0, OPT.bo.numSeedPoints - nSeeded);
        if nSeeded > 0
            fprintf(['Seeding GP with %d stored point(s); %d fresh seed ' ...
                     'point(s) to go.\n'], nSeeded, seedLeft);
        end
        if nSeedOut > 0
            % bayesopt warns about these itself, and still uses them.
            fprintf(['               %d of them lie outside the oblong box: ' ...
                     'they inform the surrogate but cannot be\n' ...
                     '               returned by bayesopt as its best point.\n'], ...
                nSeedOut);
        end
        if nSeeded >= OPT.bo.maxEvals
            warning('optimize_rib_oblong_maxdef:budgetExhausted', ...
                ['The store already holds %d point(s) for run "%s", ' ...
                 'which meets or exceeds OPT.bo.maxEvals = %d. bayesopt ' ...
                 'counts seeded points against its budget, so NO new ' ...
                 'evaluations will be made. Raise OPT.bo.maxEvals to ' ...
                 'extend the study.'], ...
                nSeeded, OPT.runId, OPT.bo.maxEvals);
        end

        % XConstraintFcn is evaluated on CANDIDATES, before the objective is
        % called at all, so an unfabricable design is never proposed and
        % never costs a solve. isFabricableRib is vectorized because
        % bayesopt passes a table of many candidate rows at once.
        xConFcn = @(t) isFabricableRib( ...
            dARofOblong(t.oblong, t.maxdef, mirrorAspectRatio), t.maxdef, P);

        boArgs = {'MaxObjectiveEvaluations',  OPT.bo.maxEvals, ...
                  'XConstraintFcn',           xConFcn, ...
                  'NumSeedPoints',            max(1, seedLeft), ...
                  'AcquisitionFunctionName',  OPT.bo.acquisition, ...
                  'IsObjectiveDeterministic', OPT.bo.deterministic, ...
                  'NumCoupledConstraints',    numCon, ...
                  'AreCoupledConstraintsDeterministic', false(1, numCon), ...
                  'PlotFcn',                  OPT.bo.plotFcn, ...
                  'Verbose',                  1};
        if nSeeded > 0
            boArgs = [boArgs, {'InitialX', initX, ...
                               'InitialObjective', initObj, ...
                               'InitialConstraintViolations', initCon}];
        elseif OPT.bo.evalStartFirst && x0ok && x0inBox
            % InitialX without InitialObjective: bayesopt evaluates it first.
            boArgs = [boArgs, {'InitialX', ...
                table(oblong_0, maxdef_0, ...
                      'VariableNames', {'oblong', 'maxdef'})}];
        end

        boResults = bayesopt(@bayesObjective, optVars, boArgs{:});

        % Take bayesopt's feasible best as the reported point and let the
        % resolve-best block below reconcile it against bestSoFar, which
        % both optimizers maintain identically. JBest = NaN routes it there.
        try
            tBest = bestPoint(boResults);
            xBest = [tBest.oblong, tBest.maxdef];
        catch
            warning('optimize_rib_oblong_maxdef:noFeasibleBO', ...
                ['bayesopt found no point satisfying the constraints. ' ...
                 'Falling back to the best evaluated point.']);
            xBest = x0;
        end
        if isempty(xBest)
            xBest = x0;
        end
        JBest = NaN;

        % Keep the study object with the run's other artifacts.
        try
            save([OPT.rootLoc, 'bayesopt_results.mat'], 'boResults');
        catch ME
            warning('optimize_rib_oblong_maxdef:boSaveFailed', ...
                'Could not save the bayesopt results object: %s', ME.message);
        end
end

%% ===================== RESOLVE BEST POINT =====================
best_oblong            = xBest(1);
best_maxdef            = xBest(2);
best_defectAspectRatio = dARofOblong(best_oblong, best_maxdef, mirrorAspectRatio);

% Undo the log scaling. A J at or above the penalty floor means fminsearch
% finished on an out-of-bounds, gated or failed vertex, which carries no
% fitness. The bayesopt branch sets JBest = NaN, which lands in the same
% place: the incumbent is resolved from bestSoFar below.
if isnan(JBest) || JBest >= OPT.Jpenalty
    fitnessBest = NaN;
else
    fitnessBest = 10^(-JBest);
end
% bestSoFar also holds the stored incumbent of a resumed run, which
% fminsearch's own final vertex may not have revisited.
if isfinite(bestSoFar.fitness) && ...
        (~isfinite(fitnessBest) || fitnessBest <= 0 || bestSoFar.fitness > fitnessBest)
    best_defectAspectRatio = bestSoFar.defectAspectRatio;
    best_maxdef            = bestSoFar.maxdef;
    best_oblong            = bestSoFar.oblong;
    fitnessBest            = bestSoFar.fitness;
end

% Coupling of the REPORTED design. bestSoFar carries it whenever the
% incumbent is the design reported; otherwise the optimizer's own endpoint
% is, and the store has its row.
bestCpl = cplSummary(NaN, NaN, NaN);
if bestSoFar.defectAspectRatio == best_defectAspectRatio && ...
        bestSoFar.maxdef == best_maxdef
    bestCpl = cplSummary(bestSoFar.gOM, bestSoFar.gMB, bestSoFar.gPE);
elseif OPT.useStore
    hitBest = store.lookup(OPT.runId, [best_defectAspectRatio, best_maxdef]);
    if ~isempty(hitBest)
        bestCpl = cplSummary(hitBest.gOM, hitBest.gMB, hitBest.gPE);
    end
end
bestCpl.rerun = cplSummary(NaN, NaN, NaN);   % filled after the BEST re-run

stats = struct('nEval', evalCount, 'nSolved', nSolved, ...
               'nCacheHits', nCacheHits, 'nInfeasible', nInfeasible, ...
               'bestCpl', bestCpl);

fprintf('\n==================== OPTIMIZATION COMPLETE ====================\n');
fprintf('Best oblong            = %.4f\n', best_oblong);
fprintf('Best maxdef            = %.4f\n', best_maxdef);
fprintf('Best defectAspectRatio = %.6f  (centre-cell w/t, derived)\n', best_defectAspectRatio);
fprintf('Best fitness           = %.4g\n', fitnessBest);
fprintf('Best gOM               = %s\n', cplText(bestCpl));
fprintf('Highest eval index     = %d\n', evalCount);
fprintf('Solved this session    = %d\n', nSolved);
fprintf('Rejected unfabricable  = %d (no solve attempted)\n', nInfeasible);
if OPT.useStore
    fprintf('Served from store      = %d (solves skipped)\n', nCacheHits);
    fprintf('Run store              -> %s [%s]\n', store.dbPath, store.backend);
end
fprintf('Log written            -> %s\n', logFile);

%% ===================== CONVERGENCE PLOTS =====================
if isempty(history.eval)
    warning('optimize_rib_oblong_maxdef:noEvals', ...
        'No evaluations were made or replayed; skipping convergence plots.');
else
    % Gated (fit == 0) and failed (fit == NaN) points cannot go on a log axis.
    validFit = history.fitness;
    validFit(~isfinite(validFit) | validFit <= 0) = NaN;

    runningBest = nan(size(validFit));
    bestYet     = -Inf;
    for iHist = 1:numel(validFit)
        if isfinite(validFit(iHist)) && validFit(iHist) > bestYet
            bestYet = validFit(iHist);
        end
        if isfinite(bestYet); runningBest(iHist) = bestYet; end
    end

    figConv = figure('Name', [OPT.optimizer, ' convergence (rib)'], 'Color', 'w', ...
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
        % mostly gated/failed is obvious rather than showing as gaps.
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
    % meaningful; bayesopt jumps wherever the acquisition function points.
    isSimplex = strcmpi(OPT.optimizer, 'neldermead');
    nexttile;
    legendText = {'evals', 'start', 'best', 'bounds'};
    if isSimplex
        plot(history.oblong, history.maxdef, '-', 'Color', [0.7 0.7 0.7]);
        hold on;
        legendText = [{'path'}, legendText];
    end
    scatter(history.oblong, history.maxdef, 34, history.eval, 'filled');
    hold on;
    cb = colorbar;
    cb.Label.String = 'evaluation number';
    plot(oblong_0, maxdef_0, 'ks', 'MarkerSize', 11, 'LineWidth', 1.4);
    plot(best_oblong, best_maxdef, 'r*', 'MarkerSize', 15, 'LineWidth', 1.6);
    plot([OPT.oblong_min, OPT.oblong_max, OPT.oblong_max, ...
          OPT.oblong_min, OPT.oblong_min], ...
         [OPT.maxdef_min, OPT.maxdef_min, OPT.maxdef_max, ...
          OPT.maxdef_max, OPT.maxdef_min], 'k--', 'LineWidth', 1);
    if ~all(scanOK(:))
        % Edge of the fabricable region inside the box, from the start-up
        % scan: the corners beyond it are never solved.
        contour(scanOB, scanMD, double(scanOK), [0.5 0.5], ':', ...
            'LineColor', [0.85 0.33 0.10], 'LineWidth', 1.2);
        legendText = [legendText, {'fabricable limit'}];
    end
    grid on;
    xlabel('oblong');
    ylabel('maxdef');
    if isSimplex
        title('Simplex trajectory');
    else
        title('Sampled designs');
    end
    legend(legendText, 'Location', 'best');

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
if isfield(Pfinal, 'fileBase'); Pfinal = rmfield(Pfinal, 'fileBase'); end
if OPT.shortFileBase
    Pfinal.fileBase = shortFileBase('BEST_', best_oblong, best_maxdef);
end

bestLoc = [OPT.rootLoc, 'BEST', filesep];
if OPT.dryRun
    fprintf('\nDry run: skipping the COMSOL re-run of the best design.\n');
    dsBest = dryRunDs(Pfinal, best_defectAspectRatio, OPT);
else
    fprintf('\nRe-running best design with plots -> %s\n', bestLoc);
    [dsBest, ~] = RunNanobeamFEM(Pfinal, bestLoc);
end

% Coupling of the re-run itself. It is a fresh solve, so it need not
% reproduce the search-time numbers above exactly (mesh noise).
[gOMre, gMBre, gPEre] = local_getGOM(dsBest);
stats.bestCpl.rerun = cplSummary(gOMre, gMBre, gPEre);
fprintf('BEST re-run gOM        = %s\n', cplText(stats.bestCpl.rerun));

%% ===================== CLOSE OUT THE RUN =====================
if OPT.useStore
    % Disarm the interrupted-marker BEFORE writing 'complete': clearing an
    % onCleanup fires its callback, so the other order would let the guard
    % overwrite the completion status on the way out.
    clear runGuard;
    store.finishRun(OPT.runId, 'complete');
    fprintf(['\nRun "%s" marked complete. To extend it, raise ' ...
             'OPT.bo.maxEvals (or OPT.nm.MaxFunEvals) and re-run\nwith the ' ...
             'same OPT.runId -- the %d stored point(s) are reused for free.\n'], ...
        OPT.runId, store.evalCount(OPT.runId));
end

%% ===================== NESTED OBJECTIVE =====================
    function J = safeObjective(x)
    % Wrapper around objective() so that nothing can abort the optimization.
        try
            J = objective(x);
        catch ME
            warning('optimize_rib_oblong_maxdef:objError', ...
                'objective failed at x=[%.4f, %.4f]: %s', x(1), x(2), ME.message);
            J = 2 * OPT.Jpenalty;
        end
    end

    function [obj, con] = bayesObjective(t)
    % bayesopt adapter. bayesopt passes a one-row TABLE, and wants the
    % objective plus a vector of coupled-constraint values (<= 0 means
    % satisfied). Penalties are NOT folded into the objective: the Nelder-
    % Mead ladder (1000, 2000) next to feasible values around -16 would
    % wreck the GP's length scale.
    %   usable solve : -1 on success, +1 when there is no usable mode
    %   Q_mech gate  : 1 - Qmech/QmechMin   (only with OPT.useQmechGate)
    % A NaN objective marks an evaluation bayesopt could not use.
        R = evalPoint([t.oblong, t.maxdef]);

        usable = isfinite(R.fitRaw) && R.fitRaw > 0;
        if usable
            obj = -log10(R.fitRaw);
            conUsable = -1;
        else
            obj = NaN;
            conUsable = 1;
        end

        if OPT.useQmechGate
            if isfinite(R.Qmech)
                conQ = 1 - R.Qmech/OPT.QmechMin;
            else
                conQ = 1;   % no Q_mech at all counts as violating the gate
            end
            con = [conQ, conUsable];
        else
            con = conUsable;
        end
    end

    function J = objective(x)
    % Nelder-Mead adapter: sloped box barrier, then the penalty ladder.
    % The out-of-bounds penalty grows with the distance outside the box
    % (normalized by each parameter's range), so the simplex always has a
    % direction back into the feasible set. No FEM call is made.
        obRange = OPT.oblong_max - OPT.oblong_min;
        mdRange = OPT.maxdef_max - OPT.maxdef_min;
        viol = max([(OPT.oblong_min - x(1)) / obRange, ...
                    (x(1) - OPT.oblong_max) / obRange, ...
                    (OPT.maxdef_min - x(2)) / mdRange, ...
                    (x(2) - OPT.maxdef_max) / mdRange, ...
                    0]);
        if viol > 0
            J = OPT.Jpenalty * (2 + viol);
            return;
        end

        R = evalPoint(x);
        J = R.J;
    end

    function R = evalPoint(x)
    % Shared evaluation core for BOTH optimizers: x = [oblong, maxdef]. It
    % serves the point from the store if it was already solved, otherwise
    % runs RunNanobeamFEM (or the dry-run stand-in), then logs to the CSV,
    % the store and the plotting history.
    %   R.fit    - GATED fitness (0 below OPT.QmechMin, gate on) -- Nelder-Mead
    %   R.fitRaw - UNGATED fitness                               -- bayesopt
    %   R.Qmech  - mechanical Q estimate                 -- the BO constraint
    %   R.J      - Nelder-Mead penalty-ladder objective
        oblong_i = x(1);
        maxdef_i = x(2);
        % Centre-cell w/t: derived, not searched. It is the pre-filter's
        % input and, with maxdef, the run store's key (see logEval below).
        defectAspectRatio_i = dARofOblong(oblong_i, maxdef_i, mirrorAspectRatio);
        R = struct('fit', NaN, 'fitRaw', NaN, 'Qmech', NaN, ...
                   'J', 2*OPT.Jpenalty, 'cached', false, ...
                   'fabricable', true);

        % --- pre-filter: cheapest rejection there is ----------------------
        % bayesopt also has this as its XConstraintFcn, so under 'bayesopt'
        % this branch is a backstop rather than the main line.
        [fabOK, fabMargin, fabWhy] = isFabricableRib( ...
            defectAspectRatio_i, maxdef_i, P);
        if ~fabOK
            nInfeasible = nInfeasible + 1;
            % Graded like the out-of-bounds barrier, and on the same tier.
            % Normalizing by 50 nm keeps the violation O(1).
            viol   = min(10, -fabMargin / 50e-9);
            R.J    = OPT.Jpenalty * (2 + viol);
            R.fabricable = false;
            fprintf(['  eval  --: oblong=%.4f  maxdef=%.4f  (dAR=%.4f)  ' ...
                     '[UNFABRICABLE: %s, %.1f nm short]\n'], ...
                oblong_i, maxdef_i, defectAspectRatio_i, fabWhy{1}, -fabMargin*1e9);
            return;
        end

        % --- resume fast-forward: was this exact point already solved? ---
        % Matched on (dAR, maxdef), the store's key, so a row written by a
        % session that searched dAR is found as well.
        if OPT.useStore && OPT.resume
            hit = store.lookup(OPT.runId, [defectAspectRatio_i, maxdef_i]);
            if ~isempty(hit)
                % Re-derive J rather than trusting the stored objective_j:
                % the row may have been written by the OTHER optimizer.
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
                nCacheHits = nCacheHits + 1;
                fprintf(['  eval %3d: oblong=%.4f  maxdef=%.4f  ' ...
                         '[CACHED  fitness=%.4g  %s]\n' ...
                         '         gOM=%.2f kHz (MB %.2f + PE %.2f)\n'], ...
                    hit.eval_idx, oblong_i, maxdef_i, hit.fitness, ...
                    hit.status, hit.gOM*1e-3, hit.gMB*1e-3, hit.gPE*1e-3);

                if isfinite(hit.fitness) && hit.fitness > bestSoFar.fitness
                    bestSoFar.fitness = hit.fitness;
                    bestSoFar.oblong  = oblong_i;
                    bestSoFar.maxdef  = maxdef_i;
                    bestSoFar.defectAspectRatio = defectAspectRatio_i;
                    bestSoFar.gOM = hit.gOM;
                    bestSoFar.gMB = hit.gMB;
                    bestSoFar.gPE = hit.gPE;
                end

                % Replayed points join the history so the convergence plots
                % show the whole search, not just this session's part of it.
                history.eval(end+1, 1)    = hit.eval_idx;
                history.dAR(end+1, 1)     = defectAspectRatio_i;
                history.maxdef(end+1, 1)  = maxdef_i;
                history.oblong(end+1, 1)  = oblong_i;
                history.fitness(end+1, 1) = hit.fitness;
                history.Qmech(end+1, 1)   = hit.Qmech;
                history.J(end+1, 1)       = R.J;

                % Keep the counter ahead of every index already stored.
                evalCount = max(evalCount, hit.eval_idx);
                return;
            end
        end

        evalCount = evalCount + 1;
        nSolved   = nSolved + 1;

        % --- the search variables go into P unchanged ---
        Pe = P;
        Pe.maxdef = maxdef_i;
        Pe.oblong = oblong_i;
        sigCtr  = 1 - maxdef_i;
        tCtr_i  = t_0  * sigCtr^(1 - Pe.oblong);
        wCtr_i  = wm_0 * sigCtr^(1 + Pe.oblong);

        % Unique folder per evaluation -- RunNanobeamFEM never reuses a cache.
        evLoc = sprintf('%seval_%04d_ob%.4f_md%.4f%s', ...
            OPT.rootLoc, evalCount, Pe.oblong, maxdef_i, filesep);

        % Strip any stale fileBase so it is regenerated for THIS design.
        if isfield(Pe, 'fileBase'); Pe = rmfield(Pe, 'fileBase'); end
        if OPT.shortFileBase
            Pe.fileBase = shortFileBase('', Pe.oblong, maxdef_i);
        end

        fprintf(['  eval %3d: oblong=%.4f  maxdef=%.4f  (dAR=%.4f, ' ...
                 't_ctr=%.1f nm, w_ctr=%.1f nm)  -> %s\n'], ...
            evalCount, Pe.oblong, maxdef_i, defectAspectRatio_i, ...
            tCtr_i*1e9, wCtr_i*1e9, evLoc);

        wM_i     = NaN;
        gOM_i    = NaN;
        gMB_i    = NaN;
        gPE_i    = NaN;
        Qmech_i  = NaN;
        Qopt_i   = NaN;
        lambda_i = NaN;
        fit      = NaN;
        fitRaw   = NaN;
        status_i = 'failed';

        try
            if OPT.dryRun
                % The real geometry arithmetic still runs, so a pre-filter
                % that disagreed with the builder would show up as 'failed'.
                createRibCavityGeom(Pe);
                ds = dryRunDs(Pe, defectAspectRatio_i, OPT);
            else
                [ds, ~] = RunNanobeamFEM(Pe, evLoc);
            end

            % RunNanobeamFEM catches its own errors and returns an empty ds.
            if isempty(ds) || ~isstruct(ds)
                error('optimize_rib_oblong_maxdef:emptyDs', ...
                    'RunNanobeamFEM returned empty ds');
            end

            % --- check for localised mechanical modes ---
            hasMech = isfield(ds, 'mfem') && isfield(ds.mfem, 'freqs') && ...
                      ~isempty(ds.mfem.freqs);
            if ~hasMech
                warning('optimize_rib_oblong_maxdef:noMechModes', ...
                    'eval %d (oblong=%.4f, maxdef=%.4f): no mechanical modes found', ...
                    evalCount, Pe.oblong, maxdef_i);
                status_i = 'no_mech_modes';
            else
                status_i = 'ok';
            end

            % --- mechanical frequency (prefer the coupled-mode result) ---
            if isfield(ds, 'cpl') && isfield(ds.cpl, 'mechFreq') && ...
                    ~isempty(ds.cpl.mechFreq)
                wM_i = max(real(ds.cpl.mechFreq));
            elseif hasMech
                wM_i = max(real(ds.mfem.freqs));
            end

            % --- check for optical modes ---
            if Pe.solveOpt
                hasOpt = isfield(ds, 'ofem') && isfield(ds.ofem, 'freqs') && ...
                         ~isempty(ds.ofem.freqs);
                if ~hasOpt && ~OPT.dryRun
                    warning('optimize_rib_oblong_maxdef:noOptModes', ...
                        'eval %d (oblong=%.4f, maxdef=%.4f): no optical modes found', ...
                        evalCount, Pe.oblong, maxdef_i);
                    if strcmp(status_i, 'ok')
                        status_i = 'no_opt_modes';
                    else
                        status_i = [status_i, '+no_opt_modes'];
                    end
                end
            end

            % --- g_OM and its mechanism breakdown (signed parts) ---
            [gOM_i, gMB_i, gPE_i] = local_getGOM(ds);

            % --- mechanical Q estimate (diagnostic; the gate if enabled) ---
            [Qmech_i, ~] = local_getQmechRib(ds);

            % --- optical Q and wavelength ---
            [Qopt_i, lambda_i] = local_getQopt(ds, OPT.targetLambda);

            % Both the gated and the ungated score are computed: Nelder-Mead
            % uses the first, the GP surrogate the second.
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
            warning('optimize_rib_oblong_maxdef:evalFailed', ...
                'eval %d (oblong=%.4f, maxdef=%.4f) failed: %s', ...
                evalCount, Pe.oblong, maxdef_i, ME.message);
            status_i = 'failed';
        end

        % --- flag designs rejected by the binary Q_mech gate ---
        if isfinite(fit) && fit == 0 && strcmp(status_i, 'ok')
            status_i = 'Qmech_gated';
        end

        fprintf(['         wM=%.3f GHz   gOM=%.2f kHz (MB %.2f + PE %.2f)\n' ...
                 '         Qm(est)=%.3g   Qo=%.3g   lam=%.1f nm   ' ...
                 'fitness=%.4g   [%s]\n'], ...
            wM_i*1e-9, gOM_i*1e-3, gMB_i*1e-3, gPE_i*1e-3, ...
            Qmech_i, Qopt_i, lambda_i, fit, status_i);

        % Append to CSV immediately so the log survives a mid-run crash.
        % LSiV is not computed for the rib (P.calcS = 0); the column is kept
        % so the file reads like the hole optimizer's.
        fidLog = fopen(logFile, 'a');
        fprintf(fidLog, ['%d,%.6f,%.6f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,', ...
                         '%.6g,%.6g,%.6g,%s,%.6f,%.2f,%.2f\n'], ...
            evalCount, Pe.oblong, maxdef_i, ...
            wM_i*1e-9, gOM_i*1e-3, gMB_i*1e-3, gPE_i*1e-3, NaN, ...
            Qmech_i, Qopt_i, lambda_i, fit, status_i, ...
            defectAspectRatio_i, tCtr_i*1e9, wCtr_i*1e9);
        fclose(fidLog);

        % --- form the negated, log-scaled objective ---
        %   fit  > 0   -> usable design
        %   fit == 0   -> solved, but rejected by the Q_mech gate
        %   fit  = NaN -> no usable mode / solve failed
        if isfinite(fit) && fit > 0
            J = -log10(fit);
            if fit > bestSoFar.fitness
                bestSoFar.fitness = fit;
                bestSoFar.oblong  = Pe.oblong;
                bestSoFar.maxdef  = maxdef_i;
                bestSoFar.defectAspectRatio = defectAspectRatio_i;
                bestSoFar.gOM = gOM_i;
                bestSoFar.gMB = gMB_i;
                bestSoFar.gPE = gPE_i;
            end
        elseif isfinite(fit)
            J = OPT.Jpenalty;
        else
            J = 2 * OPT.Jpenalty;
        end

        % --- record the search trajectory for the convergence plots ---
        history.eval(end+1, 1)    = evalCount;
        history.dAR(end+1, 1)     = defectAspectRatio_i;
        history.maxdef(end+1, 1)  = maxdef_i;
        history.oblong(end+1, 1)  = Pe.oblong;
        history.fitness(end+1, 1) = fit;
        history.Qmech(end+1, 1)   = Qmech_i;
        history.J(end+1, 1)       = J;

        % --- persist immediately so an interrupted run resumes from here ---
        % rec.x is the STORE's key, [defectAspectRatio, maxdef], not the
        % search vector: OptimStore (shared with the hole optimizer) writes
        % x(1) to its 'dar' column and builds param_key from x. Keeping that
        % convention is what records dar AND oblong for every row and keeps
        % studies stored by the dAR-searching version resumable.
        if OPT.useStore
            try
                store.logEval(OPT.runId, evalCount, struct( ...
                    'x',        [defectAspectRatio_i, maxdef_i], ...
                    'oblong',   Pe.oblong, ...
                    'wM',       wM_i,     'gOM',      gOM_i, ...
                    'gMB',      gMB_i,    'gPE',      gPE_i, ...
                    'LSiV',     NaN,      'Qmech',    Qmech_i, ...
                    'Qopt',     Qopt_i,   'lambdaNm', lambda_i, ...
                    'fitness',  fit,      'fitnessRaw', fitRaw, ...
                    'J',        J,        'status',     status_i, ...
                    'evLoc',    evLoc));
            catch ME
                % Persistence must never kill a study that is already
                % several COMSOL-hours deep.
                warning('optimize_rib_oblong_maxdef:storeFailed', ...
                    'Could not persist eval %d: %s', evalCount, ME.message);
            end
        end

        R.fit    = fit;
        R.fitRaw = fitRaw;
        R.Qmech  = Qmech_i;
        R.J      = J;
    end

end  % function optimize_rib_oblong_maxdef

%% ===================== LOCAL HELPERS =====================

function out = mergeStruct(base, over, where)
%MERGESTRUCT  Copy fields of 'over' into 'base'; nested structs merge field
%   by field. A field that 'base' does not have is an error: a mistyped
%   option would otherwise be silently ignored.
out = base;
fn = fieldnames(over);
for k = 1:numel(fn)
    f = fn{k};
    if ~isfield(base, f)
        error('optimize_rib_oblong_maxdef:unknownOption', ...
            'Unknown option %s.%s. Valid fields: %s.', ...
            where, f, strjoin(fieldnames(base), ', '));
    end
    if isstruct(base.(f)) && isstruct(over.(f)) && isscalar(over.(f))
        out.(f) = mergeStruct(base.(f), over.(f), [where, '.', f]);
    else
        out.(f) = over.(f);
    end
end
end

function s = shortFileBase(prefix, oblong, maxdef)
%SHORTFILEBASE  Short per-evaluation file base (see OPT.shortFileBase).
s = sprintf('rib_%sob%.4f_md%.4f', prefix, oblong, maxdef);
end

function dAR = dARofOblong(oblong, maxdef, mirrorAspectRatio)
%DAROFOBLONG  Centre-cell aspect ratio w/t of a hole-scheme rib taper.
%   dAR = (w_mirror/t) * (1 - maxdef)^(2*oblong). Element-wise; the inverse
%   of oblongOfDAR. This is isFabricableRib's first input and x(1) of the
%   run store's key.
dAR = mirrorAspectRatio .* (1 - maxdef).^(2*oblong);
end

function oblong = oblongOfDAR(dAR, maxdef, mirrorAspectRatio)
%OBLONGOFDAR  Taper exponent that gives the centre-cell aspect ratio dAR.
%   oblong = log(dAR/(w_mirror/t)) / (2*log(1 - maxdef)). Element-wise; used
%   to bring stored (dAR, maxdef) points into the search variables.
oblong = log(dAR ./ mirrorAspectRatio) ./ (2*log(1 - maxdef));
end

function ds = dryRunDs(Pe, defectAspectRatio, OPT)
%DRYRUNDS  Analytic stand-in for RunNanobeamFEM (OPT.dryRun). NOT physics:
%   a smooth peak at (aspect ratio 8, maxdef 0.20), which for the default
%   mirror cell (w/t = 1300/280) is oblong = -1.2192, with the same fields,
%   units and orders of magnitude the real ds has, so every code path of the
%   optimizer can be exercised without COMSOL. Q_mech drops below 1e7 away
%   from maxdef = 0.19, so the optional gate has something to reject.
u  = log(defectAspectRatio/8)/0.6;
v  = (Pe.maxdef - 0.20)/0.07;
fM = OPT.targetFreq*(1 + 0.5*(0.20 - Pe.maxdef));
Qm = 3e7*exp(-((Pe.maxdef - 0.19)/0.09)^2);
wM = 2*pi*fM;
g  = 2e5*exp(-(u^2 + v^2)/2);

ds = struct();
ds.P = Pe;
% COMSOL eigenvalue convention: lambda = -i*omega, Q = |Im|/(2|Re|).
ds.mfem.PVals    = -wM/(2*Qm) + 1i*wM;
ds.mfem.freqs    = fM;
ds.mfem.locInd   = 1;
ds.mfem.locRatio = 0.8;
ds.ofem.Q         = 1e6*exp(-u^2/2);
ds.ofem.lambda    = OPT.targetLambda*1e-9*(1 + 0.3*(Pe.maxdef - 0.20));
ds.ofem.QAll      = ds.ofem.Q;
ds.ofem.lambdaAll = ds.ofem.lambda;
ds.ofem.freqs     = 299792458/ds.ofem.lambda;
ds.cpl.gMax     = g;
ds.cpl.gMBmax   = 1.4*g;
ds.cpl.gPEmax   = -0.4*g;
ds.cpl.mechFreq = fM;
ds.cpl.mSol     = 1;
ds.cpl.optWvl   = ds.ofem.lambda;
end

%% ===================== LOCAL FITNESS HELPERS =====================
%
% ACTIVE fitness  : fitness_ribCoop  (wired up via OPT.fitnessFcn above).
% ALTERNATIVE     : fitness_maxGOM. To switch, edit the OPT.fitnessFcn and
%                   OPT.fitnessFcnRaw lines in the knobs section, e.g.
%                       OPT.fitnessFcn    = @(ds) fitness_maxGOM(ds);
%                       OPT.fitnessFcnRaw = OPT.fitnessFcn;
%
% Anything returning NaN / non-finite is treated as an unusable point;
% anything returning 0 as a solved design rejected by the Q_mech gate.
%
% local_getQopt and local_getGOM are copies of the functions of the same
% name in optimize_oblong_maxdef.m; keep them in step.

function fit = fitness_ribCoop(ds, targetLambdaNm, tolNm, QmechMin, fWindow, fTol)
%FITNESS_RIBCOOP  Cooperativity-like figure of merit for the rib cavity.  [ACTIVE]
%
%       fit = g_OM^2 * Q_opt * L(lambda_opt) * W(f_mech)
%
%   L(lambda) = exp(-((lambda - targetLambda)/tolNm)^2), a Gaussian, so tolNm
%   is the 1/e half-width.
%   W(f) = 1 when fWindow is empty or f lies inside [fWindow(1) fWindow(2)],
%   else exp(-(df/fTol)^2) with df the distance outside the window. f is the
%   mechanical mode of the best g_OM pair.
%
%   Q_mech gate (the hole optimizer's fitness_coopProduct has it always on):
%     QmechMin = -Inf       -> no gate; Q_mech is not needed at all
%     Q_mech <  QmechMin    -> fit = 0    (design rejected outright)
%     Q_mech unavailable    -> fit = NaN  (cannot be judged)
%   Q_mech comes from local_getQmechRib: a PML value if ds.mfem.QAll exists,
%   otherwise an estimate from the low-reflecting beam end.
%
%   Unlike fitness_coopProduct there is no partial score: without g_OM, an
%   optical Q and a wavelength the result is NaN, so every finite value is
%   on the same scale. A usable design never scores exactly 0 (the result is
%   floored at realmin), so 0 always means "gated".
fit = NaN;
[Qopt, lambdaNm] = local_getQopt(ds, targetLambdaNm);
gOM = local_getGOM(ds);
if isnan(gOM) || gOM <= 0 || isnan(Qopt) || isnan(lambdaNm)
    return;
end
[Qmech, fMech] = local_getQmechRib(ds);

% --- binary Q_mech acceptance gate (off for QmechMin = -Inf) ---
if QmechMin > -Inf
    if isnan(Qmech); return; end
    if Qmech < QmechMin
        fit = 0;
        return;
    end
end

score = gOM^2 * Qopt * exp(-((lambdaNm - targetLambdaNm) / tolNm)^2);

% --- mechanical-frequency window ---
if ~isempty(fWindow)
    if isnan(fMech); return; end
    dF = max([fWindow(1) - fMech, fMech - fWindow(2), 0]);
    score = score * exp(-(dF / fTol)^2);
end

fit = max(score, realmin);
end

function [Qmech, fMech, src] = local_getQmechRib(ds)
%LOCAL_GETQMECHRIB  Mechanical Q of the cavity mode, PML or estimated.
%
%   Mode: the mechanical mode of the best g_OM pair (ds.cpl.mSol, an
%   absolute solution number, CalcGOM.m:327); without it, the localized mode
%   with the highest localization ratio.
%   Q:    ds.mfem.QAll when it exists (src 'pml': a mechanical PML was
%         solved; not possible for the rib today). Otherwise from the complex
%         eigenvalue, Q = |Im|/(2|Re|) (src 'endBC'), the same expression
%         SolveNanobeamFEM.m:397 uses. With P.mEndBC = 'lowreflecting' the
%         real part is the energy leaving through the beam end; a lossless
%         eigenvalue (fixed end) gives NaN, not Inf.
Qmech = NaN;  fMech = NaN;  src = 'none';
if ~isfield(ds, 'mfem') || ~isfield(ds.mfem, 'freqs') || isempty(ds.mfem.freqs)
    return;
end
mfem   = ds.mfem;
nModes = numel(mfem.freqs);

idx = [];
if isfield(ds, 'cpl') && isfield(ds.cpl, 'mSol') && isscalar(ds.cpl.mSol) && ...
        ds.cpl.mSol >= 1 && ds.cpl.mSol <= nModes && ds.cpl.mSol == round(ds.cpl.mSol)
    idx = ds.cpl.mSol;
elseif isfield(mfem, 'locInd') && ~isempty(mfem.locInd) && ...
        isfield(mfem, 'locRatio') && ~isempty(mfem.locRatio) && ...
        max(mfem.locInd) <= min(numel(mfem.locRatio), nModes)
    [~, k] = max(mfem.locRatio(mfem.locInd));
    idx    = mfem.locInd(k);
end
if isempty(idx); return; end

fMech = real(mfem.freqs(idx));
if isfield(mfem, 'QAll') && numel(mfem.QAll) >= idx
    Qmech = mfem.QAll(idx);
    src   = 'pml';
elseif isfield(mfem, 'PVals') && numel(mfem.PVals) >= idx
    ev = mfem.PVals(idx);
    if real(ev) ~= 0 && imag(ev) ~= 0
        Qmech = abs(imag(ev)) / (2*abs(real(ev)));
        src   = 'endBC';
    end
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
%   gOM - max |g_OM| in Hz. NaN if not computed (P.calcG = 0).
%   gMB - moving-boundary part, Hz, SIGNED.
%   gPE - photoelastic part, Hz, SIGNED.
%   The parts are recorded by CalcGOM at the same mode pair as the maximum
%   and returned signed on purpose: they routinely partially cancel.
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

function s = cplSummary(gOM, gMB, gPE)
%CPLSUMMARY  Coupling of one design as a struct, all in Hz.
%   gOM is the net magnitude; gMB (moving boundary) and gPE (photoelastic)
%   are the SIGNED parts of the same mode pair, so gMB + gPE = +/- gOM.
%   cancelRatio = (|gMB| + |gPE|)/|gOM| as in CalcGOM: 1 when the two
%   mechanisms reinforce, above 1 when they oppose. NaN where not computed.
s = struct('gOM', gOM, 'gMB', gMB, 'gPE', gPE, ...
           'cancelRatio', (abs(gMB) + abs(gPE)) / max(abs(gOM), eps));
end

function txt = cplText(s)
%CPLTEXT  One-line console form of a cplSummary struct, in kHz.
txt = sprintf('%.2f kHz (MB %.2f + PE %.2f), cancelRatio %.2f', ...
    s.gOM*1e-3, s.gMB*1e-3, s.gPE*1e-3, s.cancelRatio);
if s.cancelRatio > 1.5
    txt = [txt, '  [MB and PE oppose]'];
end
end

function fit = fitness_maxGOM(ds)
%FITNESS_MAXGOM  Maximize bare optomechanical coupling.  [ALTERNATIVE]
%   Returns max |g_OM| in Hz, NaN if g_OM was not computed. Ignores the
%   optical Q and wavelength entirely -- a quick "where is the coupling
%   strongest" scan only.
fit = local_getGOM(ds);
if ~(fit > 0)
    fit = NaN;
end
end
