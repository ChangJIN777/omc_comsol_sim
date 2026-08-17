%% BAYESOPT_BOOMERANG
% Bayesian optimization of the boomerang unit cell for a complete mechanical
% bandgap centred on a target frequency.
%
% This script is the Bayesian-optimization counterpart to two existing
% scripts in this directory:
%   * optimization_client.m     - supplied the objective (gap ratio weighted by
%                                 a Gaussian penalty on the mid-gap frequency)
%                                 but drove it with fminsearch.
%   * sweep_boomerang_code.m    - supplied the canonical boomerang P-struct and
%                                 the parameter ranges used for the grid sweeps.
%
% Why bayesopt instead of fminsearch: each objective evaluation is a COMSOL
% eigenfrequency solve costing minutes, so the evaluation budget - not the
% convergence rate - is the binding constraint. bayesopt builds a Gaussian
% process surrogate and spends each solve where it is most informative, it
% handles integer-valued variables natively, it prunes fabrication-infeasible
% candidates *before* spending a solve on them (XConstraintFcn), and it can be
% resumed after an interruption. fminsearch does none of these and, being a
% local simplex method, is additionally sensitive to the flat/zero regions this
% objective produces wherever no complete gap exists.
%
% DESIGN VARIABLES (all in integer nanometres - see note on caching below)
%   a   lattice constant; side of the rhombic cell. Sets the cell footprint.
%   r   hole ARM LENGTH, cell centre to arm tip. Not a cell dimension.
%   w   hole ARM WIDTH - the narrowest etched feature. Not a cell dimension.
%   th  full slab thickness in z
% Fillet radii r1/r2 are held fixed, matching sweep_boomerang_code.m.
%
% The cell is the primitive rhombus of a hexagonal lattice (side a), and the
% hole is a three-pointed star: three arms, each w wide and r long, radiating
% from the cell centre at 120 deg spacing. See test_Boomerang.m for the full
% geometry write-up. Beware that "unit cell width/height" comments elsewhere
% in this directory mislabel w and r - they are hole dimensions, and `a` alone
% sets the cell footprint.
%
% Variables are declared as integers in nm rather than continuous metres on
% purpose. solveBands names its output files by rounding every dimension to
% whole nm ('%.0f') and short-circuits the solve when that file already
% exists. Continuous variables would therefore alias many distinct GP query
% points onto one cached result, silently feeding the surrogate stale data.
% On an integer-nm grid the rounding is exact, so a repeated point is a true
% cache hit and existing sweep data under datLoc is reused for free.
%
% OBJECTIVE (minimized)
%   objective = -(gapSize/midGap) * exp(-((f_target - midGap)/sigma)^2)
% evaluated over every complete gap returned by findGaps, keeping the best.
% Zero means "no usable complete gap", which is the worst attainable value
% since the fitness is non-negative by construction.
%
% SOLVER BACKEND / DRY-RUN MODE
% cfg.solverBackend selects what actually produces the band structure:
%   'comsol'     the real COMSOL eigenfrequency solve. Minutes per evaluation,
%                and the ONLY backend whose numbers mean anything. Default.
%   'surrogate'  a closed-form ANALYTIC fake band structure, microseconds per
%                evaluation (surrogateBoomerangBands). It exists to debug the
%                LOOP - the objective plumbing, XConstraintFcn, the iteration
%                log, the checkpoint and all four summary figures get exercised
%                end to end in seconds instead of hours. It is not physics.
%   'stub'       no bands at all: every evaluation reports a solver failure
%                (NaN). That is the other path worth testing - the skip/failure
%                branch the success path never touches.
% Anything else is an error, never a silent fall-through: a typo'd backend name
% must not be able to start a multi-hour real study, nor to silently fake one.
%
% The string-valued backend (rather than a bare boolean dryRun flag) follows
% python-scripts/src/objective.py, which selects with mech_backend='comsol' vs
% 'surrogate_stub' and raises ValueError on anything else. A string scales past
% two tiers and self-documents wherever it is printed, logged or saved; 'stub'
% is this file's shorter name for that file's 'surrogate_stub'.
%
% CACHE ISOLATION - the one hazard that matters here.
% solveBands caches by FILENAME: it skips the solve entirely whenever
% [datLoc fileBase '_bds.mat'] exists, and boomerangParams reproduces the sweep
% scripts' filename byte-for-byte precisely so that real sweep data is reused.
% A dry run that wrote those same filenames would therefore make every later
% REAL run at those geometries silently load fabricated band data - invisible,
% permanent corruption of the one thing this directory exists to produce. Four
% independent guards prevent it, any one of which would be sufficient:
%   (1) a separate output folder, DRYRUN_<backend>_<date>/, so a real run never
%       even looks in the directory the synthetic files live in;
%   (2) a DRYRUN filename prefix on every synthetic _bds.mat, so the two could
%       not collide even if the folders were merged by hand;
%   (3) a provenance marker (ds.isSynthetic) inside every synthetic .mat, which
%       assertBandDataProvenance checks on EVERY cache load - a real run that is
%       somehow handed synthetic data errors out instead of using it;
%   (4) cfg.dryRunSaveBands = 0 writes no .mat at all, for the paranoid.
% The same separation applies to the iteration log, the checkpoint, the results
% .mat and every figure, all of which are named through cfg.filePrefix.
%
% Synthetic results are additionally labelled everywhere a human might meet
% them: the console banner, a 'backend' column in the log, a DRYRUN_README.txt
% dropped in the output folder, the figure names, the figure headlines and the
% summary text panel.
%
% REQUIREMENTS
%   * Statistics and Machine Learning Toolbox (bayesopt, optimizableVariable)
%   * COMSOL with LiveLink for MATLAB, already running, before this script is
%     started - solveBands -> runBands_2D calls the COMSOL Java API.
%
% USAGE
%   >> run('bayesopt_boomerang.m')
% then, to spend more solves on the same surrogate:
%   >> results = resume(results,'MaxObjectiveEvaluations',20);
% To debug the loop without COMSOL, set cfg.solverBackend = 'surrogate' (or
% 'stub') in the configuration block below and run exactly the same way.
%
% OUTPUTS written under cfg.datLoc
%   *_bds.mat                    per-evaluation band structure (via solveBands)
%   bayesopt_boomerang_log.txt   one tab-separated row per evaluation
%   bayesopt_boomerang_state.mat checkpoint of the BayesianOptimization object,
%                                rewritten after every iteration
%   bayesopt_boomerang_results.mat   final results + cfg + best design
%   bayesopt_boomerang_summary.*     composite review figure (see below)
%   bayesopt_boomerang_geometry.*    best unit cell, top view
%   bayesopt_boomerang_bestbands.*   band structure of the best design
%   bayesopt_boomerang_progress.*    convergence + design-space slices
% In a dry run every one of these moves to a DRYRUN_<backend>_<date>/ folder and
% gains a _DRYRUN infix (cfg.filePrefix), plus a DRYRUN_README.txt beside them.
% Each figure is written as both .png and .fig, matching what solveBands does
% with its per-evaluation band diagrams.
%
% The figure section is modelled on _run_characterization in
% python-scripts/scripts/run_opt_comsol.py - one composite multi-panel figure
% carrying the headline numbers in its title, plus standalone single-topic
% figures beside it. It reads only data already on disk (the cached _bds.mat of
% the winning design) and never calls COMSOL, so it still works if the LiveLink
% connection has since dropped.

clear all; clc; close all;                                                  %#ok<CLALL>

%% ------------------------------------------------------------------------
%  Configuration
%  ------------------------------------------------------------------------
cfg = struct();

% --- solver backend / dry run --------------------------------------------
% What actually produces the band structure for one design point. A STRING
% rather than a boolean dryRun flag, mirroring mech_backend in
% python-scripts/src/objective.py: three tiers do not fit in a boolean, a fourth
% would not either, and the name self-documents everywhere it is printed, logged
% or saved rather than becoming an anonymous "1".
%
%   'comsol'     the real solve. Minutes per evaluation, and the only backend
%                whose numbers mean anything.
%   'surrogate'  analytic FAKE bands (surrogateBoomerangBands), microseconds per
%                evaluation. Correctly shaped and plausibly trending, physically
%                meaningless. Exercises everything downstream of the solve:
%                gapFitness, the log, the checkpoint, all four summary figures.
%   'stub'       no bands at all; every evaluation reports a solver failure
%                (NaN). Exercises the failure/skip path, which the surrogate -
%                by succeeding - never reaches. objective.py calls this tier
%                'surrogate_stub'.
%
% Read the CACHE ISOLATION note in the header block before using either cheap
% backend. Short version: dry-run output lives in its own folder, under its own
% filename prefix, carrying an in-file provenance marker, so no real run can
% ever load it.
% cfg.solverBackend = 'comsol';
cfg.solverBackend = 'surrogate';

% Fail loud on an unknown backend BEFORE anything is derived from the name -
% the same intent as `raise ValueError(f"unknown mech_backend {mech_backend}")`
% in objective.py. Without this, a typo would fall through the strcmp below to
% isDryRun = true and quietly produce a whole folder of fabricated results;
% mistype it the other way and you would be waiting on real COMSOL solves you
% did not ask for. Neither failure is acceptable, so it is an error, here, first.
assertKnownBackend(cfg.solverBackend);

% Single source of truth for "is this a dry run". Everything downstream tests
% this flag instead of re-comparing strings, so there is exactly one place where
% the real backend is distinguished from the cheap ones.
cfg.isDryRun = ~strcmp(cfg.solverBackend,'comsol');

% --- dry-run options (all ignored unless cfg.isDryRun) -------------------
cfg.dryRunDelay = 0;        % artificial seconds per evaluation. 0 = as fast as
                            % MATLAB can go, which is the point of the mode.
                            % 0.5-2 s is useful when the thing being debugged is
                            % the live bayesopt plots themselves, which are hard
                            % to watch fill in at 10 000 evaluations a minute.
cfg.dryRunFailEvery = 0;    % 0 = off. N > 0 makes any design whose
                            % a+r+w+th (in nm) is divisible by N report a
                            % simulated solver failure, so the NaN branch of
                            % boomerangObjective, the 'solver-error' log status
                            % and the red error markers on the convergence panel
                            % all get exercised. Keyed on the DESIGN, never on
                            % an evaluation counter: bayesopt is told
                            % 'IsObjectiveDeterministic',true, and a counter
                            % would make a repeated design point return a
                            % different value, violating that contract.
cfg.dryRunSaveBands = 1;    % 1 writes the synthetic ds to a _bds.mat under the
                            % dry-run folder. That is what lets a dry run also
                            % exercise the cache-hit branch and
                            % loadCachedBandData, i.e. the code paths that read
                            % band data back off disk. Set 0 to touch nothing
                            % but the log; the visualization then degrades to
                            % its 'band-data-unavailable' path, which is itself
                            % worth testing.
cfg.dryRunPrefname = 'DRYRUN';  % filename prefix stamped onto synthetic band
                            % data. See boomerangParams for why it is applied by
                            % hand rather than left to solveBands' P.prefname.

% --- objective definition -------------------------------------------------
cfg.targetFreq = 7e9;      % target mechanical mid-gap frequency [Hz]
cfg.sigma      = 3e9;       % width of the Gaussian frequency penalty [Hz]

% --- fabrication tolerance ------------------------------------------------
% Minimum feature size, applied to the hole arm width w - the narrowest
% etched feature in the tri-arm geometry. See boomerangFabConstraint for why
% the `a - r` term from boomerang_optimize_sweep_diamond.m is deliberately not
% reproduced. The in-cell guard is no longer the analytic r < sqrt(3)*a/4 -
% that assumed a centred hole; it is now the measured armsOverhang test.
% Held in integer nm so the comparison is exact. Storing it in metres and
% comparing w*1e-9 >= 50e-9 would be a round-off coin flip on designs sitting
% exactly on the limit, e.g. w = 50 nm, because the computed product and the
% parsed literal need not agree to the last bit.
cfg.minFeatureNm = 50;      % [nm]  minimum ETCHED linewidth (arm width w)

% The other half of the process window: the minimum SOLID feature, i.e. the
% thinnest dielectric wall left standing. That wall runs between this cell's
% hole and a neighbouring cell's, so it exists only once the lattice is applied
% and cannot be read off a, r and w - calcFillingFactor measures it as the
% closest approach between the hole boundary and its six nearest translates.
%
% It matters because it is NOT redundant with minFeatureNm. At the current
% test_Boomerang geometry the wall is 148 nm against w = 125 nm, so w binds -
% but the wall shrinks fast as r grows toward the containment limit and can
% drop below w well before an arm actually crosses the cell edge.
%
% Set 0 to skip the check, which also skips the expensive measurement.
cfg.minSolidFeatureNm = 50;    % [nm]

% Extra clearance demanded between the hole and the cell boundary, on top of
% the exact non-intersection test in boomerangFabConstraint. 0 means "must not
% cross", which the armsOverhang test already guarantees; a positive value
% demands the hole stay that far clear. This is the constraint that ties r and
% w together - the farthest point of the hole is at hypot(w/2, r) from the hole
% centre, so a wider arm costs edge clearance exactly as a longer one does.
% Set 0 to skip; combined with cfg.minSolidFeatureNm = 0 that skips the
% expensive geometry pass entirely.
cfg.holeEdgeMarginNm = 0;      % [nm]

% --- search space, in integer nm (ranges from sweep_boomerang_code.m) -----
cfg.bounds.a  = [600 1000];
cfg.bounds.r  = [100  250];   % sweep used r in [0.6*250, 250] nm
cfg.bounds.w  = [ 50  120];

% --- fixed geometry ------------------------------------------------------
% th is FIXED, not searched. It used to be a fourth optimizable variable over
% cfg.bounds.th = [275 350]; holding it constant drops the search to three
% dimensions, which combined with the filling-factor constraint below leaves
% effectively two free directions - a much easier surface for a GP to learn on
% the evaluation budget available. Fix it at whatever thickness the process
% actually delivers; there is no point optimizing a dimension you cannot set.
cfg.th = 300e-9;            % full slab thickness in z [m]

% --- filling-factor constraint -------------------------------------------
% Restrict the search to designs whose in-plane air/dielectric area ratio sits
% in a band around cfg.fillingFactor. Measured by calcFillingFactor from the
% same polygons buildBoomerangUnitCell hands COMSOL, so the constraint and the
% geometry cannot disagree.
%
% An exact equality is unusable: (a,r,w) live on an integer-nm grid, so the
% set where the ratio hits a value exactly has measure zero and every candidate
% would be pruned. The tolerance turns the equality into a thin feasible shell,
% two-dimensional in a three-dimensional box.
%
% Set cfg.fillingFactor = [] to disable and recover the previous behaviour.
% Choose the target WITH the bounds in mind - the startup feasibility scan
% below reports what fraction of the box survives, and refuses to start a study
% that has almost nothing left to search.
cfg.fillingFactor    = 0.49;   % target area(air)/area(dielectric), or []
cfg.fillingFactorTol = 0.1;   % accept |ff - target| <= this
cfg.fillingFactorMinFeasible = 0.02;   % refuse to start below this fraction
cfg.r1 = 10e-9;             % fillet radius at the INNER corners, where the
                            % three arms meet near the cell centre
cfg.r2 = 10e-9;             % fillet radius at the OUTER arm tips

% --- solver fidelity (kept identical to sweep_boomerang_code.m) ----------
cfg.kpts     = 9;           % k-points EXCLUDING gamma
cfg.nbands   = 10;
cfg.meshSize = 4;
cfg.maxDof   = 3e6;

% --- optimizer budget ---------------------------------------------------
cfg.maxEvaluations = 40;    % total COMSOL solves bayesopt may spend
cfg.numSeedPoints  = 8;     % random design points before the GP takes over

% --- I/O ---------------------------------------------------------------
% Paths are built with fullfile, never by pasting a separator in by hand, so the
% same script produces a real nested directory on Windows, macOS and Linux
% alike. This script was the first to do so; the rest of the repository used to
% hard-code '.\test\...' and its output was therefore only correct on Windows:
% on macOS and Linux a backslash is an ordinary filename character, not a
% separator, so those scripts created a single file literally called
% "test\boomerang_sweep\08132026\..." instead of a folder tree - a name Windows
% cannot even represent, which is what broke `git pull` on the Windows machine.
% Those call sites have since been converted to fullfile as well; keep it that
% way. fullfile emits filesep for the host platform and
% collapses duplicate separators, which is what makes the outputs readable
% wherever they are opened.
%
% The one deliberate exception is the TRAILING separator on cfg.datLoc, which is
% kept because solveBands needs it: it normalises a separator into a local
% variable but then writes the .mat with the raw P.datLoc, so without one the
% band data lands in the parent folder under a mangled name. filesep is used for
% that too rather than a literal '\'.
currentDate = datestr(now,'mmddyyyy');                                      %#ok<TNOW1,DATST>

% Guard (1) of the four cache-isolation guards listed in the header block: a dry
% run gets its OWN dated folder, named for the backend that produced it. A real
% run therefore never looks inside the directory the synthetic _bds.mat files
% live in, which is what makes it structurally impossible for solveBands' filename
% cache to serve fabricated bands to a real study - the two namespaces do not
% intersect at all, rather than merely being unlikely to collide.
%
% cfg.filePrefix exists so the per-study artifacts follow the same switch. For
% 'comsol' it is exactly the literal that used to be hard-coded at each use
% site, so a real study writes byte-identical filenames to before this option
% existed. That equality is deliberate and worth preserving if these lines are
% ever edited: it is what makes "dry run off" mean "unchanged".
if cfg.isDryRun
    datedFolder    = ['DRYRUN_',cfg.solverBackend,'_',currentDate];
    cfg.filePrefix = 'bayesopt_boomerang_DRYRUN';
else
    datedFolder    = currentDate;
    cfg.filePrefix = 'bayesopt_boomerang';
end
% Trailing filesep appended deliberately - see the note above on solveBands.
cfg.datLoc    = [fullfile('.','test','boomerang_bayesopt',datedFolder),filesep];
cfg.logPath   = fullfile(cfg.datLoc,[cfg.filePrefix,'_log.txt']);
cfg.statePath = fullfile(cfg.datLoc,[cfg.filePrefix,'_state.mat']);

cfg.plotgeom          = 0;  % 1 to plot geometry each evaluation (slow, noisy)
cfg.savebndplot       = 1;  % 1 to save a band diagram per evaluation
cfg.closeSolveFigures = 1;  % close figures opened by the solve, keeping the
                            % bayesopt live plots alive

% --- post-run summary figures --------------------------------------------
% Drawn once, after bayesopt returns, from data already on disk. Behind a flag
% so an unattended/headless run can skip them, and so a plotting problem can be
% switched off without touching the study itself.
cfg.makeSummaryFigures  = 1;    % 0 to skip the post-run visualization entirely
cfg.figResolution       = 150;  % PNG export resolution [dpi]. 150 matches the
                                % house style of run_opt_comsol.py: readable
                                % when pasted into a report, not enormous.
cfg.figNPeriods         = 3;    % lattice periods tiled in the geometry panel.
                                % A single cell hides the inter-hole ligament -
                                % the feature that actually sets the gap - so
                                % the default matches the n_periods=3 used by
                                % _draw_geometry in run_opt_comsol.py.
cfg.closeSummaryFigures = 0;    % 0 leaves the saved figures on screen, which is
                                % what you want after babysitting a multi-hour
                                % run. Set 1 for unattended runs; the close is
                                % then done on the handles created here, never
                                % with a bare "close all" - that would also
                                % destroy the bayesopt live plots.

if ~exist(cfg.datLoc,'dir')
    mkdir(cfg.datLoc);
end

% Drop a plain-text warning next to the synthetic data and shout at the console.
% Both are for the same reader: whoever opens this folder in a month, having
% forgotten which run was which. The folder name already says DRYRUN, but a name
% is easy to skim past and a file called DRYRUN_README.txt is not.
if cfg.isDryRun
    writeDryRunMarker(cfg);
    printDryRunBanner(cfg);
end

initIterationLog(cfg.logPath);

%% ------------------------------------------------------------------------
%  Filling-factor feasibility scan
%  ------------------------------------------------------------------------
% A coarse grid over the bounds, scored through the SAME constraint function
% bayesopt will use, to answer one question before any solve time is spent: is
% there anything left to search?
%
% This exists because the failure mode it catches is silent and expensive.
% bayesopt does not report "your constraint pruned everything" - it keeps
% drawing candidates, finds almost none admissible, and either stalls or spends
% the whole budget in a sliver of the box. The filling-factor band makes that
% easy to hit by accident: the ratio spans roughly 0.03 to 0.37 over the shipped
% bounds with a median near 0.09, so a target of 0.25 leaves under 1% of the box
% at +/- 0.01 while a target of 0.15 leaves about 7%.
%
% Coarse on purpose. This is a go/no-go check, not a survey - a fine grid would
% cost more than it tells you, since each point runs the polygon reconstruction.
if ~isempty(cfg.fillingFactor)
    nScan = 12;                       % 12^3 = 1728 candidates, seconds
    [aG,rG,wG] = ndgrid( ...
        round(linspace(cfg.bounds.a(1), cfg.bounds.a(2), nScan)), ...
        round(linspace(cfg.bounds.r(1), cfg.bounds.r(2), nScan)), ...
        round(linspace(cfg.bounds.w(1), cfg.bounds.w(2), nScan)));
    scanTbl  = table(aG(:),rG(:),wG(:),'VariableNames',{'a','r','w'});
    scanOK   = boomerangFabConstraint(scanTbl,cfg);
    fracOK   = nnz(scanOK)/numel(scanOK);

    fprintf(['Filling-factor constraint: target %.4f +/- %.4f\n' ...
             '  feasible fraction of the search box: %.2f%% (%d of %d scanned)\n'], ...
        cfg.fillingFactor,cfg.fillingFactorTol,100*fracOK,nnz(scanOK),numel(scanOK));

    if fracOK == 0
        error('bayesopt_boomerang:fillingFactorInfeasible', ...
            ['No point on a %d^3 scan of the bounds satisfies the filling-' ...
             'factor band %.4f +/- %.4f, so bayesopt would have nothing to ' ...
             'search and the study was NOT started.\nEither widen ' ...
             'cfg.fillingFactorTol, move cfg.fillingFactor toward what the ' ...
             'bounds can actually produce, or widen cfg.bounds.'], ...
            nScan,cfg.fillingFactor,cfg.fillingFactorTol);
    elseif fracOK < cfg.fillingFactorMinFeasible
        % Warn rather than error: a thin shell is workable if it is where the
        % design has to live, it just needs more seed points to find.
        warning('bayesopt_boomerang:fillingFactorTight', ...
            ['Only %.2f%% of the search box satisfies the filling-factor ' ...
             'band, below the %.2f%% comfort threshold. bayesopt can work ' ...
             'here but its seed points are drawn from the box, so most will ' ...
             'be rejected and the GP may start with very little. Consider ' ...
             'raising cfg.numSeedPoints, widening cfg.fillingFactorTol, or ' ...
             'tightening cfg.bounds around the feasible shell.'], ...
            100*fracOK,100*cfg.fillingFactorMinFeasible);
    end
end

%% ------------------------------------------------------------------------
%  Optimization variables
%  ------------------------------------------------------------------------
optVars = [ ...
    optimizableVariable('a', cfg.bounds.a, 'Type','integer'), ...
    optimizableVariable('r', cfg.bounds.r, 'Type','integer'), ...
    optimizableVariable('w', cfg.bounds.w, 'Type','integer')];
% th is deliberately absent - see cfg.th above.

objFcn   = @(x) boomerangObjective(x,cfg);
xConFcn  = @(x) boomerangFabConstraint(x,cfg);
outFcn   = @(res,state) checkpointState(res,state,cfg.statePath);

%% ------------------------------------------------------------------------
%  Run the optimization
%  ------------------------------------------------------------------------
fprintf('Starting Bayesian optimization: %d evaluations, target %.2f GHz\n', ...
    cfg.maxEvaluations, cfg.targetFreq/1e9);
fprintf('Results -> %s\n',cfg.datLoc);

results = bayesopt(objFcn, optVars, ...
    'MaxObjectiveEvaluations',  cfg.maxEvaluations, ...
    'NumSeedPoints',            cfg.numSeedPoints, ...
    'AcquisitionFunctionName',  'expected-improvement-plus', ...
    'IsObjectiveDeterministic', true, ...   % a COMSOL solve is repeatable
    'XConstraintFcn',           xConFcn, ...
    'OutputFcn',                outFcn, ...
    'PlotFcn',                  {@plotMinObjective,@plotObjective}, ...
    'Verbose',                  1);

%% ------------------------------------------------------------------------
%  Report
%  ------------------------------------------------------------------------
xBest = results.XAtMinObjective;
fprintf('\n===== Best boomerang unit cell =====\n');
fprintf('  a  = %d nm\n', xBest.a);
fprintf('  r  = %d nm\n', xBest.r);
fprintf('  w  = %d nm\n', xBest.w);
fprintf('  th = %d nm  (FIXED, not optimized)\n', round(cfg.th*1e9));
fprintf('  objective = %.6f  (fitness = %.6f)\n', ...
    results.MinObjective, -results.MinObjective);
fprintf('  iteration log: %s\n', cfg.logPath);

save(fullfile(cfg.datLoc,[cfg.filePrefix,'_results.mat']),'results','cfg','xBest');

% Repeated at the END as well as the start. The start-of-run banner has scrolled
% far off the top of the console by now, and "best boomerang unit cell" above is
% exactly the block someone would screenshot - so the disclaimer has to be next
% to it, not thousands of lines earlier. cfg is saved inside the .mat too, so the
% backend travels with the results.
if cfg.isDryRun
    printDryRunBanner(cfg);
end

%% ------------------------------------------------------------------------
%  Visualization
%  ------------------------------------------------------------------------
% Deliberately placed AFTER the save above. A study that has already spent
% hours of COMSOL time must not lose its results to a plotting bug, so the
% .mat is on disk before anything is drawn and everything below is
% best-effort: this block is wrapped, each figure is wrapped inside
% visualizeBestDesign, and each individual panel is wrapped again. The worst
% case is a warning and a missing panel, never a lost study.
if cfg.makeSummaryFigures
    try
        visualizeBestDesign(results,xBest,cfg);
    catch ME
        warning('bayesopt_boomerang:visualizationFailed', ...
            'Post-run visualization failed (results were already saved): %s', ...
            ME.message);
    end
end

%% ========================================================================
%  Local functions
%  ========================================================================

function tf = boomerangFabConstraint(x,cfg)
%BOOMERANGFABCONSTRAINT Deterministic feasibility test on the design table.
%
% bayesopt calls this with a table of MANY candidate rows and expects a
% logical column vector, so it must be vectorized - do not assume height 1.
% Points failing here are never passed to the objective, which is the whole
% point: an infeasible geometry costs zero COMSOL time instead of a wasted
% solve plus a penalty value that distorts the surrogate.

% Compared in integer nm, the same units the design variables are declared
% in, so a design exactly on the fabrication limit is not decided by
% floating-point round-off.
%
% Five tests, ordered cheapest and most-selective first. That ordering is
% load-bearing, not tidiness: tests (4) and (5) need calcFillingFactor's O(n^2)
% wall measurement, and running it only on rows that already passed the
% filling-factor band - a few percent of the box - is what makes it affordable
% inside a function bayesopt calls on thousands of candidates per iteration.
%
% (1) Minimum feature. For the tri-arm hole the narrowest feature is simply
%     the arm width w, since that is an etched slit. Note that the `a - r`
%     term used as a fab proxy in boomerang_optimize_sweep_diamond.m is not
%     applied here: it has no geometric meaning for this shape. The diamond
%     ligament between neighbouring holes is set by the lattice and the arm
%     geometry together, and over these bounds it never drops below ~210 nm,
%     while `a - r` reports 350 nm and up - so that term would neither bind
%     nor measure anything real.
%
% (2) Arms must stay inside the cell, or the hole punches through the
%     boundary and the geometry build produces a shape the periodic BCs no
%     longer describe. This USED to be the analytic test r < sqrt(3)*a/4,
%     which is the containment limit for a hole sitting on the cell
%     CENTROID. That stopped being the geometry the builder produces when
%     resolveHoleCentreFrac introduced the default 0.4 diagonal offset, so
%     the test is now the measured calcFillingFactor armsOverhang flag -
%     exact for any hole position, and automatically correct if the offset
%     is retuned again.
%
% (3) Filling factor within cfg.fillingFactorTol of cfg.fillingFactor, when
%     that target is set. This is what turns a three-variable box into a
%     roughly two-dimensional shell, and it is measured from the same
%     polygons COMSOL will be handed rather than from a closed form - the
%     naive 3*w*r overstates the hole area by ~6% because the three arms
%     overlap at the centre.
%
% (4) Minimum SOLID feature >= cfg.minSolidFeatureNm. The etched side is
%     covered by (1); this is the dielectric side - the thinnest wall left
%     standing, between this cell's hole and a neighbouring cell's. Measured
%     over the six nearest lattice translates, since the wall does not exist
%     until the lattice is applied.
%
%     Worth knowing: over the SHIPPED bounds this never binds. The tightest
%     corner (a = 600, r = 250, w = 120 nm) leaves a 107 nm wall against the
%     50 nm limit, and the current test_Boomerang geometry leaves 148 nm. It
%     is a safety net that becomes live if the bounds widen or r is pushed
%     toward the containment limit - not a constraint shaping today's search.
%     Set cfg.minSolidFeatureNm = 0 to skip it and its cost.
%
% (5) Hole-to-cell clearance >= cfg.holeEdgeMarginNm. Test (2) is already an
%     exact non-intersection test, so this only adds a margin; 0 makes it a
%     no-op. It constrains r and w JOINTLY - the farthest point of the hole
%     sits at hypot(w/2, r) from the hole centre, so fattening an arm spends
%     edge clearance exactly as lengthening one does, and neither variable
%     alone expresses the limit.
holeWidthNm = x.w;                      % narrowest etched feature
tf = holeWidthNm(:) >= cfg.minFeatureNm;

% (2) and (3) both need the reconstructed geometry, so they share one call per
%     row. 'FeatureSizes',false skips calcFillingFactor's O(n^2) wall
%     measurement, which is far too slow here - bayesopt evaluates this on
%     thousands of candidate rows per iteration. The areas, the filling factor
%     and armsOverhang all survive the fast path.
%
%     Rows already failing (1) are skipped: the geometry work is the expensive
%     part and their verdict cannot change.
useFill  = ~isempty(cfg.fillingFactor);
% Whether the expensive metrics are needed at all. Both thresholds default to
% something meaningful, but a caller can zero them out and skip the cost.
needSlow = (cfg.minSolidFeatureNm > 0) || (cfg.holeEdgeMarginNm > 0);
for ii = 1:height(x)
    if ~tf(ii), continue, end

    Pi = struct('a',double(x.a(ii))*1e-9, ...
                'w',double(x.w(ii))*1e-9, ...
                'r',double(x.r(ii))*1e-9);
    if isfield(cfg,'holeCentreFrac'), Pi.holeCentreFrac = cfg.holeCentreFrac; end

    try
        ffi = calcFillingFactor(Pi,'FeatureSizes',false);
    catch
        % A degenerate candidate (hole swallowing the cell) throws rather than
        % returning a ratio. That is infeasible, not a crash of the study.
        tf(ii) = false;
        continue
    end

    % (2) Arms must stay inside the cell. MEASURED, not the old analytic
    %     r < sqrt(3)*a/4 test - that assumed the hole sits on the cell
    %     centroid, which stopped being true when resolveHoleCentreFrac
    %     introduced the default 0.4 offset. armsOverhang compares the hole
    %     area before and after the clip, so it is exact for any hole position.
    if ffi.armsOverhang
        tf(ii) = false;
        continue
    end

    % (3) Filling factor inside the requested band.
    if useFill && abs(ffi.fillingFactor - cfg.fillingFactor) > cfg.fillingFactorTol
        tf(ii) = false;
        continue
    end

    % (4) and (5) need calcFillingFactor's EXPENSIVE metrics, so they run last
    %     and only on rows that survived everything above. That ordering is the
    %     whole reason this is affordable: the filling-factor band admits only a
    %     few percent of the box, so the O(n^2) wall measurement runs on a few
    %     percent of the candidates rather than all of them. Cheapest and most
    %     selective first.
    if ~needSlow, continue, end

    try
        ffs = calcFillingFactor(Pi,'FeatureSizes',true);
    catch
        tf(ii) = false;
        continue
    end

    % (4) Minimum SOLID feature: the thinnest dielectric wall in the tiled
    %     structure, running between this cell's hole and a neighbouring
    %     cell's. Test (1) covers the etched side - the narrowest slit, = w -
    %     but a process has two limits and this is the other one. It is not
    %     derivable from a, r, w: the wall only exists once the lattice is
    %     applied, and it collapses as r grows toward the containment limit,
    %     well before armsOverhang trips. Nothing checked it until now.
    if ffs.minSolidFeature*1e9 < cfg.minSolidFeatureNm
        tf(ii) = false;
        continue
    end

    % (5) Clearance between the hole and the cell edge. armsOverhang in (2) is
    %     already an exact "does not intersect" test, so this only adds a
    %     MARGIN: with cfg.holeEdgeMarginNm = 0 it is a no-op and (2) alone
    %     decides. The margin binds r and w together rather than either alone,
    %     because the farthest point of the hole sits at hypot(w/2, r) from the
    %     hole centre - widening the arms eats edge clearance exactly as
    %     lengthening them does. Measured from the boundary, so it stays correct
    %     for any holeCentreFrac.
    if ffs.holeToCellMargin*1e9 < cfg.holeEdgeMarginNm
        tf(ii) = false;
    end
end
tf = tf(:);
end

% -------------------------------------------------------------------------

function [objective,coupledConstraints,userData] = boomerangObjective(x,cfg)
%BOOMERANGOBJECTIVE Solve one boomerang unit cell and score its bandgap.
%
% x is a single-row table supplied by bayesopt. Returns the negated fitness,
% since bayesopt minimizes.

coupledConstraints = [];    % all constraints are deterministic (XConstraintFcn)
userData = struct();

% The backend is recorded on every point bayesopt keeps, exactly as the Python
% result record carries mech_backend. results.UserDataTrace{k}.backend then
% answers "was this number real?" for any evaluation in a saved study, without
% having to remember what cfg said at the time.
userData.backend = cfg.solverBackend;

P = boomerangParams(x,cfg);

fprintf('\n----- evaluating a=%dnm r=%dnm w=%dnm (th=%dnm fixed) -----\n', ...
    x.a, x.r, x.w, round(cfg.th*1e9));
if cfg.isDryRun
    fprintf('  [DRY RUN] backend = %s : SYNTHETIC result, not physics\n', ...
        cfg.solverBackend);
end

tEval = tic;
[gapData,wasCached] = solveBoomerangBands(P,cfg);
evalTime = toc(tEval);

if isempty(gapData)
    % A genuine solver failure, as opposed to a valid solve that simply found
    % no gap. NaN tells bayesopt to treat this point as an error rather than
    % as a real objective value of zero, so the GP is not taught that this
    % region is merely mediocre.
    objective = NaN;
    userData.status = 'solver-error';
    logIteration(cfg,x,struct('midGap',NaN,'gapSize',NaN,'gapRatio',NaN, ...
        'penalty',NaN,'fitness',NaN),objective,wasCached,evalTime);
    return
end

[fitness,detail] = gapFitness(gapData,cfg);
objective = -fitness;

userData.status   = detail.status;
userData.midGap   = detail.midGap;
userData.gapSize  = detail.gapSize;
userData.gapRatio = detail.gapRatio;

fprintf('  mid-gap   = %s GHz\n', num2str(detail.midGap/1e9,'%.3f'));
fprintf('  gap ratio = %s %%\n',  num2str(detail.gapRatio*100,'%.3f'));
fprintf('  penalty   = %s\n',     num2str(detail.penalty,'%.4f'));
fprintf('  objective = %.6f (%s, %.2f mins)\n', ...
    objective, detail.status, evalTime/60);

logIteration(cfg,x,detail,objective,wasCached,evalTime);
end

% -------------------------------------------------------------------------

function P = boomerangParams(x,cfg)
%BOOMERANGPARAMS Build the P struct for one design point.
%
% Ported from sweep_boomerang_code.m/sweep_boomerang_thickness so the physics
% settings are identical to the existing sweeps; a/r/w now come from the
% optimizer instead of being hard-coded, while th is fixed at cfg.th.

nm = 1e-9;

P = struct();

%% unit cell params
P.xsect    = 'rect';
P.beamMat  = 'diamond';
P.celltype = 'boomerang';
P.unitcell = 'hexagonal';
P.a  = double(x.a)*nm;      % lattice constant; side of the rhombic cell
P.w  = double(x.w)*nm;      % hole arm width (narrowest etched feature)
P.r  = double(x.r)*nm;      % hole arm length, cell centre to arm tip
P.th = cfg.th;              % full slab thickness in z - FIXED, see cfg.th
P.r1 = cfg.r1;              % fillet radius, inner corners (arms meet at centre)
P.r2 = cfg.r2;              % fillet radius, outer arm tips
P.nperiod    = 1;           % no. of periods to simulate for
P.holeatedge = 0;           % 1/0 for hole at edge/center of unit cell

P.kpts   = cfg.kpts;        % no. of k-points, EXCLUDING gamma point
P.nbands = cfg.nbands;      % no. of bands to solve for

P.solveasym       = 1;      % 1 to solve for antisymmetric bands
P.completeBandGaps = 1;     % 1 to plot complete bandgaps (across all symmetries)
P.plotgeom    = cfg.plotgeom;
P.savedat     = 1;          % must stay 1 - this is what makes caching work
P.savebndplot = cfg.savebndplot;
P.saveplots   = 0;          % 1 to save displacement and strain profiles
P.saveMPH     = 0;
P.bandStruct_2D = 1;        % 1 to simulate 2D band structures

%% mechanical simulation parameters
P.mbeveny       = 0;        % 1 to find even mechanical mode about y
P.TwoSymPlanes  = 0;
P.zSymCondition = 1;        % symmetry condition is wrt the z plane
P.mbevenz       = 1;        % 1 to find even mechanical mode about z
P.freq          = 0;        % 0 for bandstructure simulations
P.meshSize      = cfg.meshSize;
P.fixed_bc      = 0;        % 1 to fix the xz boundaries at y = +/- w/2

P.anisoMat = 1;
P.rxtal    = 45;            % ccw rotation of elasticity matrix [deg] from the
                            % <100> in-plane direction about the <100> normal

%% degree-of-freedom cap, to bound the solve time
P.max_dof = cfg.maxDof;

P.datLoc = cfg.datLoc;

% Set fileBase explicitly, reproducing byte-for-byte the string solveBands
% would build for celltype 'boomerang' (including its missing '_' after the r
% field). Two reasons: it lets this script predict the cache path so a
% short-circuited solve can be recovered, and it keeps the filenames
% interchangeable with data already produced by the sweep scripts.
P.fileBase = ['boomerang_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
    'r_',num2str(P.r*1e9,'%.0f'),'nm',...
    'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
    'th_',num2str(P.th*1e9,'%.0f'),'nm_', ...
    'r1_',num2str(P.r1*1e9,'%.0f'),'nm_', ...
    'r2_',num2str(P.r2*1e9,'%.0f'),'nm_'];

% Guard (2) of the cache-isolation guards. In a dry run the SAME geometry gets a
% different filename, so a synthetic _bds.mat and a real one for one design can
% coexist in one folder without either shadowing the other - which matters
% because guard (1) (the separate folder) survives only until someone tidies up
% by hand, and file moves are exactly the operation nobody thinks twice about.
%
% solveBands has this mechanism already: it prepends P.prefname to the fileBase
% it builds. It only does so when P.fileBase is UNSET, though, and the whole
% point of the block above is that this script sets fileBase itself. So
% P.prefname is set for anything that reads it, and the same transformation
% solveBands would have applied ([prefname '_' fBase]) is applied here by hand.
% Keep the two spellings identical if either ever changes.
if cfg.isDryRun
    P.prefname = cfg.dryRunPrefname;
    P.fileBase = [P.prefname,'_',P.fileBase];
end
end

% -------------------------------------------------------------------------

function [gapData,wasCached] = solveBoomerangBands(P,cfg)
%SOLVEBOOMERANGBANDS Produce the gap struct for one design, via cfg.solverBackend.
%
% This is the single seam between the optimization loop and whatever is standing
% in for the physics. Everything above it - the objective, gapFitness, the log,
% the checkpoint, the figures - is identical in every mode, which is the only
% reason a dry run is worth anything as a test: if the cheap path went round the
% loop instead of through it, passing the dry run would prove nothing.
%
% Contract, unchanged from before this switch existed: gapData is the ds.full
% struct (fields midGap, gapSize), or [] to mean "this evaluation failed";
% wasCached is true when the result came off disk rather than being computed.
%
% The otherwise-branch is intentionally an error and not a warning-plus-default.
% Defaulting to 'comsol' would make a typo cost hours of solve time nobody asked
% for; defaulting to 'surrogate' would fabricate a study. Both are worse than
% stopping. assertKnownBackend already ran at configuration time, so this is the
% second line of defence, for a cfg that was hand-edited or loaded from an old
% results .mat afterwards.

% A cfg saved before this option existed has no solverBackend field, and
% `resume` on such a checkpoint re-resolves the stored objective handle against
% the CURRENT file - so the old cfg meets the new code. Left unguarded that is a
% bare "Unrecognized field name 'solverBackend'" thrown from inside a resumed
% study, which says nothing about the cause. Note this is not the unknown-backend
% case and is not defaulted for the same reason: a missing field means the caller
% has not decided, and deciding on their behalf is what the whole switch exists to
% avoid.
if ~isfield(cfg,'solverBackend')
    error('bayesopt_boomerang:missingBackend', ...
        ['cfg has no solverBackend field. This cfg predates the dry-run ' ...
         'option - most likely a checkpoint or results .mat written by an ' ...
         'earlier version of this script. Re-run the script (which rebuilds ' ...
         'cfg and the objective handle) rather than resuming that state; ' ...
         'existing _bds.mat files under cfg.datLoc will be reused, so no ' ...
         'solve time is lost.']);
end

switch cfg.solverBackend
    case 'comsol'
        [gapData,wasCached] = solveBandsViaComsol(P,cfg);
    case 'surrogate'
        [gapData,wasCached] = solveBandsViaSurrogate(P,cfg);
    case 'stub'
        [gapData,wasCached] = solveBandsViaStub();
    otherwise
        error('bayesopt_boomerang:unknownBackend', ...
            ['Unknown cfg.solverBackend ''%s''. Expected ''comsol'', ' ...
             '''surrogate'' or ''stub''.'],char(string(cfg.solverBackend)));
end

% Optional pacing, cheap backends only. A real solve already takes minutes, so
% slowing it down further would be nothing but cruelty; a surrogate evaluation
% takes microseconds, which is too fast to watch anything happen.
if cfg.isDryRun && cfg.dryRunDelay > 0
    pause(cfg.dryRunDelay);
end
end

% -------------------------------------------------------------------------

function [gapData,wasCached] = solveBandsViaComsol(P,cfg)
%SOLVEBANDSVIACOMSOL Call solveBands and always return a usable gap struct.
%
% gapData is the ds.full struct (fields midGap, gapSize), or [] if the solve
% failed. wasCached is true when the result came from an existing _bds.mat.
%
% This wrapper exists because solveBands has two distinct return shapes: on a
% fresh solve it returns ds with a .full field, but when
% [datLoc fileBase '_bds.mat'] already exists it prints "Data folder exists in
% working directory" and returns a stub with only empty .sym/.asym fields.
% Reading ds.full off that stub throws "Unrecognized field name". Over an
% optimization run - where re-visiting a design point is routine - that would
% abort the whole study, so the cached .mat is reloaded instead.

gapData   = [];
wasCached = false;

matPath = fullfile(P.datLoc,[P.fileBase,'_bds.mat']);

% Answer the cache question HERE rather than leaving it to solveBands.
% solveBands decides it with strcmp(P.datLoc(end),'\'), so on a platform whose
% filesep is '/' it appends a literal backslash to the path it then uses for its
% own dir() probe and mkdir. On macOS and Linux that yields a stray directory
% named '\' and a probe that can never match the file solveBands itself wrote,
% so every design would be re-solved however many times it was visited - the
% cache would look present and do nothing. Owning the decision here makes the
% cache behave identically on all three platforms, and on Windows it changes
% nothing except skipping one call into solveBands that would have returned its
% stub anyway.
if isfile(matPath)
    wasCached = true;
    gapData   = loadBandDataFromCache(matPath,cfg);
    if ~isempty(gapData)
        fprintf('  (reused cached band structure)\n');
    end
    return
end

figsBefore = findall(groot,'Type','figure');

try
    ds = solveBands(P);
catch ME
    warning('bayesopt_boomerang:solveFailed', ...
        'solveBands failed for %s: %s', P.fileBase, ME.message);
    closeNewFigures(figsBefore,cfg);
    return
end

closeNewFigures(figsBefore,cfg);

if isfield(ds,'full') && isstruct(ds.full)
    gapData = ds.full;
    return
end

% solveBands returned its cache-hit stub even though nothing was on disk when we
% looked a moment ago - a concurrent run, or a filename we failed to predict.
% Read whatever is there now rather than throwing the evaluation away.
wasCached = true;
if ~isfile(matPath)
    warning('bayesopt_boomerang:noGapData', ...
        'solveBands returned no .full field and no cache at %s', matPath);
    return
end
gapData = loadBandDataFromCache(matPath,cfg);
end

% -------------------------------------------------------------------------

function gapData = loadBandDataFromCache(matPath,cfg)
%LOADBANDDATAFROMCACHE Read ds.full out of a cached _bds.mat, provenance-checked.

gapData = [];
cached  = load(matPath,'ds');

% Guard (3): every cache load is checked for the synthetic provenance marker
% before its numbers are believed. Guards (1) and (2) already make this
% unreachable in normal use; it is here because they are both conventions about
% filenames, and a convention is only as strong as the last person who moved a
% file. This check is not - it inspects the data itself.
if isfield(cached,'ds')
    assertBandDataProvenance(cached.ds,matPath,cfg);
end

if isfield(cached,'ds') && isfield(cached.ds,'full')
    gapData = cached.ds.full;
else
    warning('bayesopt_boomerang:badCache', ...
        'Cached file %s contains no ds.full', matPath);
end
end

%% ========================================================================
%  Dry-run backends - local functions
%  ========================================================================
% Everything below this banner is inert unless cfg.solverBackend is set to one
% of the cheap tiers. None of it is reachable from a 'comsol' run: the switch in
% solveBoomerangBands is the only caller, apart from the provenance check, which
% is reachable from a real run precisely so that it can refuse synthetic data.

function [gapData,wasCached] = solveBandsViaSurrogate(P,cfg)
%SOLVEBANDSVIASURROGATE Cheap tier 1: analytic fake bands, full success path.
%
% Deliberately shaped like solveBandsViaComsol, including the cache-hit branch,
% because the point is to walk the same plumbing rather than to shortcut it. The
% synthetic .mat it writes lives under the dry-run folder and carries the DRYRUN
% filename prefix and the in-file provenance marker - see the CACHE ISOLATION
% note in the header block.

gapData   = [];
wasCached = false;

% Optional simulated failure, checked FIRST so that a "failed" design never
% reaches the disk. A failure that still left a cache file behind would heal
% itself on the next visit to the same point, which is the opposite of
% deterministic.
if simulatedSolverFailure(P,cfg)
    warning('bayesopt_boomerang:simulatedFailure', ...
        ['SIMULATED solver failure (cfg.dryRunFailEvery = %d) for %s. ' ...
         'This is a dry-run test of the NaN path, not a real problem.'], ...
        cfg.dryRunFailEvery,P.fileBase);
    return
end

matPath = fullfile(P.datLoc,[P.fileBase,'_bds.mat']);

% Cache hit. Worth reproducing rather than skipping: the `cached` column of the
% log, the (:) shapes that come back off disk rather than out of memory, and
% loadCachedBandData's assumptions about what a saved ds looks like are all
% things a dry run should be able to check.
if cfg.dryRunSaveBands && isfile(matPath)
    cached = struct();
    try
        cached = load(matPath,'ds');
    catch ME
        warning('bayesopt_boomerang:badDryRunCache', ...
            'Could not load synthetic cache %s: %s',matPath,ME.message);
    end
    if isfield(cached,'ds')
        assertBandDataProvenance(cached.ds,matPath,cfg);
        if isfield(cached.ds,'full') && isstruct(cached.ds.full)
            gapData   = cached.ds.full;
            wasCached = true;
            fprintf('  (reused cached SYNTHETIC band structure)\n');
            return
        end
    end
end

ds = surrogateBoomerangBands(P,cfg);
gapData = ds.full;

if cfg.dryRunSaveBands
    try
        save(matPath,'ds');
    catch ME
        % Not fatal: the surrogate result is already in hand, and the only thing
        % lost is the cache-hit exercise on a later visit to this design.
        warning('bayesopt_boomerang:dryRunSave', ...
            'Could not write synthetic band data to %s: %s',matPath,ME.message);
    end
end
end

% -------------------------------------------------------------------------

function [gapData,wasCached] = solveBandsViaStub()
%SOLVEBANDSVIASTUB Cheap tier 2: no bands at all, every evaluation "fails".
%
% The counterpart of mech_backend='surrogate_stub' in
% python-scripts/src/objective.py, which returns NaN and skips the stage rather
% than estimating it. Two cheap tiers are not redundant: 'surrogate' walks the
% SUCCESS path (gapFitness scores a real-shaped gap, figures get drawn from data)
% while this one walks the FAILURE path - objective = NaN, status
% 'solver-error', the red error markers on the convergence panel, and the
% band-data-unavailable degradation inside visualizeBestDesign. A dry run that
% only ever succeeded would leave all of that untested, and those are the paths
% that only ever execute on a bad day, i.e. the ones most likely to be broken.
%
% Takes no arguments on purpose: there is nothing about the design or the
% configuration that could change the answer.

gapData   = [];
wasCached = false;
fprintf(['  [stub backend] band solve SKIPPED - reporting a solver failure ' ...
         '(objective = NaN)\n']);
end

% -------------------------------------------------------------------------

function ds = surrogateBoomerangBands(P,cfg)
%SURROGATEBOOMERANGBANDS Analytic stand-in for a band structure. NOT PHYSICS.
%
% =========================================================================
%  THIS FUNCTION DOES NOT COMPUTE A BAND STRUCTURE. It evaluates a handful of
%  closed-form expressions chosen to LOOK like one. The frequencies it returns
%  are not approximations of the boomerang cell's modes, they are not converged,
%  they are not even derived from an elastic model - there is no eigenproblem
%  anywhere below. Nothing that comes out of here may be used to judge a design.
%  Its only purpose is to make the optimization loop run in microseconds so the
%  loop itself can be debugged. Same disclaimer, same reasons, as the docstring
%  of python-scripts/src/optical_surrogate.py ("a cheap, dependency-free
%  estimator ... use it only as a cheap pre-screen and to exercise the loop") -
%  except that this one is weaker still, being a fit to nothing at all.
% =========================================================================
%
% WHAT IT MUST GET RIGHT, and why:
%
% 1. SHAPE. Returns the same ds a real 2-sector solve does: ds.sym and ds.asym
%    each with F [nk x nbands] in Hz and k_norm [nk x 1] running 0->3, and
%    ds.full with F, midGap and gapSize. nk = 3*kpts+1, matching the
%    Gamma-M-K-Gamma circuit runBands_2D walks plus the wrapped final Gamma. Get
%    this wrong and the dry run tests the surrogate instead of the pipeline.
%
% 2. DETERMINISM. No rand, no clock, no counters - the value depends on
%    (a,r,w,th) and on cfg.kpts/cfg.nbands and on nothing else. bayesopt is told
%    'IsObjectiveDeterministic',true, and a surrogate that answered differently
%    on a repeated design point would violate that contract and quietly corrupt
%    the GP fit, which is exactly the class of bug a dry run exists to find
%    rather than to introduce.
%
% 3. A NON-TRIVIAL OPTIMUM. A flat or monotone surface would make the
%    convergence trace and the four design slices meaningless as tests - a plot
%    of nothing renders perfectly well. So the fitness has an interior maximum in
%    all four variables, produced by five competing gap-width factors:
%      * air filling fraction  fill  -> best near 0.0823
%      * slab aspect ratio     th/a  -> best near 0.390
%      * arm aspect ratio      w/r   -> best near 0.425
%      * hole reach            r/a   -> best near 0.250
%      * absolute thickness    th    -> best near 312 nm  (see the note at the
%                                       gapStrength expression for why one
%                                       non-dimensionless factor is needed)
%    and one scaling law:
%      * frequencies ~ 1/a           -> the Gaussian target-frequency penalty in
%                                       gapFitness has real work to do, and it
%                                       FIGHTS the gap-width factors rather than
%                                       agreeing with them.
%    Every one of those optima is placed so the peak lands at the CENTRE of the
%    default bounds: a = 800, r = 200, w = 85, th ~ 312 nm, fitness ~ 0.24,
%    mid-gap 13.0 GHz. Measured on a grid scan of the feasible box, the peak beats
%    the best design available on each bound FACE by 2.5x (a), 1.9x (r), 2.4x (w)
%    and 1.5x (th) - so all four design slices show a real interior peak rather
%    than a plateau running into a wall. About 64% of the feasible box has NO
%    complete gap at all, so the gapless branch of gapFitness is exercised too.
%    If you change cfg.bounds, cfg.targetFreq or cfg.sigma, re-check that the
%    optimum is still interior: a surrogate whose optimum has drifted onto a
%    bound tests considerably less than it appears to.
%
% 4. SELF-CONSISTENCY. midGap and gapSize are NOT written by hand: the
%    synthesized band matrix is passed through the real findGaps, the same way
%    solveBands does it. Hand-written gap values could disagree with the bands
%    plotted beside them, and would leave findGaps itself untested.
%
% The trends are chosen to be plausible rather than correct: bigger cells give
% lower frequencies, more air softens the medium, and there is a best filling
% fraction beyond which the gap closes again. Real boomerang cells do behave
% qualitatively like that. The numbers still mean nothing.

nm   = 1e9;
aNm  = P.a*nm;
rNm  = P.r*nm;
wNm  = P.w*nm;
thNm = P.th*nm;

% --- three dimensionless shape descriptors -------------------------------
% Rhombic primitive cell of side a, 60 deg interior angle.
cellArea = (sqrt(3)/2)*aNm^2;

% Union of the three w-by-r arms meeting at 120 deg. The -0.75*w^2 is a crude
% correction for their mutual overlap at the cell centre; the exact constant is
% irrelevant, all that is needed is a quantity that rises with w and r and falls
% with a. (drawBoomerangCell builds the same shape exactly, with polyshape - not
% reused here because the objective path should not depend on graphics.)
holeArea   = max(3*wNm*rNm - 0.75*wNm^2, 0);
fill       = min(max(holeArea/cellArea, 0), 0.85);   % air filling fraction
slabAspect = thNm/aNm;      % slab thickness against cell size
armAspect  = wNm/rNm;       % how stubby each arm is
holeReach  = rNm/aNm;       % how far the arms get across the cell

% --- frequency scale: everything scales as 1/a ---------------------------
% cRef is NOT a material property. It is a fitting constant, picked so that the
% engineered gap lands on the default cfg.targetFreq of 13 GHz at a = 800 nm,
% which is the middle of the default a range - so the frequency penalty pulls the
% optimum to the centre of the box rather than to a bound. (Measured: 12.995 GHz
% at the intended optimum.) It happens to be of the same order as a shear speed
% in diamond; that is a coincidence, not a justification.
cRef = 9650;                            % [m/s]
cEff = cRef*sqrt(1-fill);                % more air -> softer -> slower
f1   = cEff/(2*aNm*1e-9);                % [Hz] first-rung frequency scale, ~1/a

% --- gap strength, in [0,1] ---------------------------------------------
% Product of four Gaussians in four DIFFERENT dimensionless ratios, so the
% optimum is interior in every variable and the surface is smooth - a GP has
% something to fit and a finite-difference probe of the loop would behave.
% gapStrength = 1 is the widest gap the surrogate can produce; below about 0.34
% the gap closes entirely (see the ladder below), which is what gives the dry
% run a genuine no-complete-gap region to fall into - about 64% of the feasible
% box, on a grid scan.
%
% Every optimum below is set so that the peak lands at the CENTRE of each of the
% default bounds (a = 800, r = 200, w = 85, th ~ 312 nm). That is not cosmetic: an
% optimum near a bound makes the "best design on that bound face" almost as good
% as the peak, and the corresponding design slice then looks flat and reports
% "[at bound]". Centring the peak is what buys the 1.5-2.5x margins quoted above.
%
% Four RATIOS rather than two or three because the smaller sets left a RIDGE:
% fill and armAspect between them can be held at their optima while (a,r,w) all
% slide together, so the optimum was a curve and the w design slice came out
% nearly flat. holeReach pins r against a and closes it. Four constraints on four
% variables also means no design satisfies all of them exactly, so the peak is a
% genuine compromise rather than a point where everything happens to agree.
%
% The FIFTH factor is different in kind and frankly arbitrary: a preference for
% an absolute thickness near 320 nm. It is here for one reason. Every
% dimensionless ratio is blind to scaling a and th together, so with ratios alone
% a design pushed to the th bound simply rescales a to keep th/a optimal and pays
% nothing but the (broad, sigma = 5 GHz) frequency penalty - measured margin over
% the th bound face was 1.05x, i.e. the th slice was very nearly featureless. One
% absolute length scale breaks that degeneracy and takes the margin to ~2.6x.
% Note the consequence: unlike the four ratios, this factor does not follow
% cfg.th, so moving that fixed thickness far from ~320 nm will push the optimum
% onto a bound. Re-run a grid scan if you move them.
gapStrength = exp(-((fill       - 0.0823)/0.065)^2) ...
            * exp(-((slabAspect - 0.390 )/0.110)^2) ...
            * exp(-((armAspect  - 0.425 )/0.200)^2) ...
            * exp(-((holeReach  - 0.250 )/0.065)^2) ...
            * exp(-((thNm       - 312   )/75   )^2);
gapStrength = min(max(gapStrength,0),1);

% --- ladder of pass bands ------------------------------------------------
% The spectrum is built as a ladder of 2*nbands non-dispersive "pass bands",
% each given a k dependence afterwards. Rung centres go as m^0.55 so successive
% bands crowd together at high order, the way a real folded spectrum does; a
% linear ladder would put the top band at 30x the first and squash the
% interesting region of the band plot into the bottom 3% of the axes.
nBandTot = 2*cfg.nbands;
m        = 1:nBandTot;
mGap     = 3;                   % the gap is opened between rungs 3 and 4

cLad = m.^0.55;                 % rung centres, in units of f1
hLad = (0.08 + 0.07*(1-gapStrength)).*cLad;     % rung half-widths

% Rungs below the engineered gap are deliberately made wide enough to overlap
% each other. Without this the ladder's own spacing leaves an accidental gap
% between rungs 1 and 2 at every design, gapFitness would always find SOMETHING,
% and the 'no-complete-gap' status would be unreachable in a dry run. Widening
% the low rungs is also the honest choice physically: the acoustic branches near
% Gamma are the ones that really do overlap.
hLad(m < mGap) = 3*hLad(m < mGap);

% Open the gap by lifting every rung above mGap. Half-widths were computed from
% the UNSHIFTED centres on purpose, so the shift moves the gap edges and not the
% band widths.
cLad(m > mGap) = cLad(m > mGap) + 0.55*gapStrength;

% --- k dependence --------------------------------------------------------
% k_norm runs 0->3 over Gamma-M-K-Gamma, so the dispersion must return to its
% starting value at k=3. A single cosine would do that but looks obviously fake;
% adding a second harmonic bends the curves without breaking periodicity. It is
% then renormalised to span exactly [-1,1] so each rung's extremes are exactly
% its centre +/- its half-width - which is what makes the resulting gaps
% predictable, and therefore testable, in a harness.
nk    = 3*cfg.kpts + 1;
kNorm = linspace(0,3,nk)';
theta = 2*pi*kNorm/3;
raw   = 0.8*cos(theta) + 0.2*cos(2*theta);      % range [-0.6, 1]
shape = (raw + 0.6)/0.8 - 1;                    % range [-1, 1] exactly

% Alternating curvature so neighbouring bands bend opposite ways, as folded
% bands do. Outer product: [nk x 1] * [1 x nBandTot] -> [nk x nBandTot].
curvature = (-1).^m;
F = f1*(repmat(cLad,nk,1) + shape*(hLad.*curvature));

% --- split into symmetry sectors ----------------------------------------
% Odd rungs to the even-about-z sector, even rungs to the odd-about-z sector, so
% each gets exactly cfg.nbands columns and the two edges of the engineered gap
% come from DIFFERENT sectors. That last detail matters: it makes the gap a
% genuinely complete gap - empty in both sectors - rather than one that only
% looks complete because a single sector was plotted.
symSector        = struct();
symSector.k_norm = kNorm;
symSector.F      = F(:,1:2:end);

asymSector        = struct();
asymSector.k_norm = kNorm;
asymSector.F      = F(:,2:2:end);

% Per-sector gaps, then the complete gaps, exactly as solveBands does it -
% including the [sym asym] column order of full.F.
[symSector.midGap, symSector.gapSize]  = findGaps(symSector);
[asymSector.midGap,asymSector.gapSize] = findGaps(asymSector);

fullBands    = struct();
fullBands.F  = [symSector.F, asymSector.F];
[fullBands.midGap,fullBands.gapSize] = findGaps(fullBands);

ds      = struct();
ds.sym  = symSector;
ds.asym = asymSector;
ds.full = fullBands;

% Provenance marker - guard (3). Travels inside the data, so it survives being
% renamed, moved, emailed or loaded by a script that has never heard of dry runs.
% assertBandDataProvenance refuses to let a real run consume anything carrying
% it. Kept as plain fields on ds rather than a nested struct so that a bare
% `load(f); ds` at the command line shows the warning immediately.
ds.isSynthetic      = true;
ds.syntheticBackend = cfg.solverBackend;
ds.syntheticWarning = ['SYNTHETIC band data from surrogateBoomerangBands in ' ...
    'bayesopt_boomerang.m. This was NEVER solved in COMSOL. The frequencies ' ...
    'are analytic placeholders for loop debugging and mean nothing physically. ' ...
    'Do not use for design, publication, or comparison with real sweeps.'];
ds.syntheticCreated = char(datetime('now'));
end

% -------------------------------------------------------------------------

function tf = simulatedSolverFailure(P,cfg)
%SIMULATEDSOLVERFAILURE Deterministic pseudo-failure for exercising the NaN path.
%
% Keyed on the DESIGN, not on an evaluation counter or a random draw. bayesopt is
% told 'IsObjectiveDeterministic',true, so the same design must always produce
% the same answer; a counter would make the second visit to a design succeed
% where the first failed, teaching the GP that the objective is noisy and
% inventing a bug the real solver does not have.
%
% The sum a+r+w+th is a poor hash but a perfectly good one here: designs are
% integers in nm, the bayesopt trace visits them in no particular order, so
% "every Nth value of the sum" scatters failures through the run without
% clustering them anywhere the optimizer would notice.

tf = false;
if ~cfg.isDryRun || cfg.dryRunFailEvery <= 0
    return
end
nm  = 1e9;
key = round(P.a*nm) + round(P.r*nm) + round(P.w*nm) + round(P.th*nm);
tf  = mod(key,round(cfg.dryRunFailEvery)) == 0;
end

% -------------------------------------------------------------------------

function assertKnownBackend(backend)
%ASSERTKNOWNBACKEND Reject an unrecognised cfg.solverBackend, loudly.
%
% The MATLAB equivalent of `raise ValueError(f"unknown mech_backend ...")` in
% python-scripts/src/objective.py, and it is an error rather than a warning for
% the same reason: there is no safe default. Falling back to 'comsol' would spend
% hours of solve time on a typo; falling back to 'surrogate' would hand back a
% folder full of fabricated results that look exactly like a study. Stopping is
% the only option that cannot mislead.

known = {'comsol','surrogate','stub'};
if ~(ischar(backend) || (isstring(backend) && isscalar(backend))) ...
        || ~any(strcmp(char(backend),known))
    error('bayesopt_boomerang:unknownBackend', ...
        ['Unknown cfg.solverBackend. Expected one of: %s.\n' ...
         'Got: %s\n' ...
         '  ''comsol''    - the real COMSOL solve (the only real physics)\n' ...
         '  ''surrogate'' - analytic fake bands, for debugging the loop\n' ...
         '  ''stub''      - no bands at all, every evaluation returns NaN'], ...
        strjoin(cellfun(@(s) ['''',s,''''],known,'UniformOutput',false),', '), ...
        formatBackendForError(backend));
end
end

% -------------------------------------------------------------------------

function txt = formatBackendForError(backend)
%FORMATBACKENDFORERROR Describe whatever was found in cfg.solverBackend.
%
% Split out so assertKnownBackend's message can name the offending value without
% assuming it is printable. The whole point of that function is to catch a value
% nobody expected, and a sprintf('%s',...) on a struct or a cell would throw from
% inside the error handler - replacing a clear diagnostic with a confusing one.

try
    if ischar(backend) || (isstring(backend) && isscalar(backend))
        txt = ['''',char(backend),''''];
    else
        txt = ['a ',class(backend),' of size ',mat2str(size(backend))];
    end
catch
    txt = '<unprintable>';
end
end

% -------------------------------------------------------------------------

function tf = isSyntheticBandData(ds)
%ISSYNTHETICBANDDATA True if a ds struct carries the dry-run provenance marker.
%
% Written to answer "false" for anything that is not clearly marked, including
% [], a non-struct, and every real _bds.mat ever written by solveBands - none of
% which have the field. The asymmetry is the safe one: an unmarked file is
% treated as real (and a real run proceeds), a marked file is treated as
% synthetic (and a real run stops).

tf = isstruct(ds) && isscalar(ds) && isfield(ds,'isSynthetic') ...
    && ~isempty(ds.isSynthetic) && all(logical(ds.isSynthetic(:)));
end

% -------------------------------------------------------------------------

function assertBandDataProvenance(ds,matPath,cfg)
%ASSERTBANDDATAPROVENANCE Refuse to mix synthetic and real band data.
%
% Guard (3) of the four cache-isolation guards, and the only one that inspects
% the data rather than trusting a filename convention. Called on every path that
% reads a _bds.mat back off disk.
%
% Synthetic data reaching a real run is an ERROR, not a warning. The failure mode
% being prevented is silent and permanent: fabricated frequencies entering a real
% study's GP, its log and its results .mat, indistinguishable afterwards from
% solved ones. Stopping a study that can be resumed from its checkpoint is a far
% smaller loss than finishing one whose numbers cannot be trusted.
%
% The reverse - real data reaching a dry run - is only a warning. It wastes the
% dry run (it is now testing whatever that file contains) but corrupts nothing,
% and it usually means cfg.datLoc was pointed at a real folder by hand, which is
% a thing someone might do on purpose to test the visualization against real
% bands.

if isSyntheticBandData(ds)
    if ~cfg.isDryRun
        error('bayesopt_boomerang:syntheticDataInRealRun', ...
            ['REFUSING to use %s: it carries the synthetic dry-run marker ' ...
             '(ds.isSynthetic), so these frequencies were never solved.\n' ...
             'A real study must not consume dry-run output. Delete or move ' ...
             'the file, or point cfg.datLoc at a clean folder, and resume ' ...
             'from cfg.statePath.'],matPath);
    end
elseif cfg.isDryRun
    warning('bayesopt_boomerang:realDataInDryRun', ...
        ['%s carries no synthetic marker, so it looks like REAL solved data ' ...
         'being read by a dry run (cfg.solverBackend = %s). The dry run is ' ...
         'now exercising that file rather than the surrogate.'], ...
        matPath,cfg.solverBackend);
end
end

% -------------------------------------------------------------------------

function printDryRunBanner(cfg)
%PRINTDRYRUNBANNER Unmissable console notice that nothing here is physics.
%
% Deliberately loud and deliberately repeated (start and end of the run). The
% console is where the numbers are first read, usually while they are still
% scrolling, and the single most expensive mistake this whole mode can cause is
% someone believing one of them.

bar = repmat('=',1,74);
fprintf('\n%s\n',bar);
fprintf('  *** DRY RUN - cfg.solverBackend = ''%s'' ***\n',cfg.solverBackend);
fprintf('  NO COMSOL SOLVE IS PERFORMED. Every frequency, gap, fitness and\n');
fprintf('  figure produced by this run is SYNTHETIC and means NOTHING\n');
fprintf('  physically. The mode exists to debug the optimization loop.\n');
fprintf('  Output is isolated under: %s\n',cfg.datLoc);
fprintf('  Set cfg.solverBackend = ''comsol'' for a real study.\n');
fprintf('%s\n',bar);
end

% -------------------------------------------------------------------------

function writeDryRunMarker(cfg)
%WRITEDRYRUNMARKER Leave a plain-text warning inside the dry-run output folder.
%
% For the reader who finds this folder in a month with no memory of the run and
% no console scrollback. The folder name and the filenames already say DRYRUN,
% but names get skimmed; a file that has to be actively ignored is better. Kept
% to plain text on purpose - readable from Finder, Explorer, git and a terminal
% without MATLAB.
%
% Best-effort: a failure here must not cost the run, since the run is a debugging
% exercise and the data it produces is worthless by design anyway.

markerPath = fullfile(cfg.datLoc,'DRYRUN_README.txt');
fid = fopen(markerPath,'wt+');
if fid < 0
    warning('bayesopt_boomerang:dryRunMarker', ...
        'Could not write the dry-run marker %s',markerPath);
    return
end
fprintf(fid,'SYNTHETIC DATA - DO NOT TRUST ANYTHING IN THIS FOLDER\r\n');
fprintf(fid,'====================================================\r\n\r\n');
fprintf(fid,'Written by bayesopt_boomerang.m running with\r\n');
fprintf(fid,'    cfg.solverBackend = ''%s''\r\n',cfg.solverBackend);
fprintf(fid,'on %s.\r\n\r\n',char(datetime('now')));
fprintf(fid,'No COMSOL solve was performed. Every band structure, gap,\r\n');
fprintf(fid,'fitness value, log row and figure in this folder was produced\r\n');
fprintf(fid,'by an analytic placeholder (surrogateBoomerangBands), or is a\r\n');
fprintf(fid,'deliberate failure (the ''stub'' backend). The numbers are not\r\n');
fprintf(fid,'approximations of the real cell - they are not physics at all.\r\n\r\n');
fprintf(fid,'This folder exists only to debug the optimization loop:\r\n');
fprintf(fid,'objective plumbing, constraints, logging, checkpointing and the\r\n');
fprintf(fid,'post-run figures, without spending hours of COMSOL time.\r\n\r\n');
fprintf(fid,'Every _bds.mat here additionally carries ds.isSynthetic = true,\r\n');
fprintf(fid,'and a real run refuses to load a file marked that way. Do not\r\n');
fprintf(fid,'strip that field, and do not move these files into a real\r\n');
fprintf(fid,'results folder.\r\n');
fclose(fid);
end

% -------------------------------------------------------------------------

function bannerLines = dryRunBannerText(cfg)
%DRYRUNBANNERTEXT One-line figure banner, or {} outside a dry run.
%
% Returns a cell array so callers can concatenate it onto an existing headline
% with no conditional of their own, and so that returning {} in the normal case
% leaves a real run's figures byte-for-byte as they were.

if ~cfg.isDryRun
    bannerLines = {};
    return
end
bannerLines = {sprintf(['*** DRY RUN: SYNTHETIC %s DATA - NOT PHYSICS - ' ...
    'DO NOT USE FOR DESIGN ***'],upper(cfg.solverBackend))};
end

% -------------------------------------------------------------------------

function [headlineText,fracHeight] = withDryRunHeadline(headlineText,fracHeight,cfg)
%WITHDRYRUNHEADLINE Prefix a figure headline with the synthetic-data banner.
%
% Also grows the strip addFigureHeadline reserves at the top of the figure, since
% the headline just gained a line and the whole reason that strip is sized
% explicitly is that a headline which outgrows it collides with the axes titles
% underneath.
%
% Returns its inputs untouched outside a dry run, including leaving a char
% headline as a char rather than promoting it to a 1x1 cellstr - so a real run's
% figures are unchanged rather than merely equivalent.

if ~cfg.isDryRun
    return
end
if ~iscell(headlineText)
    headlineText = {headlineText};
end
headlineText = [dryRunBannerText(cfg), headlineText];
fracHeight   = fracHeight + 0.035;
end

% -------------------------------------------------------------------------

function rgb = dryRunHeadlineColor(cfg)
%DRYRUNHEADLINECOLOR Red headline text in a dry run, and NO opinion otherwise.
%
% Colour is doing real work here rather than decoration: a PNG pasted into an
% email or a lab notebook arrives with no folder name, no filename and no console
% output attached, so the figure itself has to carry the disclaimer. Red plus the
% banner line survives that trip; a subtly different title does not.
%
% Returns [] - not black - outside a dry run, so that addFigureHeadline leaves the
% Color property untouched and a real run's figures are unchanged. See the note
% there for why "black" and "the default" are not the same colour.

if cfg.isDryRun
    rgb = [0.75 0 0];
else
    rgb = [];
end
end

%% ========================================================================
%  Shared helpers - local functions
%  ========================================================================
% Backend-independent: these run identically whichever backend produced the
% numbers, which is the property that makes the dry run a test of them.

function closeNewFigures(figsBefore,cfg)
%CLOSENEWFIGURES Close only the figures a solve opened.
%
% A plain "close all" would also kill the bayesopt live plots, so the figure
% list is differenced against a snapshot taken before the solve.

if ~cfg.closeSolveFigures
    return
end
figsNew = setdiff(findall(groot,'Type','figure'),figsBefore);
if ~isempty(figsNew)
    close(figsNew);
end
end

% -------------------------------------------------------------------------

function [fitness,detail] = gapFitness(gapData,cfg)
%GAPFITNESS Score a set of complete bandgaps against the target frequency.
%
%   fitness = max over complete gaps of
%             (gapSize/midGap) * exp(-((f_target - midGap)/sigma)^2)
%
% Every complete gap is scored and the best is kept. optimization_client.m
% instead tried to *select* one gap first, via
%     gap_ind = find(abs(target_freq-mFreqs) < mGaps);
%     if isempty(gap_ind)
%         gap_ind = find(min(abs(target_freq-mFreqs) - mGaps/2));
%     end
% which had two defects. The first line can match several gaps at once,
% leaving mGap/mFreq as vectors and making F a vector - fminsearch then errors
% or silently misbehaves. The second line is a no-op masquerading as an argmin:
% find(min(v)) evaluates the scalar min and returns 1 if it is nonzero and
% empty if it is exactly zero, so it always selects gap 1 (or nothing)
% regardless of which gap is closest. Scoring all gaps and taking the max
% removes the need for that selection step entirely.

detail = struct('midGap',NaN,'gapSize',NaN,'gapRatio',0, ...
    'penalty',0,'fitness',0,'status','no-complete-gap');

if ~isfield(gapData,'gapSize') || isempty(gapData.gapSize)
    fitness = 0;
    fprintf('  no complete gaps\n');
    return
end

midGaps  = gapData.midGap(:);
gapSizes = gapData.gapSize(:);

% Guard against a degenerate mid-gap frequency before dividing.
valid = midGaps > 0 & gapSizes > 0;
if ~any(valid)
    fitness = 0;
    return
end
midGaps  = midGaps(valid);
gapSizes = gapSizes(valid);

gapRatios = gapSizes ./ midGaps;
penalties = exp(-((cfg.targetFreq - midGaps)/cfg.sigma).^2);
scores    = gapRatios .* penalties;

[fitness,iBest] = max(scores);

detail.midGap   = midGaps(iBest);
detail.gapSize  = gapSizes(iBest);
detail.gapRatio = gapRatios(iBest);
detail.penalty  = penalties(iBest);
detail.fitness  = fitness;

% Flag whether the target frequency actually falls inside the winning gap,
% which is the distinction optimization_client.m printed as a message.
if abs(cfg.targetFreq - detail.midGap) < detail.gapSize/2
    detail.status = 'target-in-gap';
else
    detail.status = 'nearest-gap';
end
end

% -------------------------------------------------------------------------

function initIterationLog(logPath)
%INITITERATIONLOG Create the tab-separated evaluation log with a header row.
%
% The trailing 'backend' column is the log's answer to "is this row physics?".
% A dry run writes its log to a different folder entirely, so the column is not
% what keeps the two apart - but a log row that has been copied, pasted or
% readtable'd out of its folder has lost every other clue, and a table of gap
% ratios with no provenance is precisely the artifact that gets believed.

if isfile(logPath)
    return
end
fid = fopen(logPath,'wt+');
if fid < 0
    warning('bayesopt_boomerang:logOpen','Could not create log %s',logPath);
    return
end
fprintf(fid,['a_nm\tr_nm\tw_nm\tth_nm\tmidGap_GHz\tgapSize_GHz\t' ...
    'gapRatio\tpenalty\tfitness\tobjective\tcached\tevalTime_min\tstatus\t' ...
    'backend\r\n']);
fclose(fid);
end

% -------------------------------------------------------------------------

function logIteration(cfg,x,detail,objective,wasCached,evalTime)
%LOGITERATION Append one evaluation to the tab-separated log.
%
% Written and closed per evaluation, rather than buffered, so a run killed
% mid-study still leaves a complete record of everything solved so far.

fid = fopen(cfg.logPath,'at+');
if fid < 0
    warning('bayesopt_boomerang:logAppend','Could not append to %s',cfg.logPath);
    return
end
status = 'solver-error';
if isfield(detail,'status')
    status = detail.status;
end
% Backend last, matching the header written by initIterationLog. Appended rather
% than inserted so an existing log from before this column existed still lines up
% for its first 13 fields.
fprintf(fid,'%d\t%d\t%d\t%d\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%d\t%.3f\t%s\t%s\r\n', ...
    x.a, x.r, x.w, round(cfg.th*1e9), ...
    detail.midGap/1e9, detail.gapSize/1e9, detail.gapRatio, ...
    detail.penalty, detail.fitness, objective, wasCached, evalTime/60, status, ...
    logBackendLabel(cfg));
fclose(fid);
end

% -------------------------------------------------------------------------

function label = logBackendLabel(cfg)
%LOGBACKENDLABEL Value for the log's backend column.
%
% A dry-run row says e.g. 'surrogate-SYNTHETIC' rather than plain 'surrogate'.
% The suffix is redundant to anyone who knows what the backends are, and that is
% the point: the reader who does not know is the one at risk, and 'SYNTHETIC'
% needs no glossary. A real row is the bare 'comsol', so nothing is added to the
% format a real study writes beyond the backend name itself.

label = cfg.solverBackend;
if cfg.isDryRun
    label = [label,'-SYNTHETIC'];
end
end

% -------------------------------------------------------------------------

function stop = checkpointState(results,state,statePath)
%CHECKPOINTSTATE Save the BayesianOptimization object after every iteration.
%
% With multi-minute COMSOL evaluations, losing a partially finished study to a
% crashed LiveLink connection is expensive. The checkpoint can be reloaded and
% extended:
%   S = load(statePath); results = resume(S.results,'MaxObjectiveEvaluations',10);

stop = false;
switch state
    case {'iteration','done'}
        try
            save(statePath,'results');
        catch ME
            warning('bayesopt_boomerang:checkpoint', ...
                'Could not checkpoint to %s: %s',statePath,ME.message);
        end
end
end

%% ========================================================================
%  Visualization - local functions
%  ========================================================================
% Structure mirrors _run_characterization in run_opt_comsol.py: one composite
% figure assembled from single-purpose panel painters, then the same painters
% reused at full size in standalone figures. Keeping each panel in its own
% function is what makes the per-panel try/catch cheap, and it means a panel can
% be re-styled without touching the layout code.

function visualizeBestDesign(results,xBest,cfg)
%VISUALIZEBESTDESIGN Post-run summary figures for a finished study.
%
% Reloads the winning design's cached band structure and draws:
%   1. a composite review figure (geometry + bands + convergence + summary
%      text + design-space slices),
%   2. standalone geometry, band-structure and progress figures.
%
% Nothing here talks to COMSOL. boomerangParams rebuilds exactly the P the
% winning evaluation used, and P.fileBase therefore resolves to the _bds.mat
% solveBands already wrote for it - which is the whole reason fileBase is set
% explicitly in boomerangParams rather than left to CreateFileBase.

P  = boomerangParams(xBest,cfg);
ds = loadCachedBandData(P,cfg);

% Re-score the winning gap with the SAME helper the objective used instead of
% re-deriving the formula here. A second copy of the fitness expression would
% be free to drift away from the one that actually drove the search, and a
% summary panel that disagrees with the optimizer is worse than no panel.
if ~isempty(ds) && isfield(ds,'full') && isstruct(ds.full)
    [~,detail] = gapFitness(ds.full,cfg);
else
    % No cached bands: still report what the results object knows, and label
    % the gap numbers as unavailable rather than inventing zeros.
    detail = struct('midGap',NaN,'gapSize',NaN,'gapRatio',NaN, ...
        'penalty',NaN,'fitness',-results.MinObjective, ...
        'status','band-data-unavailable');
end

fprintf('\n----- writing summary figures -----\n');

figs = gobjects(0);     % handles we created, so closeSummaryFigures can close
                        % exactly these and nothing else

% --- composite review figure --------------------------------------------
try
    figs(end+1) = buildCompositeFigure(results,xBest,P,ds,detail,cfg);
catch ME
    warning('bayesopt_boomerang:compositeFigure', ...
        'Could not build the composite figure: %s',ME.message);
end

% --- standalone geometry figure -----------------------------------------
% Single annotated cell next to the tiled lattice: the first shows what the
% four design variables mean, the second shows the ligament between holes,
% which is what actually opens or closes the gap.
try
    figG = newSummaryFigure(summaryFigureName('boomerang geometry',cfg),12,5.5);
    figs(end+1) = figG;
    % Default (loose) Padding here, not 'compact': both panels carry two-line
    % axes titles, and compact padding shrinks the strip reserved for the
    % layout title until it prints on top of them.
    tlG = tiledlayout(figG,1,2,'TileSpacing','compact');
    axG1 = nexttile(tlG);
    tryPanel(@() drawBoomerangCell(axG1,P,1,true),'geometry (single cell)');
    axG2 = nexttile(tlG);
    tryPanel(@() drawBoomerangCell(axG2,P,cfg.figNPeriods,false),'geometry (tiled)');
    % The geometry sketch itself is honest in every mode - polyshape draws the
    % real cell from the real dimensions, with or without COMSOL. It still gets
    % the banner, because the DESIGN it is showing is the winner of a synthetic
    % search, and that is the part someone would act on.
    [headG,fracG] = withDryRunHeadline(sprintf( ...
        'Boomerang unit cell  |  a = %.0f   r = %.0f   w = %.0f   th = %.0f nm', ...
        P.a*1e9,P.r*1e9,P.w*1e9,P.th*1e9),0.09,cfg);
    addFigureHeadline(figG,tlG,headG,fracG,dryRunHeadlineColor(cfg));
    saveFigurePair(figG,fullfile(cfg.datLoc,[cfg.filePrefix,'_geometry']),cfg);
catch ME
    warning('bayesopt_boomerang:geometryFigure', ...
        'Could not build the geometry figure: %s',ME.message);
end

% --- standalone band structure figure -----------------------------------
try
    figB = newSummaryFigure(summaryFigureName('boomerang best bands',cfg),8,6);
    figs(end+1) = figB;
    % 'Parent' name-value rather than axes(figB): the positional form was only
    % added in R2016a and this is the one place a plain single axes is wanted.
    axB = axes('Parent',figB);
    % No addFigureHeadline on this one - it is a bare axes, not a tiledlayout, so
    % there is no OuterPosition to shrink. plotBestBands carries the dry-run
    % banner in its axes title instead.
    tryPanel(@() plotBestBands(axB,ds,cfg,detail),'band structure');
    saveFigurePair(figB,fullfile(cfg.datLoc,[cfg.filePrefix,'_bestbands']),cfg);
catch ME
    warning('bayesopt_boomerang:bandFigure', ...
        'Could not build the band structure figure: %s',ME.message);
end

% --- standalone progress figure -----------------------------------------
% Convergence on top of the four design-space slices, because the two are read
% together: a flat convergence curve plus a variable pinned at its bound means
% the bound is setting the answer, not the physics.
try
    figP = newSummaryFigure(summaryFigureName('boomerang optimizer progress',cfg),13,7);
    figs(end+1) = figP;
    tlP = tiledlayout(figP,2,4,'TileSpacing','compact','Padding','compact');
    axP = nexttile(tlP,[1 4]);
    tryPanel(@() plotConvergenceTrace(axP,results),'convergence');
    drawDesignSpaceRow(tlP,results,cfg);
    [headP,fracP] = withDryRunHeadline(sprintf( ...
        'Optimizer progress  |  %d evaluations  |  best fitness = %.4f', ...
        numel(results.ObjectiveTrace),-results.MinObjective),0.07,cfg);
    addFigureHeadline(figP,tlP,headP,fracP,dryRunHeadlineColor(cfg));
    saveFigurePair(figP,fullfile(cfg.datLoc,[cfg.filePrefix,'_progress']),cfg);
catch ME
    warning('bayesopt_boomerang:progressFigure', ...
        'Could not build the progress figure: %s',ME.message);
end

% Close only the handles created above. Never "close all" - the bayesopt live
% plot windows are still open at this point and are worth keeping, which is the
% same reasoning behind closeNewFigures.
if cfg.closeSummaryFigures && ~isempty(figs)
    close(figs(isgraphics(figs)));
end
end

% -------------------------------------------------------------------------

function ds = loadCachedBandData(P,cfg)
%LOADCACHEDBANDDATA Reload the _bds.mat solveBands wrote for one design.
%
% Returns [] rather than throwing whenever the data cannot be had, so a missing
% or half-written cache degrades to a figure with an empty band panel instead of
% taking down the whole post-processing step. Reasons the file can legitimately
% be absent: the winning evaluation errored out, datLoc was moved mid-run, or
% the best point came from a resumed study whose data lives under another date
% folder.

ds = [];
matPath = fullfile(P.datLoc,[P.fileBase,'_bds.mat']);

if ~isfile(matPath)
    warning('bayesopt_boomerang:noBandCache', ...
        'No cached band data for the best design at %s',matPath);
    return
end

try
    S = load(matPath,'ds');
catch ME
    warning('bayesopt_boomerang:badBandCache', ...
        'Could not load %s: %s',matPath,ME.message);
    return
end

if ~isfield(S,'ds')
    warning('bayesopt_boomerang:badBandCache', ...
        '%s contains no variable ds',matPath);
    return
end

% Guard (3) again, on the figure path this time. The figures are the artifact most
% likely to outlive the folder they were written in, so band data reaching them
% is checked exactly as band data reaching the objective is. In a real run this
% throws rather than drawing a plausible-looking band diagram from fabricated
% frequencies - a figure is believed far more readily than a .mat file.
assertBandDataProvenance(S.ds,matPath,cfg);

ds = S.ds;
fprintf('  band data <- %s\n',matPath);
end

% -------------------------------------------------------------------------

function tryPanel(painterFcn,label)
%TRYPANEL Run one panel painter, downgrading any failure to a warning.
%
% Panels are independent, so a plot call that fails on one of them (an empty
% band matrix, a NaN-only trace, a graphics quirk on an older release) should
% cost that panel and nothing else. The axes is left blank, which is a visible
% and honest signal that something was missing.

try
    painterFcn();
catch ME
    warning('bayesopt_boomerang:panelFailed', ...
        'Could not draw the %s panel: %s',label,ME.message);
end
end

% -------------------------------------------------------------------------

function figH = newSummaryFigure(name,widthIn,heightIn)
%NEWSUMMARYFIGURE Create a figure sized for a deterministic PNG export.
%
% PaperPosition is set explicitly rather than using PaperPositionMode='auto'.
% 'auto' inherits the ON-SCREEN size, so the identical script would emit
% different pixel dimensions on a laptop and on a large monitor, and would
% silently shrink whenever the window got clipped by the screen. Pinning
% PaperPosition makes every export widthIn*cfg.figResolution by
% heightIn*cfg.figResolution pixels regardless of the machine it ran on.
% The 100 px/in on-screen size is only for previewing.

figH = figure('Name',name,'NumberTitle','off','Color','w', ...
    'Units','pixels','Position',[60 60 round(widthIn*100) round(heightIn*100)]);

% Force a light theme. From R2025a MATLAB themes the AXES background too, and
% figure('Color','w') does not override it - so on a machine set to dark mode
% these figures export as white-bordered black panels, and the near-black
% 'sym' band curves become invisible. Wrapped because the property does not
% exist on earlier releases, where light was the only behaviour anyway.
try
    figH.Theme = 'light';
catch
    % pre-R2025a: no themes, nothing to force
end

figH.PaperUnits    = 'inches';
figH.PaperPosition = [0 0 widthIn heightIn];
figH.PaperSize     = [widthIn heightIn];
end

% -------------------------------------------------------------------------

function figName = summaryFigureName(baseName,cfg)
%SUMMARYFIGURENAME Figure window title, marked when the data behind it is fake.
%
% The window title bar is where a figure is identified while it is still on
% screen, before anyone has looked at the file it came from - so it is worth
% marking. Returns baseName untouched outside a dry run.

figName = baseName;
if cfg.isDryRun
    figName = ['[DRY RUN - SYNTHETIC] ',baseName];
end
end

% -------------------------------------------------------------------------

function addFigureHeadline(figH,tl,headlineText,fracHeight,textColor)
%ADDFIGUREHEADLINE Suptitle above a tiledlayout, guaranteed not to overlap it.
%
% The obvious route - title(tl,...) - reserves exactly ONE line of space at the
% top of the layout and takes no account of the axes titles inside the top row.
% With axis-equal panels the axes box sits hard against the top of its tile, its
% own title spills up into that single line, and the two strings print on top of
% each other. Whether they collide then depends on a font-size coincidence.
%
% Shrinking the layout's OuterPosition and drawing the headline as a
% figure-level annotation instead makes the separation geometric: the layout
% physically cannot reach into the reserved strip. annotation() also predates
% tiledlayout by many releases, so nothing new is required.
%
% headlineText may be a char row or a cell array of lines.
%
% textColor is optional, and when it is absent or empty the Color property is not
% set AT ALL rather than being set to black. That distinction is not pedantic: an
% annotation textbox does not default to pure black (it inherits a dark grey), so
% passing [0 0 0] "for the default" visibly changes every headline. Caught by a
% pixel-diff of the exported PNGs against the pre-change version - 18 000 pixels
% differed on the composite, maximum channel difference 33 - which is exactly the
% class of silent change this option is not allowed to make. Omitting the property
% reproduces the previous output bit for bit.

colorArgs = {};
if nargin >= 5 && ~isempty(textColor)
    colorArgs = {'Color',textColor};
end

tl.OuterPosition = [0, 0, 1, 1-fracHeight];
annotation(figH,'textbox',[0, 1-fracHeight, 1, fracHeight], ...
    'String',headlineText, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'EdgeColor','none','FitBoxToText','off', ...
    'FontWeight','bold','FontSize',11,'Interpreter','none', ...
    colorArgs{:});
end

% -------------------------------------------------------------------------

function saveFigurePair(figH,pathNoExt,cfg)
%SAVEFIGUREPAIR Write one figure as .png and .fig, matching solveBands.
%
% solveBands saves each band diagram as both, so the summary figures follow
% suit: the .png for reports and email, the .fig so axes can be zoomed,
% re-scaled and re-styled months later without re-solving anything.
%
% print(...,'-dpng','-r<dpi>') is used instead of saveas because saveas always
% writes at screen resolution and ignores cfg.figResolution.
%
% pathNoExt already comes from fullfile at the call sites, so it carries the
% host platform's separators and no directory part needs adding here. Only the
% extension is appended, which is a suffix rather than a path join and so is one
% of the few places plain concatenation is correct.

try
    print(figH,[pathNoExt,'.png'],'-dpng',['-r',num2str(cfg.figResolution)]);
    savefig(figH,[pathNoExt,'.fig']);
    fprintf('  figure    -> %s.png / .fig\n',pathNoExt);
catch ME
    warning('bayesopt_boomerang:figureSave', ...
        'Could not save %s: %s',pathNoExt,ME.message);
end
end

% -------------------------------------------------------------------------

function figH = buildCompositeFigure(results,xBest,P,ds,detail,cfg)
%BUILDCOMPOSITEFIGURE One-page review of the whole study.
%
% LAYOUT - 3 rows x 4 columns of tiles:
%   row 1   geometry (2 tiles)      band structure (2 tiles)
%   row 2   convergence (2 tiles)   summary text (2 tiles)
%   row 3   fitness vs a | vs r | vs w | vs th
%
% Why 3x4 and not the Python script's 2x3: this study has five topics, not six,
% and one of them (the design-space exploration) is naturally four small axes
% sharing a y-axis. A 3x4 grid gives every large panel two full tiles - so the
% band structure and the geometry are both big enough to read - and lets the
% four slices sit in their own row without a nested tiledlayout. Nested layouts
% would work but add an API dependency for no visual gain here.
%
% Panels are painted through tryPanel so one bad panel cannot cost the figure.

% 16x10 rather than 16x9: with three rows, 9 in leaves the text panel too short
% for its ~20 lines at a legible font size.
figH = newSummaryFigure(summaryFigureName('bayesopt boomerang summary',cfg),16,10);
tl = tiledlayout(figH,3,4,'TileSpacing','compact','Padding','compact');

% Row 1: what was designed, and how it behaves.
% doAnnotate is FALSE here: at cfg.figNPeriods the dimension labels are drawn on
% one cell out of nine and collide with each other in a half-width tile. The
% numbers they would carry are already in the summary panel and the suptitle, so
% this panel's job is purely the shape and the inter-hole ligament. The
% standalone geometry figure carries the annotated single cell.
axGeo = nexttile(tl,[1 2]);
tryPanel(@() drawBoomerangCell(axGeo,P,cfg.figNPeriods,false),'geometry');

axBnd = nexttile(tl,[1 2]);
tryPanel(@() plotBestBands(axBnd,ds,cfg,detail),'band structure');

% Row 2: how the optimizer got there, and the numbers in text form.
axCon = nexttile(tl,[1 2]);
tryPanel(@() plotConvergenceTrace(axCon,results),'convergence');

axSum = nexttile(tl,[1 2]);
tryPanel(@() writeSummaryPanel(axSum,xBest,detail,cfg,results),'summary text');

% Row 3: where in the design space the solves were actually spent.
drawDesignSpaceRow(tl,results,cfg);

% suptitle equivalent, carrying the headline metrics the way the Python script's
% fig.suptitle does - so the figure is self-describing once it has left the
% folder it was written into.
[headC,fracC] = withDryRunHeadline({ ...
    sprintf('Best boomerang unit cell   |   fitness = %.4f   (objective = %.4f)', ...
        -results.MinObjective,results.MinObjective), ...
    sprintf(['a=%d  r=%d  w=%d  th=%d nm   |   mid-gap = %.3f GHz   ' ...
             'gap ratio = %.2f%%   target = %.2f GHz   |   %d evaluations   |   %s'], ...
        xBest.a,xBest.r,xBest.w,round(cfg.th*1e9), ...
        detail.midGap/1e9,detail.gapRatio*100,cfg.targetFreq/1e9, ...
        numel(results.ObjectiveTrace),detail.status)},0.075,cfg);
addFigureHeadline(figH,tl,headC,fracC,dryRunHeadlineColor(cfg));

saveFigurePair(figH,fullfile(cfg.datLoc,[cfg.filePrefix,'_summary']),cfg);
end

% -------------------------------------------------------------------------

function drawDesignSpaceRow(tl,results,cfg)
%DRAWDESIGNSPACEROW Fill the next row of tiles with one slice per variable.
%
% Variable names are read off XTrace rather than hard-coded so this keeps
% working if the search space is ever extended or reordered. Capped at four
% because that is what the callers' layouts reserve; a fifth variable would
% otherwise silently overflow into a tile another panel owns.

try
    varNames = results.XTrace.Properties.VariableNames;
catch
    warning('bayesopt_boomerang:noXTrace', ...
        'results.XTrace is not a table - skipping the design-space slices');
    return
end

nSlice = min(numel(varNames),4);
for k = 1:nSlice
    axK = nexttile(tl);
    % varNames{k} is captured by value when the handle is created, so each
    % iteration gets its own variable name.
    tryPanel(@() plotDesignSlice(axK,results,cfg,varNames{k}), ...
        ['design slice ',varNames{k}]);
end
end

% -------------------------------------------------------------------------

function drawBoomerangCell(ax,P,nPeriods,doAnnotate)
%DRAWBOOMERANGCELL Top view of the boomerang cell, drawn without COMSOL.
%
% Pure MATLAB on purpose: the point of this panel is to check at a glance that
% the winning design is the shape you think it is, and launching a COMSOL
% geometry build to answer that would make the panel unusable after the
% LiveLink session ends - and would fail on any machine without a licence.
%
% GEOMETRY (transcribed from buildBoomerangUnitCell.m, which is the builder
% runBands_2D actually calls for celltype 'boomerang'):
%   * cell    - primitive rhombus of the hexagonal lattice, side a, 60 deg
%               interior angle, vertices (0,0) (a/2,a*sqrt(3)/2)
%               (3a/2,a*sqrt(3)/2) (a,0). `a` alone sets the footprint.
%   * centre  - the rhombus centroid, (3a/4, a*sqrt(3)/4).
%   * hole    - three-pointed star: three rectangles, each w wide and r long,
%               radiating from the centre at theta = 90, 210, 330 deg. Each arm
%               spans radially from the centre out to r, matching the COMSOL
%               rectangle centres at (3a/4, a*sqrt(3)/4 + r/2) and the two
%               120 deg rotations of it.
%   * th      - full slab thickness in z, so it is INVISIBLE in a top view. It
%               is reported in the title instead of drawn.
%
% FILLETS (P.r1 inner corners, P.r2 outer tips) are NOT drawn. At the fixed
% 10 nm they are under 2% of the smallest arm width and would be invisible at
% this scale, while reproducing COMSOL's fillet-vertex selection in polyshape
% would be guesswork - a wrong fillet is worse than an honest sharp corner.
% The title says so, so the sketch is never mistaken for the meshed geometry.

nm = 1e9;
a  = P.a*nm;
r  = P.r*nm;
w  = P.w*nm;
th = P.th*nm;

% The closing repeat of (0,0) that the COMSOL Polygon table carries is dropped:
% polyshape closes the ring itself and warns about duplicate vertices.
cellPoly = polyshape([0, a/2, 3*a/2, a], [0, a*sqrt(3)/2, a*sqrt(3)/2, 0]);

cx = 3*a/4;                 % cell centre, and the star's origin
cy = a*sqrt(3)/4;

% One arm, pointing at theta = 90 deg (along +y) from the centre, then rotated
% by +120 and +240 deg about the centre for the other two. rotate() takes the
% angle in DEGREES, counter-clockwise, about the given reference point.
arm  = polyshape([-w/2, w/2, w/2, -w/2] + cx, [0, 0, r, r] + cy);
hole = union(union(arm,rotate(arm,120,[cx cy])),rotate(arm,240,[cx cy]));

% subtract leaves an interior ring, so the star renders as a genuine void
% rather than as a white shape painted over the slab.
solid = subtract(cellPoly,hole);

% Lattice vectors: the two rhombus edges leaving the origin.
a1 = [a, 0];
a2 = [a/2, a*sqrt(3)/2];

hold(ax,'on');

% Neighbour cells first, faded, so the reference cell reads as the foreground
% and the ligament between adjacent holes - the feature that sets the gap - is
% visible without measuring anything.
% ii/jj rather than i/j so the imaginary unit is not shadowed.
for ii = 0:nPeriods-1
    for jj = 0:nPeriods-1
        if ii == 0 && jj == 0
            continue    % reference cell is drawn last, on top
        end
        plot(ax,translate(solid,ii*a1 + jj*a2), ...
            'FaceColor',[0.80 0.86 0.93],'FaceAlpha',0.70, ...
            'EdgeColor',[0.55 0.62 0.72],'LineWidth',0.6);
    end
end
plot(ax,solid, ...
    'FaceColor',[0.64 0.79 0.94],'FaceAlpha',1.0, ...
    'EdgeColor',[0.10 0.32 0.60],'LineWidth',1.5);

if doAnnotate
    % Lattice vectors as true data-space arrows. quiver with an explicit scale
    % of 0 disables its autoscaling, which would otherwise shrink both arrows
    % to fit an invisible grid and make them lie about the lattice.
    quiver(ax,[0 0],[0 0],[a1(1) a2(1)],[a1(2) a2(2)],0, ...
        'Color',[0.15 0.15 0.15],'LineWidth',1.2,'MaxHeadSize',0.18);
    text(ax,a/2,-0.05*a,sprintf('a_1 = a = %.0f nm',a), ...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8);
    text(ax,a/4-0.03*a,a*sqrt(3)/4,'a_2', ...
        'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',8);

    % r and w are measured on the theta = 90 deg arm. Both indicator lines fall
    % inside the void, so they sit on the axes background and stay legible.
    plot(ax,cx,cy,'+','Color',[0.75 0.15 0.15],'MarkerSize',8,'LineWidth',1.0);
    plot(ax,[cx cx],[cy cy+r],'-','Color',[0.75 0.15 0.15],'LineWidth',1.2);
    text(ax,cx+0.03*a,cy+0.45*r,sprintf('r = %.0f nm',r), ...
        'Color',[0.75 0.15 0.15],'FontSize',8,'VerticalAlignment','middle');
    plot(ax,[cx-w/2 cx+w/2],(cy+0.80*r)*[1 1],'-', ...
        'Color',[0.10 0.45 0.15],'LineWidth',1.2);
    text(ax,cx+w/2+0.02*a,cy+0.80*r,sprintf('w = %.0f nm',w), ...
        'Color',[0.10 0.45 0.15],'FontSize',8,'VerticalAlignment','middle');
end

% Limits computed from the tiling rather than left to autoscale, so the aspect
% ratio stays honest and the annotations outside the cell are not clipped:
% x spans 0..1.5*a*nPeriods and y spans 0..(sqrt(3)/2)*a*nPeriods.
pad = 0.18*a;
axis(ax,'equal');
xlim(ax,[-pad, 1.5*a*nPeriods + pad]);
ylim(ax,[-pad, (sqrt(3)/2)*a*nPeriods + pad]);
grid(ax,'on');
box(ax,'on');
xlabel(ax,'x [nm]','FontSize',9);
ylabel(ax,'y [nm]','FontSize',9);

% Deliberately terse. With axis equal the axes box grows upward inside its tile,
% so the title ends up level with the enclosing tiledlayout's own centred title;
% the collision is then in x, not y, and no Padding setting fixes it. Keeping the
% title short enough that it cannot reach the middle of the figure does.
title(ax,sprintf('Top view, %dx%d periods',nPeriods,nPeriods),'FontSize',9);

% The out-of-plane and fillet caveats live inside the axes rather than in the
% title: they stay attached to the drawing they qualify, and they compete for no
% title space. They have to be stated somewhere - this sketch is a sketch, not
% the meshed geometry.
text(ax,0.98,0.03,{sprintf('th = %.0f nm out of plane',th), ...
    sprintf('fillets r1/r2 = %.0f/%.0f nm not drawn',P.r1*1e9,P.r2*1e9)}, ...
    'Units','normalized', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom', ...
    'FontSize',7,'Color',[0.40 0.40 0.40]);

set(ax,'FontSize',8);
end

% -------------------------------------------------------------------------

function plotBestBands(ax,ds,cfg,detail)
%PLOTBESTBANDS Band structure of the best design, complete gaps shaded.
%
% Both z-symmetry sectors are drawn, because with TwoSymPlanes=0 and
% solveasym=1 solveBands returns the spectrum split in two - ds.sym (even about
% z) and ds.asym (odd about z) - and a gap is only a COMPLETE gap if it is
% empty in BOTH. Plotting one sector alone would show gaps that the other
% sector fills in, which is exactly the mistake the objective avoids by scoring
% ds.full. Solid = even, dashed = odd; the shaded bands are ds.full, i.e. the
% complete gaps across both sectors.

hold(ax,'on');
legH   = gobjects(0);
legTxt = {};
fMax   = 0;

% Complete gaps FIRST: within one axes MATLAB stacks later children above
% earlier ones, so shading drawn now sits behind the dispersion curves instead
% of hiding them.
if isfield(ds,'full') && isstruct(ds.full) && isfield(ds.full,'gapSize') ...
        && ~isempty(ds.full.gapSize) && isfield(ds.full,'midGap')
    mid = ds.full.midGap(:);        % (:) - these are stored as ROW vectors
    gsz = ds.full.gapSize(:);
    nGap = min(numel(mid),numel(gsz));
    hGap = gobjects(0);
    for k = 1:nGap
        % Vertex order [0 3 3 0] x [top top bottom bottom] traces the band
        % without self-intersecting, matching the patches solveBands draws.
        % hGap is deliberately overwritten: only one handle is wanted so the
        % legend carries a single "complete gap" entry, not one per gap.
        yBand = 1e-9*(mid(k) + 0.5*[gsz(k) gsz(k) -gsz(k) -gsz(k)]);
        hGap = patch(ax,'XData',[0 3 3 0],'YData',yBand, ...
            'FaceColor',[1.00 0.82 0.40],'FaceAlpha',0.45, ...
            'EdgeColor',[0.80 0.55 0.10],'LineStyle','--','LineWidth',0.6);
    end
    if ~isempty(hGap)
        legH(end+1)   = hGap;
        legTxt{end+1} = sprintf('complete gap (%d)',nGap);
    end
end

% Even sector (sym).
if isfield(ds,'sym')
    [h,f1] = plotBandSector(ax,ds.sym,'-',[0.10 0.10 0.10],1.3);
    if ~isempty(h)
        legH(end+1)   = h;
        legTxt{end+1} = 'even about z (sym)';
    end
    fMax = max(fMax,f1);
end

% Odd sector (asym).
if isfield(ds,'asym')
    [h,f2] = plotBandSector(ax,ds.asym,'--',[0.20 0.40 0.75],1.1);
    if ~isempty(h)
        legH(end+1)   = h;
        legTxt{end+1} = 'odd about z (asym)';
    end
    fMax = max(fMax,f2);
end

% Target frequency. A plotted line rather than yline(): yline arrived in
% R2018b and its axes-first form later still, and this file has to run wherever
% the sweeps run.
hTgt = plot(ax,[0 3],cfg.targetFreq/1e9*[1 1],':', ...
    'Color',[0.75 0.15 0.15],'LineWidth',1.6);
legH(end+1)   = hTgt;
legTxt{end+1} = sprintf('target %.2f GHz',cfg.targetFreq/1e9);

% The winning mid-gap, so the panel shows WHICH gap the fitness scored - with
% several complete gaps shaded it is otherwise ambiguous.
if isfield(detail,'midGap') && isfinite(detail.midGap)
    hMid = plot(ax,[0 3],detail.midGap/1e9*[1 1],'-.', ...
        'Color',[0.10 0.55 0.20],'LineWidth',1.3);
    legH(end+1)   = hMid;
    legTxt{end+1} = sprintf('scored mid-gap %.2f GHz',detail.midGap/1e9);
end

xlim(ax,[0 3]);

% Choose the upper y-limit from the bands AND from every horizontal marker
% drawn above, not from the bands alone. Scaling to the bands only is the
% obvious thing to do and it is wrong: whenever the solved bands top out below
% the target - routine at large `a`, where nbands=15 may not reach 13 GHz - the
% target line is drawn outside the axes and silently clipped. The marker then
% contributes nothing precisely in the case it is most needed, which is telling
% you the design is nowhere near target. Better to squash the bands a little
% and keep every reference line on screen.
if fMax > 0
    fTop = fMax;
    fTop = max(fTop, cfg.targetFreq/1e9);
    if isfield(detail,'midGap') && isfinite(detail.midGap)
        fTop = max(fTop, detail.midGap/1e9);
    end
    ylim(ax,[0 1.05*fTop]);

    % If the bands really do stop well short of the target, say so in the
    % panel: a squashed band plot plus a lone line near the top is otherwise
    % easy to misread as a solved gap sitting at the target.
    if cfg.targetFreq/1e9 > 1.02*fMax
        % Sits BELOW the target line, not at the top of the axes. Whenever this
        % note appears the target line is by construction the topmost thing on
        % the plot (it is what set the limit), so anchoring the text at the top
        % would place it straight through the line it is describing.
        text(ax,0.03,0.88, ...
            sprintf('bands reach only %.2f GHz - target is above the solved range', fMax), ...
            'Units','normalized','VerticalAlignment','top', ...
            'FontSize',7.5,'Color',[0.75 0.15 0.15]);
    end
end

% k_norm runs 0->3 over Gamma-M-K-Gamma for a hexagonal cell; see the kx/ky
% expressions and their comments in runBands_2D.m. Note that solveBands labels
% the same axis Gamma-X-M-Gamma, which are the SQUARE-lattice point names - the
% labels below follow runBands_2D, which is what actually set the path.
set(ax,'XTick',0:3,'XTickLabel',{'\Gamma','M','K','\Gamma'});
grid(ax,'on');       % the x grid lines double as high-symmetry-point markers
box(ax,'on');
xlabel(ax,'wavevector along \Gamma-M-K-\Gamma','FontSize',9);
ylabel(ax,'frequency [GHz]','FontSize',9);

if isfield(detail,'gapRatio') && isfinite(detail.gapRatio)
    ttl2 = sprintf('scored gap: %.3f GHz wide, ratio %.2f%%  (%s)', ...
        detail.gapSize/1e9,detail.gapRatio*100,detail.status);
else
    ttl2 = 'no complete gap scored';
end
% The banner goes in the TITLE of this panel rather than a figure annotation,
% because plotBestBands is used both standalone (bare axes, nothing to annotate
% around) and as one tile of the composite. A title line travels with the axes in
% both cases. This is the one panel where the curves themselves are fabricated,
% so it is the one that most needs saying so - and in the standalone figure it is
% the ONLY marking inside the image, since a bare axes has no headline. Hence the
% colour, which is applied only in a dry run: title() would otherwise inherit the
% axes TitleColor, and passing an explicit black would change a real run's figure.
ttlArgs = {'FontSize',9};
if cfg.isDryRun
    ttlArgs = [ttlArgs,{'Color',dryRunHeadlineColor(cfg)}];
end
title(ax,[dryRunBannerText(cfg), ...
    {'Mechanical band structure of the best design',ttl2}],ttlArgs{:});

if ~isempty(legH)
    % legend(subset,labels) attaches to the axes owning the handles, and
    % restricting it to a subset keeps the 15 band lines per sector out of it.
    % 'best' because the dispersion curves cover most of the axes and where the
    % free space is depends entirely on the design.
    legend(legH,legTxt,'Location','best','FontSize',7);
end
set(ax,'FontSize',8);
end

% -------------------------------------------------------------------------

function [hFirst,fMax] = plotBandSector(ax,sector,lineStyle,lineColor,lineWidth)
%PLOTBANDSECTOR Plot every band of one symmetry sector.
%
% Returns the first line handle only, for a one-entry-per-sector legend, plus
% the highest finite frequency so the caller can set a y-limit that does not
% depend on a NaN. Everything is guarded: a sector struct can be empty (the
% stub solveBands returns on a cache hit) or truncated (a solve killed
% mid-sweep), and neither should abort the figure.

hFirst = gobjects(0);
fMax   = 0;

if ~isstruct(sector) || ~isfield(sector,'F') || ~isfield(sector,'k_norm') ...
        || isempty(sector.F) || isempty(sector.k_norm)
    return
end

k = sector.k_norm(:);
F = sector.F*1e-9;              % Hz -> GHz

% Defensive: runBands_2D appends the wrapped Gamma point to both F and k_norm
% at the end, so they normally match. A solve interrupted between those two
% appends would leave them one apart, and plot() would error out on the whole
% panel instead of just losing a point.
if size(F,1) ~= numel(k)
    n = min(size(F,1),numel(k));
    k = k(1:n);
    F = F(1:n,:);
end

hAll   = plot(ax,k,F,lineStyle,'Color',lineColor,'LineWidth',lineWidth);
hFirst = hAll(1);

fFinite = F(isfinite(F));
if ~isempty(fFinite)
    fMax = max(fFinite);
end
end

% -------------------------------------------------------------------------

function plotConvergenceTrace(ax,results)
%PLOTCONVERGENCETRACE Best-so-far and per-evaluation fitness vs evaluation.
%
% FITNESS (= -objective) is plotted, not the objective bayesopt minimizes, so
% that "up is better" reads the way everyone expects from a bandgap figure of
% merit. The axis label says so explicitly, because the sign convention is the
% one thing about this plot that can be misread.
%
% Two series, following _plot_progress in run_opt_comsol.py: the per-evaluation
% points show where the acquisition function chose to look (including the
% deliberate excursions into bad regions that a GP must make), and the
% non-decreasing best-so-far line shows actual progress. Reading only the
% points makes a healthy run look erratic.

objTrace = results.ObjectiveTrace(:);           % (:) - do not assume column
minTrace = results.ObjectiveMinimumTrace(:);
nEval    = numel(objTrace);
idx      = (1:nEval)';

hold(ax,'on');
legH   = gobjects(0);
legTxt = {};

hEval = plot(ax,idx,-objTrace,'o','LineStyle','none','MarkerSize',5, ...
    'MarkerFaceColor',[0.55 0.68 0.85],'MarkerEdgeColor',[0.20 0.32 0.52]);
legH(end+1)   = hEval;
legTxt{end+1} = 'each evaluation';

% Trimmed to a common length: the two traces are the same length for a
% completed study, but a study stopped by an output function can leave them one
% apart, and that must not cost the panel.
m = min(numel(minTrace),nEval);
if m > 0
    hBest = plot(ax,idx(1:m),-minTrace(1:m),'-','LineWidth',1.9, ...
        'Color',[0.85 0.33 0.10]);
    legH(end+1)   = hBest;
    legTxt{end+1} = 'best so far';
end

% Errored evaluations. boomerangObjective returns NaN for a solver failure
% rather than 0, so these points are absent from the scatter above; marking
% them on the zero line makes the difference between "no gap found" (a real
% fitness of 0) and "the solve failed" visible.
iErr = find(~isfinite(objTrace));
if ~isempty(iErr)
    hErr = plot(ax,iErr,zeros(numel(iErr),1),'x','LineStyle','none', ...
        'MarkerSize',8,'LineWidth',1.4,'Color',[0.75 0.15 0.15]);
    legH(end+1)   = hErr;
    legTxt{end+1} = 'solver error (NaN)';
end

if nEval > 0
    xlim(ax,[0.5, nEval+0.5]);
end
grid(ax,'on');
box(ax,'on');
xlabel(ax,'COMSOL evaluation #','FontSize',9);
ylabel(ax,'fitness = -objective  (higher is better)','FontSize',9);
title(ax,sprintf('Convergence: best fitness %.4f after %d evaluations', ...
    -results.MinObjective,nEval),'FontSize',9);
if ~isempty(legH)
    % 'best' rather than a fixed corner: the best-so-far curve climbs into the
    % upper right, but how much of the axes the scatter occupies varies per run.
    legend(legH,legTxt,'Location','best','FontSize',7);
end
set(ax,'FontSize',8);
end

% -------------------------------------------------------------------------

function writeSummaryPanel(ax,xBest,detail,cfg,results)
%WRITESUMMARYPANEL Monospaced text block, the analogue of ax_info in Python.
%
% Text in a fixed-width font rather than a uitable or an annotation box: the
% point is one at-a-glance record of the run that survives being screenshotted
% into a lab notebook or pasted into an email, and monospace is what keeps the
% columns of numbers lined up when it gets there.

axis(ax,'off');

nEval = numel(results.ObjectiveTrace);

% Distance from the target, which is what the Gaussian penalty is a function
% of. Reported separately because a good gap ratio at the wrong frequency and a
% poor gap ratio at the right one produce similar fitness values, and only this
% number tells them apart.
fOff = abs(cfg.targetFreq - detail.midGap)/1e9;

% A cell array of char with SCALAR x/y makes ONE multi-line text object; a
% vector x/y would instead make one object per line, all stacked at the origin.
%
% Kept to ~20 lines on purpose. The panel owns half a row of a three-row
% layout, roughly 3 in tall, and a text block taller than that is silently
% clipped by the tile below it rather than being scaled to fit - so lines are
% packed (TARGET and STUDY use trailing continuations) instead of being spread
% out with blank separators.
lines = { ...
    'BEST DESIGN   (integer nm, as searched)', ...
    sprintf('  a       = %6d nm   cell side / lattice constant',xBest.a), ...
    sprintf('  r       = %6d nm   arm length, centre to tip',xBest.r), ...
    sprintf('  w       = %6d nm   arm width (narrowest feature)',xBest.w), ...
    sprintf('  th      = %6d nm   slab thickness (FIXED, not searched)', ...
        round(cfg.th*1e9)), ...
    sprintf('  r1 / r2 = %3.0f / %.0f nm   fillets (held fixed)', ...
        cfg.r1*1e9,cfg.r2*1e9), ...
    '', ...
    sprintf('COMPLETE GAP   (%s)',detail.status), ...
    sprintf('  mid-gap   = %9.3f GHz',detail.midGap/1e9), ...
    sprintf('  gap size  = %9.3f GHz',detail.gapSize/1e9), ...
    sprintf('  gap ratio = %9.3f %%    gapSize/midGap',detail.gapRatio*100), ...
    sprintf('  penalty   = %9.4f      exp(-((f_t-f_mid)/sigma)^2)',detail.penalty), ...
    sprintf('  fitness   = %9.4f      gapRatio * penalty',detail.fitness), ...
    sprintf('  objective = %9.4f      -fitness (bayesopt minimizes)', ...
        results.MinObjective), ...
    '', ...
    sprintf('TARGET    f = %.3f GHz   sigma = %.3f GHz', ...
        cfg.targetFreq/1e9,cfg.sigma/1e9), ...
    sprintf('          |f_mid - f_target| = %.3f GHz  (frequency miss)',fOff), ...
    '', ...
    sprintf('STUDY     %d of %d evaluations, %d seed points', ...
        nEval,cfg.maxEvaluations,cfg.numSeedPoints), ...
    sprintf('          %d k-points, %d bands, mesh %d, max %g DOF', ...
        cfg.kpts,cfg.nbands,cfg.meshSize,cfg.maxDof), ...
    sprintf('          %s',cfg.datLoc)};

% The dry-run banner goes at the TOP, where it is read first - and the block has
% to get SHORTER, not longer, to fit it. Two effects were measured on the rendered
% PNG rather than guessed at, and both push the same way:
%   * the 21 lines above already fill this tile exactly, so a 22nd is not scaled
%     to fit, it is clipped by the tile below - and what it lands on is the row of
%     design-slice titles;
%   * the dry-run headline is one line taller, which shrinks the whole
%     tiledlayout and costs this tile roughly one further line.
% So all three blank separators are dropped in dry-run mode, buying three lines
% for the cost of one, which leaves the block with real slack instead of sitting
% on the limit again. The section headers are unindented and their contents are
% indented, so the structure survives losing the blank lines.
% If you add a line here, re-render and look at the PNG.
if cfg.isDryRun
    lines = [dryRunBannerText(cfg), lines(~strcmp(lines,''))];
end

% Interpreter 'none': the status strings, field names and folder paths are full
% of underscores, which the default TeX interpreter would turn into subscripts.
text(ax,0.02,0.98,lines,'Units','normalized', ...
    'FontName','FixedWidth','FontSize',8, ...
    'VerticalAlignment','top','HorizontalAlignment','left', ...
    'Interpreter','none', ...
    'BackgroundColor',[0.95 0.96 1.00],'EdgeColor',[0.72 0.76 0.86], ...
    'Margin',5);
end

% -------------------------------------------------------------------------

function plotDesignSlice(ax,results,cfg,varName)
%PLOTDESIGNSLICE Fitness against one design variable, over every evaluation.
%
% The analogue of _scatter_gap_tw in run_opt_comsol.py. Nothing is interpolated
% and no surface is fitted: the samples sit wherever the acquisition function
% chose to look, so each marker is one real COMSOL solve and the scatter is the
% honest projection of a 4-D sample set onto one axis.
%
% What this panel is FOR is spotting a variable pinned against its bound. If the
% best point sits on a bound then the BOUND, not the physics, is setting the
% answer, and the range should be widened before trusting the result - so the
% bounds are drawn as guide lines and the title flags the pinning explicitly.

fit = -results.ObjectiveTrace;
fit = fit(:);

if ~ismember(varName,results.XTrace.Properties.VariableNames)
    return
end
v = double(results.XTrace.(varName));
v = v(:);

n   = min(numel(v),numel(fit));     % guard against a truncated trace
v   = v(1:n);
fit = fit(1:n);

hold(ax,'on');
% MarkerFaceColor is set explicitly instead of passing the 'filled' flag: the
% flag sets MarkerFaceColor to 'flat', so the two together depend on which is
% applied last. Semi-transparent faces matter here because repeated designs
% (cache hits) land on exactly the same point.
scatter(ax,v,fit,26, ...
    'MarkerFaceColor',[0.29 0.45 0.69],'MarkerFaceAlpha',0.65, ...
    'MarkerEdgeColor','none');

% max() skips NaN, so an errored evaluation can never be reported as the best.
[fitBest,iBest] = max(fit);
haveBest = ~isempty(iBest) && isfinite(fitBest);
if haveBest
    plot(ax,v(iBest),fitBest,'p','LineStyle','none','MarkerSize',13, ...
        'MarkerFaceColor',[0.90 0.33 0.13],'MarkerEdgeColor','k');
end

% Bound guides drawn after the data, so the y-range has already settled; ylim is
% then re-applied because the vertical lines would otherwise be autoscaled over.
pinTxt = '';
yl = ylim(ax);
if isfield(cfg.bounds,varName)
    b = double(cfg.bounds.(varName));
    for k = 1:numel(b)
        plot(ax,[b(k) b(k)],yl,':','Color',[0.55 0.55 0.55],'LineWidth',1.0);
    end
    if numel(b) >= 2
        pad = max(0.05*(b(2)-b(1)),1);     % max(...,1) keeps xlim non-degenerate
        xlim(ax,[b(1)-pad, b(2)+pad]);     % if a bound range is ever collapsed
        % 0.5 nm tolerance: the variables are integers in nm, so anything closer
        % than half a step to a bound IS on it.
        if haveBest && (abs(v(iBest)-b(1)) < 0.5 || abs(v(iBest)-b(2)) < 0.5)
            pinTxt = '  [at bound]';
        end
    end
end
ylim(ax,yl);

grid(ax,'on');
box(ax,'on');
% Interpreter 'none' on anything carrying varName: the names happen to be
% single letters today, but a future variable with an underscore would
% otherwise be silently rendered as a TeX subscript.
xlabel(ax,sprintf('%s [nm]',varName),'FontSize',9,'Interpreter','none');
ylabel(ax,'fitness','FontSize',9);
if haveBest
    title(ax,sprintf('%s: best %g nm%s',varName,v(iBest),pinTxt), ...
        'FontSize',9,'Interpreter','none');
else
    title(ax,sprintf('%s: no finite fitness',varName), ...
        'FontSize',9,'Interpreter','none');
end
set(ax,'FontSize',8);
end
