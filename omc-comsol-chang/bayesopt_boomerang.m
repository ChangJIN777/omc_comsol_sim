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
% REQUIREMENTS
%   * Statistics and Machine Learning Toolbox (bayesopt, optimizableVariable)
%   * COMSOL with LiveLink for MATLAB, already running, before this script is
%     started - solveBands -> runBands_2D calls the COMSOL Java API.
%
% USAGE
%   >> run('bayesopt_boomerang.m')
% then, to spend more solves on the same surrogate:
%   >> results = resume(results,'MaxObjectiveEvaluations',20);
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

% --- objective definition -------------------------------------------------
cfg.targetFreq = 13e9;      % target mechanical mid-gap frequency [Hz]
cfg.sigma      = 5e9;       % width of the Gaussian frequency penalty [Hz]

% --- fabrication tolerance ------------------------------------------------
% Minimum feature size, applied to the hole arm width w - the narrowest
% etched feature in the tri-arm geometry. See boomerangFabConstraint for why
% the `a - r` term from boomerang_optimize_sweep_diamond.m is deliberately not
% reproduced, and for the separate r < sqrt(3)*a/4 in-cell guard.
% Held in integer nm so the comparison is exact. Storing it in metres and
% comparing w*1e-9 >= 50e-9 would be a round-off coin flip on designs sitting
% exactly on the limit, e.g. w = 50 nm, because the computed product and the
% parsed literal need not agree to the last bit.
cfg.minFeatureNm = 50;      % [nm]

% --- search space, in integer nm (ranges from sweep_boomerang_code.m) -----
cfg.bounds.a  = [600 1000];
cfg.bounds.r  = [150  250];   % sweep used r in [0.6*250, 250] nm
cfg.bounds.w  = [ 50  120];
cfg.bounds.th = [275  350];

% --- fixed geometry ------------------------------------------------------
cfg.r1 = 10e-9;             % fillet radius at the INNER corners, where the
                            % three arms meet near the cell centre
cfg.r2 = 10e-9;             % fillet radius at the OUTER arm tips

% --- solver fidelity (kept identical to sweep_boomerang_code.m) ----------
cfg.kpts     = 9;           % k-points EXCLUDING gamma
cfg.nbands   = 15;
cfg.meshSize = 4;
cfg.maxDof   = 3e6;

% --- optimizer budget ---------------------------------------------------
cfg.maxEvaluations = 40;    % total COMSOL solves bayesopt may spend
cfg.numSeedPoints  = 8;     % random design points before the GP takes over

% --- I/O ---------------------------------------------------------------
% NOTE the trailing '\' is required. solveBands normalises a separator into a
% local variable but then writes the .mat with the raw P.datLoc, so a missing
% separator silently drops output into the parent folder.
currentDate = datestr(now,'mmddyyyy');                                      %#ok<TNOW1,DATST>
cfg.datLoc  = ['.\test\boomerang_bayesopt\',currentDate,'\'];
cfg.logPath   = [cfg.datLoc,'bayesopt_boomerang_log.txt'];
cfg.statePath = [cfg.datLoc,'bayesopt_boomerang_state.mat'];

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
initIterationLog(cfg.logPath);

%% ------------------------------------------------------------------------
%  Optimization variables
%  ------------------------------------------------------------------------
optVars = [ ...
    optimizableVariable('a', cfg.bounds.a, 'Type','integer'), ...
    optimizableVariable('r', cfg.bounds.r, 'Type','integer'), ...
    optimizableVariable('w', cfg.bounds.w, 'Type','integer'), ...
    optimizableVariable('th',cfg.bounds.th,'Type','integer')];

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
fprintf('  th = %d nm\n', xBest.th);
fprintf('  objective = %.6f  (fitness = %.6f)\n', ...
    results.MinObjective, -results.MinObjective);
fprintf('  iteration log: %s\n', cfg.logPath);

save([cfg.datLoc,'bayesopt_boomerang_results.mat'],'results','cfg','xBest');

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
% Two separate tests:
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
% (2) Arms must stay inside the cell. The cell centre sits sqrt(3)*a/4 from
%     the nearest cell edge, and an arm reaches r from the centre, so r must
%     stay below that or the hole punches through the boundary and the
%     geometry build produces a shape the periodic BCs no longer describe.
%     The current bounds satisfy this with only ~10 nm to spare at
%     a = 600 nm, r = 250 nm, so the guard matters as soon as the ranges
%     widen.
holeWidthNm = x.w;                      % narrowest etched feature
armReachNm  = x.r;                      % centre to arm tip
cellEdgeNm  = sqrt(3)*double(x.a)/4;    % centre to nearest cell edge

tf = (holeWidthNm >= cfg.minFeatureNm) & (armReachNm < cellEdgeNm);
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

P = boomerangParams(x,cfg);

fprintf('\n----- evaluating a=%dnm r=%dnm w=%dnm th=%dnm -----\n', ...
    x.a, x.r, x.w, x.th);

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
% settings are identical to the existing sweeps; only a/r/w/th now come from
% the optimizer instead of being hard-coded.

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
P.th = double(x.th)*nm;     % full slab thickness in z
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
end

% -------------------------------------------------------------------------

function [gapData,wasCached] = solveBoomerangBands(P,cfg)
%SOLVEBOOMERANGBANDS Call solveBands and always return a usable gap struct.
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

matPath = [P.datLoc,P.fileBase,'_bds.mat'];
cacheHitExpected = isfile(matPath);

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

% No .full field: solveBands short-circuited on the existing file.
wasCached = true;
if ~cacheHitExpected
    warning('bayesopt_boomerang:noGapData', ...
        'solveBands returned no .full field and no cache at %s', matPath);
    return
end

cached = load(matPath,'ds');
if isfield(cached,'ds') && isfield(cached.ds,'full')
    gapData = cached.ds.full;
    fprintf('  (reused cached band structure)\n');
else
    warning('bayesopt_boomerang:badCache', ...
        'Cached file %s contains no ds.full', matPath);
end
end

% -------------------------------------------------------------------------

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

if isfile(logPath)
    return
end
fid = fopen(logPath,'wt+');
if fid < 0
    warning('bayesopt_boomerang:logOpen','Could not create log %s',logPath);
    return
end
fprintf(fid,['a_nm\tr_nm\tw_nm\tth_nm\tmidGap_GHz\tgapSize_GHz\t' ...
    'gapRatio\tpenalty\tfitness\tobjective\tcached\tevalTime_min\tstatus\r\n']);
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
fprintf(fid,'%d\t%d\t%d\t%d\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%d\t%.3f\t%s\r\n', ...
    x.a, x.r, x.w, x.th, ...
    detail.midGap/1e9, detail.gapSize/1e9, detail.gapRatio, ...
    detail.penalty, detail.fitness, objective, wasCached, evalTime/60, status);
fclose(fid);
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
ds = loadCachedBandData(P);

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
    figG = newSummaryFigure('boomerang geometry',12,5.5);
    figs(end+1) = figG;
    % Default (loose) Padding here, not 'compact': both panels carry two-line
    % axes titles, and compact padding shrinks the strip reserved for the
    % layout title until it prints on top of them.
    tlG = tiledlayout(figG,1,2,'TileSpacing','compact');
    axG1 = nexttile(tlG);
    tryPanel(@() drawBoomerangCell(axG1,P,1,true),'geometry (single cell)');
    axG2 = nexttile(tlG);
    tryPanel(@() drawBoomerangCell(axG2,P,cfg.figNPeriods,false),'geometry (tiled)');
    addFigureHeadline(figG,tlG,sprintf( ...
        'Boomerang unit cell  |  a = %.0f   r = %.0f   w = %.0f   th = %.0f nm', ...
        P.a*1e9,P.r*1e9,P.w*1e9,P.th*1e9),0.09);
    saveFigurePair(figG,[cfg.datLoc,'bayesopt_boomerang_geometry'],cfg);
catch ME
    warning('bayesopt_boomerang:geometryFigure', ...
        'Could not build the geometry figure: %s',ME.message);
end

% --- standalone band structure figure -----------------------------------
try
    figB = newSummaryFigure('boomerang best bands',8,6);
    figs(end+1) = figB;
    % 'Parent' name-value rather than axes(figB): the positional form was only
    % added in R2016a and this is the one place a plain single axes is wanted.
    axB = axes('Parent',figB);
    tryPanel(@() plotBestBands(axB,ds,cfg,detail),'band structure');
    saveFigurePair(figB,[cfg.datLoc,'bayesopt_boomerang_bestbands'],cfg);
catch ME
    warning('bayesopt_boomerang:bandFigure', ...
        'Could not build the band structure figure: %s',ME.message);
end

% --- standalone progress figure -----------------------------------------
% Convergence on top of the four design-space slices, because the two are read
% together: a flat convergence curve plus a variable pinned at its bound means
% the bound is setting the answer, not the physics.
try
    figP = newSummaryFigure('boomerang optimizer progress',13,7);
    figs(end+1) = figP;
    tlP = tiledlayout(figP,2,4,'TileSpacing','compact','Padding','compact');
    axP = nexttile(tlP,[1 4]);
    tryPanel(@() plotConvergenceTrace(axP,results),'convergence');
    drawDesignSpaceRow(tlP,results,cfg);
    addFigureHeadline(figP,tlP,sprintf( ...
        'Optimizer progress  |  %d evaluations  |  best fitness = %.4f', ...
        numel(results.ObjectiveTrace),-results.MinObjective),0.07);
    saveFigurePair(figP,[cfg.datLoc,'bayesopt_boomerang_progress'],cfg);
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

function ds = loadCachedBandData(P)
%LOADCACHEDBANDDATA Reload the _bds.mat solveBands wrote for one design.
%
% Returns [] rather than throwing whenever the data cannot be had, so a missing
% or half-written cache degrades to a figure with an empty band panel instead of
% taking down the whole post-processing step. Reasons the file can legitimately
% be absent: the winning evaluation errored out, datLoc was moved mid-run, or
% the best point came from a resumed study whose data lives under another date
% folder.

ds = [];
matPath = [P.datLoc,P.fileBase,'_bds.mat'];

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

function addFigureHeadline(figH,tl,headlineText,fracHeight)
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

tl.OuterPosition = [0, 0, 1, 1-fracHeight];
annotation(figH,'textbox',[0, 1-fracHeight, 1, fracHeight], ...
    'String',headlineText, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'EdgeColor','none','FitBoxToText','off', ...
    'FontWeight','bold','FontSize',11,'Interpreter','none');
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
% pathNoExt is built by concatenating onto cfg.datLoc, whose trailing separator
% is deliberate - see the note at the cfg.datLoc definition. Do not insert
% another one here.

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
figH = newSummaryFigure('bayesopt boomerang summary',16,10);
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
addFigureHeadline(figH,tl,{ ...
    sprintf('Best boomerang unit cell   |   fitness = %.4f   (objective = %.4f)', ...
        -results.MinObjective,results.MinObjective), ...
    sprintf(['a=%d  r=%d  w=%d  th=%d nm   |   mid-gap = %.3f GHz   ' ...
             'gap ratio = %.2f%%   target = %.2f GHz   |   %d evaluations   |   %s'], ...
        xBest.a,xBest.r,xBest.w,xBest.th, ...
        detail.midGap/1e9,detail.gapRatio*100,cfg.targetFreq/1e9, ...
        numel(results.ObjectiveTrace),detail.status)},0.075);

saveFigurePair(figH,[cfg.datLoc,'bayesopt_boomerang_summary'],cfg);
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
title(ax,{'Mechanical band structure of the best design',ttl2},'FontSize',9);

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
    sprintf('  th      = %6d nm   full slab thickness in z',xBest.th), ...
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
