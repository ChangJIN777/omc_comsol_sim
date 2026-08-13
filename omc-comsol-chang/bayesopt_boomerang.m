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
