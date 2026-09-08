%% BAYESOPT_CROSS
% Bayesian optimization of the diamond CROSS unit cell for a mechanical
% bandgap. The counterpart of cross_optimize_sweep_diamond.m (Nelder-Mead),
% built on the pattern in bayesopt_boomerang.m and driving the same solveBands
% pipeline that test_CrossUnitCell.m uses.
%
% WHY BAYESOPT HERE
% Each evaluation is a multi-minute COMSOL band-structure solve, so the
% EVALUATION BUDGET, not the convergence rate, is what binds. bayesopt fits a
% Gaussian-process surrogate and spends each solve where it is most
% informative. Three further properties of this objective suit it:
%   * a large part of the design box has NO complete gap, so the objective is
%     exactly 0 over a plateau that a simplex cannot descend;
%   * unfabricable candidates can be pruned by XConstraintFcn BEFORE any solve,
%     rather than paid for and then discarded;
%   * an interrupted study can be resumed from its checkpoint.
%
% GEOMETRY -- the parameter meanings are NOT what the comments in
% test_CrossUnitCell.m say. Verified against buildCrossUnitCell.m:
%   a   square cell side, the lattice constant in BOTH x and y   (:38, [a a])
%   h   LENGTH of each cross arm                          (:43 [h w], :48 [w h])
%   w   WIDTH of each cross arm                     (NOT the slab thickness)
%   th  slab THICKNESS along z                                   (:33, :65)
%   r1, r2  fillet radii, held fixed                            (:137, :147)
% The cross is SUBTRACTED ('r_ucell-r1-r2', :53): the solid is a square slab
% with a cross-shaped void.
%
% REQUIREMENTS
%   Statistics and Machine Learning Toolbox  (bayesopt, optimizableVariable)
%   COMSOL with LiveLink for MATLAB, running -- unless cfg.solverBackend is a
%   dry-run backend, which needs neither COMSOL nor a licence for it.
%   Run from omc-comsol-chang/ so solveBands and friends are on the path.
%
% USAGE
%   1. cfg.solverBackend = 'surrogate';  run. The WHOLE loop executes in
%      seconds with no COMSOL, so the plumbing is debugged first.
%   2. cfg.solverBackend = 'comsol';     run the real study.
%   3. To extend a finished study:
%        S = load(cfg.statePath);
%        results = resume(S.results, 'MaxObjectiveEvaluations', 20);
%
% See also CROSS_OPTIMIZE_SWEEP_DIAMOND, ISCROSSFABRICABLE, SURROGATECROSSBANDS.

clear all; clc; close all; %#ok<CLALL>

%% ===================== CONFIGURATION =====================
% --- solver backend -------------------------------------------------------
% A STRING, not a boolean dryRun flag: three tiers do not fit in a boolean, and
% a name self-documents everywhere it is printed, logged or saved.
%   'comsol'    - the real solveBands call. The only real physics.
%   'surrogate' - analytic fake bands (surrogateCrossBands). Debugs the loop.
%   'stub'      - no bands at all; every evaluation returns NaN.
% cfg.solverBackend = 'comsol';
cfg.solverBackend = 'surrogate';

assertKnownBackend(cfg.solverBackend);
cfg.isDryRun = ~strcmp(cfg.solverBackend, 'comsol');

cfg.dryRunDelay     = 0;          % artificial seconds per evaluation
cfg.dryRunFailEvery = 0;          % 0 = off; N > 0 fails every Nth evaluation
cfg.dryRunPrefname  = 'DRYRUN';

% --- design variables -----------------------------------------------------
% INTEGER NANOMETRES, following bayesopt_boomerang.m. Two reasons: the design
% is patterned on a lithography grid, so a continuous variable offers precision
% the process cannot deliver; and an integer grid lets bayesopt recognise a
% repeated design instead of re-solving a point 0.01 nm away. The objective
% converts back to metres.
cfg.bounds.a  = [120, 600];       % nm, square cell side
cfg.bounds.h  = [ 100, 600];       % nm, arm length
cfg.bounds.w  = [ 50, 400];       % nm, arm width
cfg.bounds.th = [ 150, 400];       % nm, slab thickness (ignored if cfg.fixTh)

% --- hold the slab thickness fixed? ---------------------------------------
% true  -> th is NOT searched; it is held at cfg.thFixed and cfg.bounds.th is
%          ignored. The study drops from four dimensions to three.
% false -> th is a fourth optimizable variable over cfg.bounds.th.
%
% Fixing it is usually the right call, for the reason bayesopt_boomerang.m
% gives for doing the same (:270-277): th is set by whatever film the process
% actually delivers, and there is no point optimizing a dimension you cannot
% choose. It also makes the surface markedly easier for a GP to learn on a
% small evaluation budget -- with 40 evaluations, three dimensions is a very
% different problem from four.
%
% Leave it false only if you genuinely can pick the thickness.
cfg.fixTh   = true;
cfg.thFixed = 250e-9;             % m, used only when cfg.fixTh is true

cfg.r1 = 10e-9;                   % fillet radii, held fixed
cfg.r2 = 10e-9;

% --- objective ------------------------------------------------------------
% NOTE ON SCALE. A ~160 nm cell in diamond is a TENS-of-GHz structure: with a
% shear velocity near 12000 m/s the Bragg frequency is v/(2a) ~ 37 GHz. A
% 10 GHz target would discard every gap this cell has.
cfg.targetFreq = 6.5e9;            % Hz, mechanical midgap target
cfg.sigma      = 1e9;             % Hz, Gaussian width of the frequency penalty
cfg.maxFreq    = 20e9;           % ignore gaps above this

% Smallest solid wall / etched gap. 10 nm, NOT the 50 nm used in the nanobeam
% code: this cell is designed an order of magnitude finer.
%
% The binding solid feature is the wall between the voids of ADJACENT cells,
% a - h, which for the nominal 160/140 nm design is 20 nm -- twice this limit,
% so the nominal starts with real margin. isCrossFabricable used to constrain
% the half-distance to the cell boundary instead, which made the nominal look
% like it sat exactly on the limit and halved the reachable box.
cfg.minFeature = 10e-9;

% --- solver settings passed through to solveBands -------------------------
cfg.kpts     = 9;                 % k-points EXCLUDING gamma
cfg.nbands   = 10;
cfg.meshSize = 4;
cfg.max_dof  = 3e6;

% --- bayesopt settings ----------------------------------------------------
cfg.maxEvaluations = 60;          % TOTAL study size. Every one is a COMSOL
                                  % solve under the 'comsol' backend -- time a
                                  % single evaluation before raising this.
cfg.numSeedPoints  = 8;
cfg.acquisition    = 'expected-improvement-plus';
% A COMSOL band-structure solve is repeatable for a given geometry, and the
% integer-nm grid removes the mesh-jitter-on-a-nudged-geometry problem that
% would otherwise argue for false here.
cfg.isDeterministic = true;
cfg.plotFcn = {@plotMinObjective, @plotObjectiveModel};

cfg.runTag = 'bayesopt_cross_run1';

%% ===================== OUTPUT LOCATIONS =====================
currentDate = datestr(now, 'mmddyyyy'); %#ok<TNOW1,DATST>

% PROVENANCE GUARD 1 of 3: synthetic results live in their own dated folder,
% so a real run never even looks in the directory the fabricated files are in.
if cfg.isDryRun
    datedFolder    = [cfg.dryRunPrefname, '_', cfg.solverBackend, '_', currentDate];
    cfg.filePrefix = [cfg.dryRunPrefname, '_'];
else
    datedFolder    = currentDate;
    cfg.filePrefix = '';
end

% PROVENANCE GUARD 2 of 3: a DRYRUN filename prefix on every synthetic
% artifact, so the two could not collide even if the folders were merged.
cfg.rootLoc = [fullfile('.', 'test', cfg.runTag, datedFolder), filesep];
if ~exist(cfg.rootLoc, 'dir')
    mkdir(cfg.rootLoc);
end
cfg.statePath = [cfg.rootLoc, cfg.filePrefix, 'bayesopt_cross_state.mat'];
cfg.itrPath   = [cfg.rootLoc, cfg.filePrefix, 'bayesopt_cross_log.txt'];

if ~exist(cfg.itrPath, 'file')
    itr = fopen(cfg.itrPath, 'wt+');
    fprintf(itr, ['eval\ta_nm\th_nm\tw_nm\tth_nm\tmidGap_Hz\tgapSize_Hz\t', ...
                  'gapFrac\tobjective\tstatus\r\n']);
    fclose(itr);
end

if cfg.isDryRun
    fprintf(['*** DRY RUN: backend = ''%s''. Every number below is ', ...
             'FABRICATED and must\n*** not be used to judge a design. ', ...
             'Output -> %s\n\n'], cfg.solverBackend, cfg.rootLoc);
end

%% ===================== FEASIBILITY SCAN =====================
% What fraction of the box is even fabricable? A mostly-infeasible box is worth
% knowing about BEFORE spending a day of COMSOL, because bayesopt draws its
% seed points from the whole box and most would be rejected.
nScan = 21;
if cfg.fixTh
    [sa, sh, sw] = ndgrid( ...
        linspace(cfg.bounds.a(1), cfg.bounds.a(2), nScan)*1e-9, ...
        linspace(cfg.bounds.h(1), cfg.bounds.h(2), nScan)*1e-9, ...
        linspace(cfg.bounds.w(1), cfg.bounds.w(2), nScan)*1e-9);
    sth  = cfg.thFixed;
    nDim = 3;
else
    [sa, sh, sw, sth] = ndgrid( ...
        linspace(cfg.bounds.a(1),  cfg.bounds.a(2),  nScan)*1e-9, ...
        linspace(cfg.bounds.h(1),  cfg.bounds.h(2),  nScan)*1e-9, ...
        linspace(cfg.bounds.w(1),  cfg.bounds.w(2),  nScan)*1e-9, ...
        linspace(cfg.bounds.th(1), cfg.bounds.th(2), nScan)*1e-9);
    nDim = 4;
end
scanOK  = isCrossFabricable(sa, sh, sw, sth, cfg.minFeature, cfg.r1, cfg.r2);
fracFab = mean(scanOK(:));
fprintf('Fabricable: %.1f%% of the %d^%d design box\n', 100*fracFab, nScan, nDim);
if fracFab == 0
    error('bayesopt_cross:boxInfeasible', ...
        ['No point on a %d^4 scan of cfg.bounds is fabricable, so the study ', ...
         'was NOT started. Widen the bounds or relax cfg.minFeature.'], nScan);
elseif fracFab < 0.10
    warning('bayesopt_cross:boxMostlyInfeasible', ...
        ['Only %.1f%% of the box is fabricable. bayesopt can work here, but ', ...
         'its seed points are drawn from the whole box, so most will be ', ...
         'rejected and the GP may start with very little. Consider ', ...
         'tightening cfg.bounds around the feasible region or raising ', ...
         'cfg.numSeedPoints.'], 100*fracFab);
end

%% ===================== RUN =====================
optVars = [ ...
    optimizableVariable('a',  cfg.bounds.a,  'Type', 'integer'), ...
    optimizableVariable('h',  cfg.bounds.h,  'Type', 'integer'), ...
    optimizableVariable('w',  cfg.bounds.w,  'Type', 'integer')];
if ~cfg.fixTh
    optVars(end+1) = optimizableVariable('th', cfg.bounds.th, 'Type', 'integer');
end
% When cfg.fixTh, th is deliberately absent from optVars -- see cfg.thFixed.
% Every consumer reads the thickness through thOf(), so nothing below has to
% know which of the two cases it is in.

crossObjective([], cfg, true);          % reset the persistent eval counter

objFcn  = @(x) crossObjective(x, cfg);
% XConstraintFcn prunes CANDIDATES, before the objective is called at all, so
% an unfabricable design never costs a solve. That is strictly cheaper than
% evaluating and penalising it.
xConFcn = @(t) isCrossFabricable(t.a*1e-9, t.h*1e-9, t.w*1e-9, thOf(t, cfg), ...
                                 cfg.minFeature, cfg.r1, cfg.r2);
outFcn  = @(res, state) checkpointState(res, state, cfg.statePath);

fprintf('Starting Bayesian optimization: %d evaluations, target %.1f GHz\n', ...
    cfg.maxEvaluations, cfg.targetFreq*1e-9);
fprintf('Backend: %s   Results -> %s\n\n', cfg.solverBackend, cfg.rootLoc);

results = bayesopt(objFcn, optVars, ...
    'MaxObjectiveEvaluations',  cfg.maxEvaluations, ...
    'NumSeedPoints',            cfg.numSeedPoints, ...
    'AcquisitionFunctionName',  cfg.acquisition, ...
    'IsObjectiveDeterministic', cfg.isDeterministic, ...
    'XConstraintFcn',           xConFcn, ...
    'OutputFcn',                outFcn, ...
    'PlotFcn',                  cfg.plotFcn, ...
    'Verbose',                  1);

%% ===================== REPORT =====================
xBest = results.XAtMinObjective;
fprintf('\n==================== BEST DESIGN ====================\n');
fprintf('  a  = %d nm\n  h  = %d nm\n  w  = %d nm\n', ...
    xBest.a, xBest.h, xBest.w);
if cfg.fixTh
    fprintf('  th = %.0f nm  (FIXED, not optimized)\n', cfg.thFixed*1e9);
else
    fprintf('  th = %d nm\n', xBest.th);
end
fprintf('  objective = %.6g   (fitness = %.6g)\n', ...
    results.MinObjective, -results.MinObjective);
fprintf('  inter-cell wall a-h = %.1f nm\n', xBest.a - xBest.h);
if cfg.isDryRun
    fprintf('  *** SYNTHETIC -- backend was ''%s'' ***\n', cfg.solverBackend);
end
fprintf('  checkpoint -> %s\n', cfg.statePath);
fprintf('  log        -> %s\n', cfg.itrPath);

try
    save([cfg.rootLoc, cfg.filePrefix, 'bayesopt_cross_results.mat'], ...
        'results', 'cfg');
catch ME
    warning('bayesopt_cross:saveFailed', ...
        'Could not save the results object: %s', ME.message);
end


%% ========================================================================
function objective = crossObjective(x, cfg, doReset)
%CROSSOBJECTIVE  One bayesopt evaluation. Returns the NEGATED fitness.
%
%   bayesopt passes a one-row TABLE with the optimizableVariable names, in
%   INTEGER NANOMETRES; everything below works in metres.
%
%   crossObjective([], cfg, true) resets the persistent counter and returns.
%   The counter is persistent and the script's own `clear all` cannot reach a
%   persistent inside a function of the file currently executing, so without
%   this a second run in the same session keeps numbering from where the last
%   one stopped.

persistent nEval
if nargin >= 3 && doReset
    nEval     = 0;
    objective = 0;
    return;
end
if isempty(nEval); nEval = 0; end
nEval = nEval + 1;

P = buildCrossP(x, cfg);

% Belt and braces: XConstraintFcn should already have pruned this, but a
% future caller might invoke the objective directly.
if ~isCrossFabricable(P.a, P.h, P.w, P.th, cfg.minFeature, cfg.r1, cfg.r2)
    objective = 0;      % no gap credit; bayesopt minimizes, best is negative
    logRow(cfg.itrPath, nEval, P, NaN, NaN, NaN, objective, 'unfabricable');
    return;
end

P.datLoc = [cfg.rootLoc, cfg.filePrefix, ...
    sprintf('eval_%04d_a%d_h%d_w%d_th%.0f', nEval, x.a, x.h, x.w, P.th*1e9), ...
    filesep];

status = 'ok';
midGap = NaN;  gapSize = NaN;  gapRat = NaN;
try
    bds = solveCrossBackend(P, cfg, nEval);

    % PROVENANCE GUARD 3 of 3: a real run that is somehow handed synthetic data
    % errors out instead of silently using it.
    if ~cfg.isDryRun && isfield(bds, 'isSynthetic') && bds.isSynthetic
        error('bayesopt_cross:syntheticInRealRun', ...
            'A real run received SYNTHETIC band data. Refusing to use it.');
    end

    % solveBands SKIPS the solve when <fileBase>_bds.mat already exists in
    % datLoc, returning a stub with .sym/.asym and NO .full. Each evaluation
    % gets its own folder to avoid that, and this catches it if one slips.
    if ~isfield(bds, 'full')
        error('bayesopt_cross:noFull', ...
            ['solveBands returned no .full field -- it took its cache ', ...
             'branch, so this design was already solved in %s'], P.datLoc);
    end

    keep = bds.full.midGap <= cfg.maxFreq;
    gs   = bds.full.gapSize(keep);
    mg   = bds.full.midGap(keep);

    if isempty(gs)
        status  = 'no_gap';
        gapSize = 0;
        midGap  = cfg.targetFreq;
        gapRat  = 0;
    else
        [gapSize, iBest] = max(gs);
        midGap = mg(iBest);
        if midGap > 0
            gapRat = gapSize / midGap;
        else
            gapRat = 0;
        end
    end
catch ME
    % A failed solve is a legitimate outcome of a design, not a reason to
    % abandon a study that may already be COMSOL-hours deep.
    status = 'failed';
    warning('bayesopt_cross:evalFailed', 'eval %d failed: %s', nEval, ME.message);
    objective = 0;
    logRow(cfg.itrPath, nEval, P, midGap, gapSize, gapRat, objective, status);
    return;
end

% bayesopt MINIMIZES, and calFitnessCross already returns "lower is better"
% (it negates internally), so it is passed straight through.
objective = calFitnessCross(midGap, gapRat, cfg.targetFreq, cfg.sigma);

fprintf(['  eval %3d: a=%d h=%d w=%d th=%.0f%s nm -> midGap=%.2f GHz  ', ...
         'gapFrac=%.4f  obj=%.6g  [%s]%s\n'], ...
    nEval, x.a, x.h, x.w, P.th*1e9, repmat('*', 1, cfg.fixTh), ...
    midGap*1e-9, gapRat, objective, status, ...
    repmat('  [SYNTHETIC]', 1, cfg.isDryRun));

logRow(cfg.itrPath, nEval, P, midGap, gapSize, gapRat, objective, status);
end

%% ========================================================================
function th = thOf(t, cfg)
%THOF  Slab thickness for a candidate, whichever mode the study is in.
%
%   The single place that knows whether th came from the optimizer or from
%   cfg.thFixed. Everything else -- the constraint function, the objective, the
%   P builder, the log -- goes through here, so adding the fixed-th option did
%   not require each of them to branch, and did not risk one of them being
%   missed. That is how a "fixed" dimension quietly starts varying again.
%
%   Returns METRES. t.th, when present, is in integer nanometres.
if cfg.fixTh
    th = cfg.thFixed * ones(height(t), 1);
else
    th = double(t.th) * 1e-9;
end
end

%% ========================================================================
function P = buildCrossP(x, cfg)
%BUILDCROSSP  Design vector (integer nm) -> the P struct solveBands expects.
P.a  = double(x.a) * 1e-9;
P.h  = double(x.h) * 1e-9;
P.w  = double(x.w) * 1e-9;
P.th = thOf(x, cfg);       % from the optimizer, or held at cfg.thFixed
P.r1 = cfg.r1;
P.r2 = cfg.r2;

P.xsect      = 'rect';
P.beamMat    = 'diamond';
P.celltype   = 'cross';
P.unitcell   = 'square';
P.nperiod    = 1;
P.holeatedge = 0;

P.kpts   = cfg.kpts;
P.nbands = cfg.nbands;

% solveBands gates on both of these before any solve and neither has a default
% (solveBands.m:138 / :186). TwoSymPlanes = 0 means a single symmetry plane,
% two sectors; zSymCondition = 1 makes that plane z, matching mbeveny/mbevenz.
P.TwoSymPlanes  = 0;
P.zSymCondition = 1;
P.solveasym     = 1;

P.completeBandGaps = 1;
P.plotgeom    = 0;
P.savedat     = 1;
P.savebndplot = 0;      % off during a search: one figure per evaluation adds up
P.saveplots   = 0;
P.saveMPH     = 0;
P.bandStruct_2D = 1;    % -> runBands_2D; also skips the P.wc title block,
                        % which reads a field nothing ever assigns
                        % (solveBands.m:478, reachable only when this is 0)

P.mbeveny  = 0;
P.mbevenz  = 1;
P.freq     = 0;
P.meshSize = cfg.meshSize;
P.fixed_bc = 0;

P.anisoMat = 1;
P.rxtal    = 45;
P.max_dof  = cfg.max_dof;
end

%% ========================================================================
function bds = solveCrossBackend(P, cfg, nEval)
%SOLVECROSSBACKEND  Backend dispatcher: real solve, surrogate, or stub.
%
%   Validating the name a SECOND time here is deliberate: the check at
%   configuration time catches a typo before any path is derived from it, this
%   one catches a cfg hand-edited or reloaded from an old results file after.
assertKnownBackend(cfg.solverBackend);

if cfg.isDryRun && cfg.dryRunFailEvery > 0 && mod(nEval, cfg.dryRunFailEvery) == 0
    error('bayesopt_cross:injectedFailure', ...
        'injected dry-run failure (dryRunFailEvery = %d)', cfg.dryRunFailEvery);
end
if cfg.isDryRun && cfg.dryRunDelay > 0
    pause(cfg.dryRunDelay);
end

switch cfg.solverBackend
    case 'comsol'
        bds = solveBands(P);
    case 'surrogate'
        bds = surrogateCrossBands(P, cfg);
    case 'stub'
        bds.sym  = struct('F', [], 'k_norm', []);
        bds.asym = struct('F', [], 'k_norm', []);
        bds.full = struct('F', [], 'midGap', [], 'gapSize', []);
        bds.isSynthetic = true;
end
end

%% ========================================================================
function f = calFitnessCross(midGap, gapRat, targetFreq, sigma)
%CALFITNESSCROSS  Objective for the mechanical bandgap. LOWER is better.
%
%   f = -( gapRat * exp(-((targetFreq - midGap)/sigma)^2) )
%
%   Negated so it can be handed straight to a minimizer, matching the sign
%   convention in boomerang_optimize_sweep_diamond.m. The frequency term is
%   Gaussian, so sigma is the 1/e half-width: landing sigma away from the
%   target costs a factor 1/e, 2*sigma costs 1/e^4. A gapless design scores
%   exactly 0, which is the worst attainable value and the plateau the GP has
%   to model.
freqPenalty = exp(-((targetFreq - midGap)/sigma)^2);
f = -(gapRat * freqPenalty);
end

%% ========================================================================
function stop = checkpointState(results, state, statePath)
%CHECKPOINTSTATE  Save the BayesianOptimization object after every iteration.
%
%   With multi-minute COMSOL evaluations, losing a partly finished study to a
%   crashed LiveLink connection is expensive. Reload and extend with:
%     S = load(statePath);
%     results = resume(S.results, 'MaxObjectiveEvaluations', 20);
stop = false;
switch state
    case {'iteration', 'done'}
        try
            save(statePath, 'results');
        catch ME
            warning('bayesopt_cross:checkpoint', ...
                'Could not checkpoint to %s: %s', statePath, ME.message);
        end
end
end

%% ========================================================================
function assertKnownBackend(backend)
%ASSERTKNOWNBACKEND  Reject an unrecognised cfg.solverBackend, loudly.
%
%   An error rather than a warning because there is no safe default. Falling
%   back to 'comsol' would spend hours of solve time on a typo; falling back to
%   'surrogate' would hand back a folder of fabricated results that look
%   exactly like a study. Stopping is the only option that cannot mislead.
known = {'comsol', 'surrogate', 'stub'};
if ~(ischar(backend) || (isstring(backend) && isscalar(backend))) ...
        || ~any(strcmp(char(backend), known))
    if isnumeric(backend) || islogical(backend)
        got = mat2str(backend);
    elseif ischar(backend) || isstring(backend)
        got = ['''', char(backend), ''''];
    else
        got = ['<', class(backend), '>'];
    end
    error('bayesopt_cross:unknownBackend', ...
        ['Unknown cfg.solverBackend. Expected one of: ''comsol'', ', ...
         '''surrogate'', ''stub''.\nGot: %s\n', ...
         '  ''comsol''    - the real COMSOL solve (the only real physics)\n', ...
         '  ''surrogate'' - analytic fake bands, for debugging the loop\n', ...
         '  ''stub''      - no bands at all, every evaluation returns NaN'], got);
end
end

%% ========================================================================
function logRow(itrPath, nEval, P, midGap, gapSize, gapRat, objective, status)
%LOGROW  Append one evaluation. Opened and closed per row so the log survives a
%   mid-run crash -- the point of keeping it separate from the checkpoint.
try
    itr = fopen(itrPath, 'at+');
    if itr < 0; return; end
    fprintf(itr, '%d\t%.1f\t%.1f\t%.1f\t%.1f\t%.6e\t%.6e\t%.6e\t%.6g\t%s\r\n', ...
        nEval, P.a*1e9, P.h*1e9, P.w*1e9, P.th*1e9, ...
        midGap, gapSize, gapRat, objective, status);
    fclose(itr);
catch
    % Logging must never abort a study that is already COMSOL-hours deep.
end
end
