function results = bayesoptUnitCell(cfg, adapter)
%BAYESOPTUNITCELL  Shared bayesopt driver for 1-D unit-cell band-gap studies.
%
%   results = bayesoptUnitCell(cfg, adapter)
%   cfg     = bayesoptUnitCell('defaults')     % driver-level defaults only
%
%   The celltype-independent part of bayesopt_cross.m, driving solveBands
%   (Pipeline A). Everything specific to a cell comes in through ADAPTER, a
%   struct returned by unitCellAdapter_hole / unitCellAdapter_rib:
%     name        'hole' | 'rib'
%     varNames    ordered cellstr of the design quantities
%     show        per variable: struct('scale',SI->display,'label',..,'unit',..)
%     defaultCfg  @()        -> complete cfg
%     prepareCfg  @(cfg)     -> cfg, after celltype consistency checks
%     derive      @(D,cfg)   -> D with dependent quantities filled in
%     fabOK       @(D,cfg)   -> logical, same size as the broadcast design
%     buildP      @(D,cfg)   -> the P struct solveBands expects (no datLoc)
%     surrogate   @(P,cfg)   -> [fCentre_Hz, gapFrac]   (NOT physics)
%
%   DESIGN VARIABLES. cfg.vars.<name> = struct('search',tf,'bounds',[lo hi],
%   'fixed',valueSI,'step',SIperGridUnit). A searched variable is an INTEGER
%   optimizableVariable over bounds, in units of step (1e-9 m for lengths, so
%   integer nm). A variable that is not searched is held at .fixed (SI).
%   designOf() is the single place that knows which case applies; it is the
%   generalisation of thOf() in bayesopt_cross.m.
%
%   BACKENDS. cfg.solverBackend = 'comsol' | 'surrogate' | 'stub'.
%     'comsol'    the real solveBands call. The only real physics.
%     'surrogate' synthetic bands from adapter.surrogate, passed through the
%                 real findGaps. Debugs the loop.
%     'stub'      no bands at all; every evaluation takes the no_gap path.
%
%   RESUME.  S = load(cfg.statePath);
%            results = resume(S.results, 'MaxObjectiveEvaluations', 20);
%   (20 ADDITIONAL evaluations.) NOT ACROSS A CONSTRAINT CHANGE: XConstraintFcn
%   is baked into the saved object. After changing a fabrication limit, a
%   bound, or which variables are searched, start a fresh study and, if
%   wanted, set cfg.seedFromState to the old state file: its evaluations are
%   re-filtered through the NEW rules and passed as InitialX. They count
%   against cfg.maxEvaluations. Only seed from a study with the same
%   objective settings (targetFreq, sigma, gapField, gapSelect, symmetry).
%
%   Requires the Statistics and Machine Learning Toolbox. Run from
%   omc-comsol-chang/. COMSOL LiveLink is needed only for 'comsol'.
%
%   See also BAYESOPT_HOLE, BAYESOPT_RIB, BAYESOPT_CROSS, SOLVEBANDS.

if nargin == 1 && (ischar(cfg) || isstring(cfg)) && strcmp(cfg, 'defaults')
    results = driverDefaults();
    return;
end
narginchk(2, 2);

%% ===================== CONFIGURATION =====================
def = driverDefaults();
fn  = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(cfg, fn{k})
        cfg.(fn{k}) = def.(fn{k});
    end
end
req     = {'vars', 'targetFreq', 'sigma', 'kpts', 'nbands', 'meshSize', 'max_dof', 'runTag'};
missing = req(~isfield(cfg, req));
if ~isempty(missing)
    error('bayesoptUnitCell:missingCfg', ...
        'cfg is missing: %s. Start from adapter.defaultCfg().', strjoin(missing, ', '));
end

assertKnownBackend(cfg.solverBackend);
cfg.solverBackend = char(cfg.solverBackend);
cfg.isDryRun = ~strcmp(cfg.solverBackend, 'comsol');
cfg.cellName = adapter.name;
if ~any(strcmp(cfg.gapSelect, {'bestScore', 'closest', 'widest'}))
    error('bayesoptUnitCell:badGapSelect', ...
        'cfg.gapSelect must be ''bestScore'', ''closest'' or ''widest''.');
end
cfg = adapter.prepareCfg(cfg);          % celltype consistency; may error

names      = adapter.varNames;
isSearched = false(1, numel(names));
for k = 1:numel(names)
    n = names{k};
    if ~isfield(cfg.vars, n)
        error('bayesoptUnitCell:missingVar', 'cfg.vars.%s is missing.', n);
    end
    v = cfg.vars.(n);
    isSearched(k) = logical(v.search);
    if isSearched(k)
        b = v.bounds;
        if numel(b) ~= 2 || any(b ~= round(b)) || b(1) >= b(2)
            error('bayesoptUnitCell:badBounds', ...
                ['cfg.vars.%s.bounds must be [lo hi] with integer lo < hi, in ', ...
                 'units of cfg.vars.%s.step.'], n, n);
        end
    end
end
cfg.searchedNames = names(isSearched);
nS = numel(cfg.searchedNames);
if nS == 0
    error('bayesoptUnitCell:nothingSearched', ...
        'No variable has cfg.vars.<name>.search = true.');
end

%% ===================== OUTPUT LOCATIONS =====================
% TIME-stamped, not date-stamped: a same-day rerun must not land on the
% first run's cached <fileBase>_bds.mat files (optimize_rib_unitCell.m:83-87).
stamp = char(datetime('now', 'Format', 'MMddyyyy_HHmmss'));

% PROVENANCE GUARD 1 of 3: synthetic results live in their own folder.
% PROVENANCE GUARD 2 of 3: and carry a DRYRUN filename prefix.
if cfg.isDryRun
    runFolder      = [cfg.dryRunPrefname, '_', cfg.solverBackend, '_', stamp];
    cfg.filePrefix = [cfg.dryRunPrefname, '_'];
else
    runFolder      = stamp;
    cfg.filePrefix = '';
end
if isempty(cfg.rootLoc)
    cfg.rootLoc = [fullfile('.', 'test', cfg.runTag, runFolder), filesep];
elseif ~strcmp(cfg.rootLoc(end), filesep)
    cfg.rootLoc = [cfg.rootLoc, filesep];       % caller-supplied (tests)
end
if ~exist(cfg.rootLoc, 'dir')
    mkdir(cfg.rootLoc);
end
baseName      = [cfg.filePrefix, 'bayesopt_', adapter.name];
cfg.statePath = [cfg.rootLoc, baseName, '_state.mat'];
cfg.itrPath   = [cfg.rootLoc, baseName, '_log.txt'];

if ~exist(cfg.itrPath, 'file')
    labels = cellfun(@(n) adapter.show.(n).label, names, 'UniformOutput', false);
    itr = fopen(cfg.itrPath, 'w');
    fprintf(itr, 'eval\t%s\tmidGap_Hz\tgapSize_Hz\tgapFrac\tobjective\tstatus\n', ...
        strjoin(labels, sprintf('\t')));
    fclose(itr);
end

if cfg.isDryRun
    fprintf(['*** DRY RUN: backend = ''%s''. Every number below is FABRICATED ', ...
             'and must\n*** not be used to judge a design. Output -> %s\n\n'], ...
        cfg.solverBackend, cfg.rootLoc);
end

%% ===================== FEASIBILITY SCAN =====================
nScan = max(2, min(cfg.nScan, floor(cfg.maxScanPoints^(1/nS))));
axesC = cell(1, nS);
for k = 1:nS
    b = cfg.vars.(cfg.searchedNames{k}).bounds;
    axesC{k} = unique(round(linspace(b(1), b(2), nScan)));
end
gridC = cell(1, nS);
[gridC{:}] = ndgrid(axesC{:});
G = cell2struct(gridC(:), cfg.searchedNames(:), 1);
scanOK  = fabOfX(G, cfg, adapter);
fracFab = mean(scanOK(:));
fprintf('Fabricable: %.1f%% of a %d-point scan of the %d-D design box\n', ...
    100*fracFab, numel(scanOK), nS);
if fracFab == 0
    error('bayesoptUnitCell:boxInfeasible', ...
        ['No scanned point of the design box is fabricable, so the study was ', ...
         'NOT started. Widen cfg.vars.*.bounds or relax the fabrication limits.']);
elseif fracFab < 0.10
    warning('bayesoptUnitCell:boxMostlyInfeasible', ...
        ['Only %.1f%% of the box is fabricable. Consider tightening the bounds ', ...
         'around the feasible region or raising cfg.numSeedPoints.'], 100*fracFab);
end

%% ===================== RUN =====================
c = cell(1, nS);
for k = 1:nS
    n = cfg.searchedNames{k};
    c{k} = optimizableVariable(n, cfg.vars.(n).bounds, 'Type', 'integer');
end
optVars = [c{:}];

unitCellObjective([], cfg, adapter, true);      % reset the persistent counter

objFcn  = @(x) unitCellObjective(x, cfg, adapter);
% XConstraintFcn receives an N-row TABLE of candidates (thousands of rows per
% iteration) and must return an N-by-1 logical. fabOfX serves it, the scan
% above and the check inside the objective.
xConFcn = @(T) fabOfX(T, cfg, adapter);
outFcn  = @(res, state) checkpointState(res, state, cfg);

boArgs = { ...
    'MaxObjectiveEvaluations',  cfg.maxEvaluations, ...
    'NumSeedPoints',            cfg.numSeedPoints, ...
    'AcquisitionFunctionName',  cfg.acquisition, ...
    'IsObjectiveDeterministic', cfg.isDeterministic, ...
    'XConstraintFcn',           xConFcn, ...
    'OutputFcn',                outFcn, ...
    'PlotFcn',                  cfg.plotFcn, ...
    'Verbose',                  cfg.verbose};
if ~isempty(cfg.seedFromState)
    [X0, obj0] = seedFromState(cfg, adapter);
    if ~isempty(X0)
        boArgs = [boArgs, {'InitialX', X0, 'InitialObjective', obj0}];
    end
end

fprintf('Starting Bayesian optimization of the %s cell: %d evaluations, target %.2f GHz\n', ...
    adapter.name, cfg.maxEvaluations, cfg.targetFreq*1e-9);
fprintf('Searched: %s   Backend: %s   Results -> %s\n\n', ...
    strjoin(cfg.searchedNames, ', '), cfg.solverBackend, cfg.rootLoc);

results = bayesopt(objFcn, optVars, boArgs{:});

%% ===================== REPORT =====================
xBest = results.XAtMinObjective;
fprintf('\n==================== BEST DESIGN (%s) ====================\n', adapter.name);
if isempty(xBest)
    warning('bayesoptUnitCell:noBest', 'bayesopt returned no best point.');
else
    Dbest = designOf(xBest, cfg, adapter);
    for k = 1:numel(names)
        n = names{k};
        if cfg.vars.(n).search
            tag = '';
        elseif isequal(Dbest.(n), cfg.vars.(n).fixed)
            tag = '  (FIXED, not optimized)';
        else
            tag = '  (derived)';
        end
        fprintf('  %-3s = %.4g %s%s\n', n, Dbest.(n)*adapter.show.(n).scale, ...
            adapter.show.(n).unit, tag);
    end
    fprintf('  objective = %.6g   (fitness = %.6g)\n', ...
        results.MinObjective, -results.MinObjective);
    if results.MinObjective >= 0
        warning('bayesoptUnitCell:noUsableGap', ...
            'No evaluated design had a usable gap (best objective is 0). See %s', ...
            cfg.itrPath);
    end
end
if cfg.isDryRun
    fprintf('  *** SYNTHETIC -- backend was ''%s'' ***\n', cfg.solverBackend);
end
fprintf('  checkpoint -> %s\n  log        -> %s\n', cfg.statePath, cfg.itrPath);

try
    save([cfg.rootLoc, baseName, '_results.mat'], 'results', 'cfg');
catch ME
    warning('bayesoptUnitCell:saveFailed', ...
        'Could not save the results object: %s', ME.message);
end
end

%% ========================================================================
function cfg = driverDefaults()
%DRIVERDEFAULTS  Knobs that do not depend on the celltype.
cfg.solverBackend   = 'surrogate';   % safe default: no COMSOL time on a typo-free first run
cfg.dryRunDelay     = 0;             % artificial seconds per evaluation
cfg.dryRunFailEvery = 0;             % N > 0: designs whose grid values sum to a
                                     % multiple of N fail (keyed on the DESIGN,
                                     % so a repeat fails again: deterministic)
cfg.dryRunPrefname  = 'DRYRUN';
cfg.gapField        = 'sym';         % which solveBands sector the gap is read from
cfg.gapSelect       = 'bestScore';   % 'bestScore' | 'closest' | 'widest'
cfg.maxFreq         = Inf;           % ignore gaps above this (Hz)
cfg.maxEvaluations  = 60;            % TOTAL study size (includes seeded points)
cfg.numSeedPoints   = 12;
cfg.acquisition     = 'expected-improvement-plus';
cfg.isDeterministic = true;          % integer grid + repeatable solve
cfg.plotFcn         = {@plotMinObjective, @plotObjectiveModel};   % [] = none
cfg.verbose         = 1;
cfg.nScan           = 21;            % scan points per searched variable (capped)
cfg.maxScanPoints   = 2e6;
cfg.rootLoc         = '';            % '' = ./test/<runTag>/<time stamp>/
cfg.seedFromState   = '';            % path to an earlier *_state.mat, or ''
end

%% ========================================================================
function objective = unitCellObjective(x, cfg, adapter, doReset)
%UNITCELLOBJECTIVE  One bayesopt evaluation. Returns the NEGATED fitness.
%
%   x is a one-row table of the searched variables, in integer grid units.
%   unitCellObjective([], cfg, adapter, true) resets the counter.
%   The counter is tied to the log file: after a MATLAB restart, or when a
%   different study is resumed in the same session, it continues from the
%   number of rows already in THAT study's log.
persistent nEval logPath
if nargin >= 4 && doReset
    nEval     = 0;
    logPath   = cfg.itrPath;
    objective = 0;
    return;
end
if isempty(nEval) || ~strcmp(logPath, cfg.itrPath)
    nEval   = countLogRows(cfg.itrPath);
    logPath = cfg.itrPath;
end
nEval = nEval + 1;

D = designOf(x, cfg, adapter);

% Belt and braces: XConstraintFcn should already have pruned this.
if ~all(adapter.fabOK(D, cfg))
    objective = 0;      % no gap credit; bayesopt minimizes, best is negative
    logRow(cfg, adapter, nEval, D, NaN, NaN, NaN, objective, 'unfabricable');
    return;
end

P = adapter.buildP(D, cfg);
% Unique folder per evaluation: solveBands skips the solve when
% <datLoc><fileBase>_bds.mat exists (solveBands.m:172). fullfile/filesep only.
P.datLoc = [cfg.rootLoc, cfg.filePrefix, sprintf('eval_%04d', nEval), ...
    designTag(x, cfg), filesep];

midGap = NaN;  gapSize = NaN;  gapRat = NaN;
try
    ds = solveBackend(P, cfg, adapter, x);

    % PROVENANCE GUARD 3 of 3.
    if ~cfg.isDryRun && isfield(ds, 'isSynthetic') && ds.isSynthetic
        error('bayesoptUnitCell:syntheticInRealRun', ...
            'A real run received SYNTHETIC band data. Refusing to use it.');
    end
    % solveBands' cache branch returns sym = [] / asym = [] and NO .full
    % (solveBands.m:608-612); the solve branch always sets .full (:281).
    if ~isfield(ds, 'full')
        error('bayesoptUnitCell:noFull', ...
            ['solveBands returned no .full field -- it took its cache ', ...
             'branch, so this design was already solved in %s'], P.datLoc);
    end
    if ~isfield(ds, cfg.gapField) || ~isstruct(ds.(cfg.gapField)) ...
            || ~isfield(ds.(cfg.gapField), 'midGap')
        error('bayesoptUnitCell:noSector', ...
            'solveBands returned no usable ds.%s (check cfg.gapField).', cfg.gapField);
    end

    % THE ONE SEAM for a different gap metric (e.g. a tracked band later).
    [midGap, gapSize, gapRat, status] = pickGap(ds.(cfg.gapField), cfg);
catch ME
    % A failed solve is a legitimate outcome of a design, not a reason to
    % abandon a study that may already be COMSOL-hours deep.
    status = 'failed';
    warning('bayesoptUnitCell:evalFailed', 'eval %d failed: %s', nEval, ME.message);
    objective = 0;
    logRow(cfg, adapter, nEval, D, midGap, gapSize, gapRat, objective, status);
    return;
end

objective = calFitness(midGap, gapRat, cfg.targetFreq, cfg.sigma);

fprintf('  eval %3d: %s -> midGap=%.3f GHz  gapFrac=%.4f  obj=%.6g  [%s]%s\n', ...
    nEval, describeDesign(D, adapter), midGap*1e-9, gapRat, objective, status, ...
    repmat('  [SYNTHETIC]', 1, cfg.isDryRun));

logRow(cfg, adapter, nEval, D, midGap, gapSize, gapRat, objective, status);
end

%% ========================================================================
function D = designOf(X, cfg, adapter)
%DESIGNOF  Grid units -> a struct of SI design quantities.
%
%   X is a table (bayesopt: one row in the objective, N rows in
%   XConstraintFcn) or a struct of equally sized arrays (the scan). A searched
%   variable becomes an array of the size it arrived in; a fixed one stays a
%   scalar and broadcasts. double() is deliberate: were a column ever a true
%   integer class, int32(500)*1e-9 would be 0 and every candidate infeasible.
names = adapter.varNames;
D = struct();
for k = 1:numel(names)
    n = names{k};
    v = cfg.vars.(n);
    if v.search
        D.(n) = double(X.(n)) * v.step;
    else
        D.(n) = v.fixed;
    end
end
D = adapter.derive(D, cfg);
end

%% ========================================================================
function ok = fabOfX(X, cfg, adapter)
%FABOFX  Fabricability of grid-unit candidates: designOf + the adapter rule.
ok = adapter.fabOK(designOf(X, cfg, adapter), cfg);
end

%% ========================================================================
function [midGap, gapSize, gapRat, status] = pickGap(sec, cfg)
%PICKGAP  Choose one gap of a sector for scoring.
%   Zero-size entries are dropped: findGaps returns them where bands touch.
status = 'ok';
mg   = sec.midGap(:);
gs   = sec.gapSize(:);
keep = gs > 0 & mg > 0 & mg <= cfg.maxFreq;
mg   = mg(keep);
gs   = gs(keep);
if isempty(gs)
    status  = 'no_gap';
    gapSize = 0;
    midGap  = cfg.targetFreq;
    gapRat  = 0;
    return;
end
ratio = gs ./ mg;
switch cfg.gapSelect
    case 'bestScore'    % the gap that maximizes the fitness itself
        [~, k] = max(ratio .* exp(-((cfg.targetFreq - mg)/cfg.sigma).^2));
    case 'closest'      % rule of optimize_hole/rib_unitCell (closestGap)
        [~, k] = min(abs(mg - cfg.targetFreq));
    case 'widest'       % rule of bayesopt_cross.m:361
        [~, k] = max(gs);
end
midGap  = mg(k);
gapSize = gs(k);
gapRat  = ratio(k);
end

%% ========================================================================
function f = calFitness(midGap, gapRat, targetFreq, sigma)
%CALFITNESS  f = -( gapRat * exp(-((targetFreq - midGap)/sigma)^2) ). LOWER is
%   better; a gapless design scores exactly 0 (as calFitnessCross).
f = -(gapRat * exp(-((targetFreq - midGap)/sigma)^2));
f = f + 0;          % -0 -> +0, so a gapless row logs as 0, not -0
end

%% ========================================================================
function ds = solveBackend(P, cfg, adapter, x)
%SOLVEBACKEND  Backend dispatcher: real solve, surrogate, or stub.
assertKnownBackend(cfg.solverBackend);      % second check, on purpose

if cfg.isDryRun && cfg.dryRunFailEvery > 0 ...
        && mod(designKey(x, cfg), cfg.dryRunFailEvery) == 0
    error('bayesoptUnitCell:injectedFailure', ...
        'injected dry-run failure (dryRunFailEvery = %d)', cfg.dryRunFailEvery);
end
if cfg.isDryRun && cfg.dryRunDelay > 0
    pause(cfg.dryRunDelay);
end

switch cfg.solverBackend
    case 'comsol'
        % Runs whether solveBands returns or throws.
        tidy = onCleanup(@() removeStrayBackslashDir(P.datLoc));
        ds = solveBands(P);
    case 'surrogate'
        [fCentre, gapFrac] = adapter.surrogate(P, cfg);
        ds = synthBands(fCentre, gapFrac, P, cfg);
        ds.surrogateTag = sprintf('SYNTHETIC band data for the %s cell - NOT physics.', ...
            adapter.name);
    case 'stub'
        empty   = struct('F', [], 'k_norm', [], 'midGap', [], 'gapSize', []);
        ds.sym  = empty;
        ds.full = empty;
        ds.isSynthetic = true;
end
end

%% ========================================================================
function ds = synthBands(fCentre, gapFrac, P, cfg)
%SYNTHBANDS  Synthetic 1-D band structure with a gap of relative width
%   gapFrac at fCentre. NOT PHYSICS. Shaped like a single-plane solveBands
%   result (solveBands.m:263-281): .sym, .asym only if P.solveasym, and .full.
%   midGap/gapSize come from the REAL findGaps, never written by hand.
k_norm = linspace(0, 1, cfg.kpts + 1).';        % Gamma-X, runBands.m:109-113
ds.sym = makeSector(k_norm, cfg.nbands, fCentre, gapFrac, 0);
[ds.sym.midGap, ds.sym.gapSize] = findGaps(ds.sym);
full.F = ds.sym.F;
if isfield(P, 'solveasym') && P.solveasym
    ds.asym = makeSector(k_norm, cfg.nbands, fCentre, gapFrac - 0.02, 0.02);
    [ds.asym.midGap, ds.asym.gapSize] = findGaps(ds.asym);
    full.F = [ds.sym.F, ds.asym.F];
end
[full.midGap, full.gapSize] = findGaps(full);
ds.full = full;
ds.isSynthetic = true;
end

%% ========================================================================
function sec = makeSector(k_norm, nb, fCentre, gapFrac, offset)
%MAKESECTOR  One synthetic sector (after surrogateCrossBands/makeSector).
%   GAPLESS CASE: the two band groups must OVERLAP. If they merely touch
%   (fLow == fHigh), findGaps reports a spurious gap between the two highest
%   LOWER bands, because its test is a strict inequality (findGaps.m:31). That
%   matters here because the objective reads ds.sym, not ds.full.
fC = fCentre * (1 + offset);
if gapFrac > 0
    fLow  = fC * (1 - gapFrac/2);
    fHigh = fC * (1 + gapFrac/2);
else
    fLow  = fC;
    fHigh = 0.97 * fC;
end
nLow   = max(1, floor(nb/2));
nUp    = nb - nLow;
ramp   = (1 - cos(pi*k_norm))/2;    % 0 at Gamma, 1 at X
span   = 1.6;                       % consecutive bands overlap by 0.6 spacings
spread = 0.8;
F = zeros(numel(k_norm), nb);
for m = 1:nLow
    F(:, m) = fLow * ((m-1) + span*ramp) / (nLow - 1 + span);
end
for j = 1:nUp
    F(:, nLow + j) = fHigh * (1 + spread * ((j-1) + span*ramp) / (nUp - 1 + span));
end
sec.F      = F;
sec.k_norm = k_norm;
end

%% ========================================================================
function [X0, obj0] = seedFromState(cfg, adapter)
%SEEDFROMSTATE  Earlier evaluations as InitialX, re-filtered by the CURRENT
%   box and fabricability rules (idea from nanobeam/seedBayesFromStore.m).
S = load(cfg.seedFromState);
if ~isfield(S, 'results')
    error('bayesoptUnitCell:badSeedFile', '%s holds no ''results''.', cfg.seedFromState);
end
if ~cfg.isDryRun && isfield(S, 'cfg') && isfield(S.cfg, 'isDryRun') && S.cfg.isDryRun
    error('bayesoptUnitCell:syntheticSeed', ...
        'Refusing to seed a real study from a DRY-RUN state file (%s).', cfg.seedFromState);
end
X0   = S.results.XTrace;
obj0 = S.results.ObjectiveTrace;
if ~isequal(sort(X0.Properties.VariableNames), sort(cfg.searchedNames))
    error('bayesoptUnitCell:seedVarMismatch', ...
        'The seed study searched {%s}; this one searches {%s}.', ...
        strjoin(X0.Properties.VariableNames, ', '), strjoin(cfg.searchedNames, ', '));
end
X0   = X0(:, cfg.searchedNames);
keep = isfinite(obj0) & fabOfX(X0, cfg, adapter);
for k = 1:numel(cfg.searchedNames)
    n = cfg.searchedNames{k};
    b = cfg.vars.(n).bounds;
    keep = keep & double(X0.(n)) >= b(1) & double(X0.(n)) <= b(2);
end
X0   = X0(keep, :);
obj0 = obj0(keep);
fprintf('Seeding from %s: %d of %d earlier evaluations kept.\n', ...
    cfg.seedFromState, nnz(keep), numel(keep));
end

%% ========================================================================
function stop = checkpointState(results, state, cfg)
%CHECKPOINTSTATE  Save the BayesianOptimization object after every iteration.
stop = false;
switch state
    case {'iteration', 'done'}
        try
            save(cfg.statePath, 'results', 'cfg');
        catch ME
            warning('bayesoptUnitCell:checkpoint', ...
                'Could not checkpoint to %s: %s', cfg.statePath, ME.message);
        end
end
end

%% ========================================================================
function assertKnownBackend(backend)
%ASSERTKNOWNBACKEND  Reject an unrecognised cfg.solverBackend. No safe default.
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
    error('bayesoptUnitCell:unknownBackend', ...
        ['Unknown cfg.solverBackend. Expected ''comsol'', ''surrogate'' or ', ...
         '''stub''.\nGot: %s'], got);
end
end

%% ========================================================================
function removeStrayBackslashDir(datLoc)
%REMOVESTRAYBACKSLASHDIR  solveBands appends a literal '\' to datLoc
%   (solveBands.m:4-8). On Windows that is the separator; elsewhere it makes an
%   empty sub-folder named '\'. Remove it (as optimize_rib_unitCell.m:417-431).
if ispc, return; end
stray = [datLoc, '\'];
if exist(stray, 'dir')
    try
        rmdir(stray);       % fails, harmlessly, if it is not empty
    catch
    end
end
end

%% ========================================================================
function tag = designTag(x, cfg)
%DESIGNTAG  '_t100_s80_L400': searched variables in integer grid units.
tag = '';
for k = 1:numel(cfg.searchedNames)
    n = cfg.searchedNames{k};
    tag = [tag, sprintf('_%s%d', n, double(x.(n)))]; %#ok<AGROW>
end
end

%% ========================================================================
function key = designKey(x, cfg)
%DESIGNKEY  Sum of the integer grid values; keys the injected dry-run failure.
key = 0;
for k = 1:numel(cfg.searchedNames)
    key = key + double(x.(cfg.searchedNames{k}));
end
end

%% ========================================================================
function s = describeDesign(D, adapter)
names = adapter.varNames;
parts = cell(1, numel(names));
for k = 1:numel(names)
    n = names{k};
    parts{k} = sprintf('%s=%.4g%s', n, D.(n)*adapter.show.(n).scale, adapter.show.(n).unit);
end
s = strjoin(parts, ' ');
end

%% ========================================================================
function n = countLogRows(itrPath)
%COUNTLOGROWS  Data rows already in the log (header excluded).
n = 0;
if exist(itrPath, 'file') ~= 2, return; end
try
    n = max(0, numel(strfind(fileread(itrPath), newline)) - 1);
catch
end
end

%% ========================================================================
function logRow(cfg, adapter, nEval, D, midGap, gapSize, gapRat, objective, status)
%LOGROW  Append one evaluation; opened and closed per row so the log survives
%   a crash. Every design quantity is logged, searched, fixed or derived.
try
    names = adapter.varNames;
    vals  = cell(1, numel(names));
    for k = 1:numel(names)
        vals{k} = sprintf('%.4f', D.(names{k})*adapter.show.(names{k}).scale);
    end
    itr = fopen(cfg.itrPath, 'a');
    if itr < 0; return; end
    fprintf(itr, '%d\t%s\t%.6e\t%.6e\t%.6e\t%.6g\t%s\n', nEval, ...
        strjoin(vals, sprintf('\t')), midGap, gapSize, gapRat, objective, status);
    fclose(itr);
catch
    % Logging must never abort a study that is already COMSOL-hours deep.
end
end
