%% CROSS_OPTIMIZE_SWEEP_DIAMOND
% Nelder-Mead optimization of the diamond CROSS unit cell for a mechanical
% bandgap, modelled on boomerang_optimize_sweep_diamond.m and driving the same
% solveBands pipeline that test_CrossUnitCell.m uses.
%
% READ THIS FIRST -- the parameter meanings are NOT what the comments in
% test_CrossUnitCell.m say. Verified against buildCrossUnitCell.m:
%
%   P.a   square unit-cell side. The lattice constant in BOTH x and y
%         (ucellPlane size is [a a]), not just x.
%   P.h   LENGTH of each cross arm            (rec_1 size [h w], rec_2 [w h])
%   P.w   WIDTH of each cross arm             (NOT the slab thickness)
%   P.th  slab THICKNESS along z              (workplane at quickz = -th/2)
%   P.r1  fillet radius 1                     (addFillet: fil1.set('radius',r1))
%   P.r2  fillet radius 2                     (addFillet: fil2.set('radius',r2))
%
% The compose formula is 'r_ucell-r1-r2', so the cross is SUBTRACTED: the
% solid is a square slab with a cross-shaped VOID through it. That is why the
% fabrication limits below are written on the remaining solid ligaments rather
% than on the arms themselves.
%
% r1/r2 are process-driven rounding, not design freedom, so they are held
% fixed and the search runs over [a, h, w, th].
%
% Requires: COMSOL with LiveLink for MATLAB running; run from omc-comsol-chang/.
%
% Usage:
%   1. Set OPT.solverBackend = 'surrogate' and run: the WHOLE loop executes in
%      seconds with no COMSOL, so the plumbing can be debugged first.
%   2. Set it back to 'comsol', run once for a single real evaluation.
%   3. Uncomment the fminsearch block in the driver for the real search.

clear all; clc; close all; %#ok<CLALL>

%% ===================== INITIAL DESIGN =====================
% Starting point from test_CrossUnitCell.m.
a0  = 160e-9;      % square cell side
h0  = 140e-9;      % cross arm length
w0  =  50e-9;      % cross arm width
th0 = 120e-9;      % slab thickness (z)

% --- hold the slab thickness fixed? ---------------------------------------
% true  -> th is NOT searched; it is held at th0 and the design vector drops
%          to three elements [a, h, w].
% false -> th is the fourth element of the design vector.
%
% Fixing it is usually the right call, for the reason bayesopt_boomerang.m
% gives for doing the same (:270-277): th is set by whatever film the process
% actually delivers, and there is no point optimizing a dimension you cannot
% choose. It also drops the simplex from four dimensions to three, which for
% Nelder-Mead means a 4-vertex simplex instead of 5 and noticeably fewer
% evaluations to contract.
%
% Set the SAME way in bayesopt_cross.m (cfg.fixTh) if you want the two
% optimizers to be searching the same space.
OPT.fixTh   = false;
OPT.thFixed = th0;         % m, used only when OPT.fixTh is true

if OPT.fixTh
    params0 = [a0, h0, w0];
else
    params0 = [a0, h0, w0, th0];
end

%% ===================== TUNABLE KNOBS =====================
OPT.r1 = 10e-9;            % fillet radius 1 -- held fixed
OPT.r2 = 10e-9;            % fillet radius 2 -- held fixed

% --- solver backend -------------------------------------------------------
% Selects what actually produces the band structure. A STRING, not a boolean
% dryRun flag: three tiers do not fit in a boolean, and a name self-documents
% everywhere it is printed, logged or saved instead of becoming an anonymous 1.
%   'comsol'    - the real solveBands call. The only real physics.
%   'surrogate' - analytic fake bands. Runs the whole loop in milliseconds so
%                 the loop itself can be debugged without COMSOL.
%   'stub'      - no bands at all; every evaluation returns NaN.
% An unrecognised name is an ERROR (assertKnownBackend), never a silent
% fall-through: defaulting to 'comsol' would spend hours on a typo, and
% defaulting to 'surrogate' would hand back a folder of fabricated results
% that look exactly like a study.
OPT.solverBackend = 'comsol';
% OPT.solverBackend = 'surrogate';

assertKnownBackend(OPT.solverBackend);
OPT.isDryRun = ~strcmp(OPT.solverBackend, 'comsol');

% --- dry-run options (ignored unless OPT.isDryRun) ------------------------
OPT.dryRunDelay     = 0;    % artificial seconds per evaluation. 0 = as fast
                            % as possible. Raise it to watch the loop live.
OPT.dryRunFailEvery = 0;    % 0 = off. N > 0 makes every Nth evaluation fail
                            % deterministically, exercising the catch path.
OPT.dryRunPrefname  = 'DRYRUN';   % stamped onto synthetic output

% --- frequency targets ----------------------------------------------------
% NOTE ON SCALE. A 160 nm cell in diamond is a TENS-of-GHz structure, not a
% 10 GHz one: with a shear velocity around 12000 m/s the Bragg frequency is
% v/(2a) = 12000/(2*160e-9) ~ 37 GHz. The earlier 10 GHz target and 30 GHz
% ceiling here were carried over from the nanobeam work and would have thrown
% away every gap this cell actually has. Re-check these if you change a.
OPT.targetFreq = 35e9;     % Hz, mechanical midgap target
OPT.sigma      = 8e9;      % Hz, Gaussian width of the frequency penalty

% Smallest solid wall / void width [m].
%
% 10 nm, NOT the 50 nm used elsewhere in this repo (CreateNanobeamGeom's
% lithography guards). This cell is designed at a far finer scale and a 50 nm
% limit rejects the nominal design outright.
%
% The binding solid feature is the wall between the voids of ADJACENT cells,
% a - h = 160 - 140 = 20 nm for the starting design -- twice this limit. An
% earlier version of isCrossFabricable constrained (a-h)/2, the half-distance
% to the cell boundary, which is not a physical feature; that made the nominal
% appear to sit exactly on the limit with no room to grow h or shrink a. It
% now has 10 nm of margin on the wall and 15 nm on the fillet.
OPT.minFeature = 10e-9;
OPT.maxFreq    = 120e9;    % ignore gaps above this (spurious high-order bands)

OPT.kpts   = 5;            % k-points EXCLUDING gamma
OPT.nbands = 18;

% Penalty returned for a design that is not worth solving. Must exceed any
% attainable fitness; calFitness returns <= 0, so any positive value works.
OPT.penalty = 100;

OPT.runTag = 'cross_optimize_run1';

%% ===================== ITERATION LOG =====================
currentDate = datestr(now, 'mmddyyyy'); %#ok<TNOW1,DATST>

% PROVENANCE GUARD 1 of 3: synthetic results live in their own dated folder,
% DRYRUN_<backend>_<date>, so a real run never even looks in the directory the
% fabricated files are in. Mixing the two would permanently corrupt the one
% thing this directory exists to produce.
if OPT.isDryRun
    datedFolder   = [OPT.dryRunPrefname, '_', OPT.solverBackend, '_', currentDate];
    OPT.filePrefix = [OPT.dryRunPrefname, '_'];
else
    datedFolder   = currentDate;
    OPT.filePrefix = '';
end

% PROVENANCE GUARD 2 of 3: a DRYRUN filename prefix on every synthetic
% artifact, so the two could not collide even if the folders were merged by
% hand. The log, the per-evaluation folders and any saved bands all go through
% OPT.filePrefix.
OPT.rootLoc = [fullfile('.', 'test', OPT.runTag, datedFolder), filesep];
if ~exist(OPT.rootLoc, 'dir')
    mkdir(OPT.rootLoc);
end
OPT.itrPath = [OPT.rootLoc, OPT.filePrefix, 'optimization_', OPT.runTag, ...
               '_', currentDate, '.txt'];

if OPT.isDryRun
    fprintf(['*** DRY RUN: backend = ''%s''. The numbers below are ', ...
             'FABRICATED and\n*** must not be used to judge a design. ', ...
             'Output -> %s\n\n'], OPT.solverBackend, OPT.rootLoc);
end

% Create the log and header only when the file does not exist yet. The
% reference script tied this to the DIRECTORY existing, so a second run on the
% same day skipped the header and appended to nothing.
if ~exist(OPT.itrPath, 'file')
    itr = fopen(OPT.itrPath, 'wt+');
    fprintf(itr, ['eval\ta\th\tw\tth\tmidGapMech\tgapSizeMech\t', ...
                  'gapFracMech\tfitness\tstatus\r\n']);
    fclose(itr);
end

%% ===================== DRIVER =====================
% Evaluation counter shared with the objective. Each evaluation gets its own
% output folder: solveBands SKIPS the solve when <fileBase>_bds.mat already
% exists in datLoc and then returns a ds with only .sym/.asym and NO .full
% (solveBands.m:551-556). An optimizer that revisits a design would read
% ds.full off that stub and error, so every call must land somewhere fresh.
cross_optimize([], OPT, true);   % reset the persistent evaluation counter

%% debug: one evaluation before committing to a search
fitness = cross_optimize(params0, OPT);
fprintf('\nSingle-evaluation fitness = %.6g\n', fitness);

% %% run the optimization
% options = optimset('PlotFcns', @optimplotfval, 'Display', 'iter', ...
%                    'MaxFunEvals', 100);
% fun = @(x) cross_optimize(x, OPT);
% [xBest, fBest] = fminsearch(fun, params0, options);
% fprintf('\nBest: a=%.1fnm h=%.1fnm w=%.1fnm  fitness=%.6g\n', ...
%     xBest(1)*1e9, xBest(2)*1e9, xBest(3)*1e9, fBest);
% if ~OPT.fixTh; fprintf('  th=%.1fnm\n', xBest(4)*1e9); end


%% ========================================================================
function fitness = cross_optimize(params, OPT, doReset)
%CROSS_OPTIMIZE  One objective evaluation: geometry -> COMSOL -> fitness.
%
%   cross_optimize([], OPT, true) resets the evaluation counter and returns.
%   The driver calls that first, because the counter is PERSISTENT and the
%   script's own `clear all` cannot reach a persistent inside a function of the
%   file currently executing: without the reset a second run in the same MATLAB
%   session carries on numbering from where the last one stopped, so the log
%   rows and the per-evaluation folder names no longer start at 1.

persistent nEval
if nargin >= 3 && doReset
    nEval   = 0;
    fitness = 0;
    return;
end
if isempty(nEval); nEval = 0; end
nEval = nEval + 1;

%% ---- unpack the design vector ----
% Length depends on OPT.fixTh, so read th from the vector only when it is
% actually in there. Indexing params(4) unconditionally would silently read
% past the end -- or worse, pick up whatever a caller happened to append.
P.a  = params(1);      % square cell side
P.h  = params(2);      % cross arm length
P.w  = params(3);      % cross arm width
if OPT.fixTh
    P.th = OPT.thFixed;    % held, not searched
else
    P.th = params(4);      % slab thickness (z)
end
P.r1 = OPT.r1;         % fillet radii, fixed
P.r2 = OPT.r2;

%% ---- fixed configuration (mirrors test_CrossUnitCell.m) ----
P.xsect      = 'rect';
P.beamMat    = 'diamond';
P.celltype   = 'cross';
P.unitcell   = 'square';
P.nperiod    = 1;
P.holeatedge = 0;

P.kpts   = OPT.kpts;
P.nbands = OPT.nbands;

% Symmetry decomposition. solveBands gates on both of these before any solve
% and neither has a default (solveBands.m:138 / :186).
P.TwoSymPlanes  = 0;    % single symmetry plane -> 2 sectors, not 4
P.zSymCondition = 1;    % that plane is z
P.solveasym     = 1;

P.completeBandGaps = 1;
P.plotgeom    = 0;      % off during a search -- one figure per evaluation adds up
P.savedat     = 1;
P.savebndplot = 1;
P.saveplots   = 0;
P.saveMPH     = 0;
P.bandStruct_2D = 1;    % -> runBands_2D; also skips the P.wc title block,
                        % which references a field nothing ever assigns
                        % (solveBands.m:478, reachable only when this is 0)

P.mbeveny  = 0;
P.mbevenz  = 1;
P.freq     = 0;         % 0 for a band-structure sweep
P.meshSize = 4;
P.fixed_bc = 0;

P.anisoMat = 1;
P.rxtal    = 45;
P.max_dof  = 3e6;

%% ---- fabrication tolerance: reject BEFORE spending a solve ----
[fabOK, fabWhy] = crossFabCheck(P, OPT.minFeature);
if ~fabOK
    fitness = OPT.penalty;
    fprintf(['eval %3d: a=%.1f h=%.1f w=%.1f th=%.1f%s nm  ', ...
             '[UNFABRICABLE: %s]\n'], nEval, P.a*1e9, P.h*1e9, ...
             P.w*1e9, P.th*1e9, repmat('*', 1, OPT.fixTh), fabWhy);
    logRow(OPT.itrPath, nEval, P, NaN, NaN, NaN, fitness, fabWhy);
    return;   % <-- the reference script set fitness here but fell through
              %     and solved anyway, overwriting it.
end

%% ---- unique output folder per evaluation ----
% See the note in the driver: a repeated design would otherwise hit
% solveBands' cache branch and come back without a .full field.
P.datLoc = [OPT.rootLoc, OPT.filePrefix, ...
    sprintf('eval_%04d_a%.0f_h%.0f_w%.0f_th%.0f', ...
    nEval, P.a*1e9, P.h*1e9, P.w*1e9, P.th*1e9), filesep];

fprintf('eval %3d: a=%.1f h=%.1f w=%.1f th=%.1f%s nm%s\n', ...
    nEval, P.a*1e9, P.h*1e9, P.w*1e9, P.th*1e9, ...
    repmat('*', 1, OPT.fixTh), ...
    repmat('  [SYNTHETIC]', 1, OPT.isDryRun));

%% ---- solve ----
status = 'ok';
midGap = NaN;  gapSize = NaN;  gapRat = NaN;
try
    bds = solveCrossBands(P, OPT, nEval);

    % PROVENANCE GUARD 3 of 3: a real run that is somehow handed synthetic
    % data errors out instead of silently using it. The marker is written by
    % the surrogate/stub backends and checked here on EVERY evaluation.
    if ~OPT.isDryRun && isfield(bds, 'isSynthetic') && bds.isSynthetic
        error('cross_optimize:syntheticInRealRun', ...
            ['A real run received SYNTHETIC band data. Refusing to use it. ', ...
             'Check OPT.solverBackend and the contents of %s'], P.datLoc);
    end

    % Guard the cache branch explicitly rather than trusting the folder to be
    % new: solveBands returns .sym/.asym only in that case.
    if ~isfield(bds, 'full')
        error('cross_optimize:noFull', ...
            ['solveBands returned no .full field -- it took its cache ', ...
             'branch, so this design was already solved in %s'], P.datLoc);
    end

    gaps = bds.full;
    % Drop gaps above maxFreq: high-order bands routinely produce large but
    % physically uninteresting gaps that would dominate the objective.
    keep = gaps.midGap <= OPT.maxFreq;
    gs   = gaps.gapSize(keep);
    mg   = gaps.midGap(keep);

    if isempty(gs)                       % isempty, not length(...)==0
        status  = 'no_gap';
        gapSize = 0;
        midGap  = OPT.targetFreq;        % neutral, keeps gapRat finite
        gapRat  = 0;
    else
        [gapSize, iBest] = max(gs);
        midGap = mg(iBest);
        if midGap > 0
            gapRat = gapSize / midGap;
        else
            gapRat = 0;                  % a zero midgap would give Inf/NaN
        end
    end
catch ME
    status  = 'failed';
    warning('cross_optimize:evalFailed', ...
        'eval %d failed: %s', nEval, ME.message);
    fitness = OPT.penalty;
    logRow(OPT.itrPath, nEval, P, midGap, gapSize, gapRat, fitness, status);
    return;
end

%% ---- fitness ----
fitness = calFitness(midGap, gapRat, OPT.targetFreq, OPT.sigma);

fprintf(['         midGap=%.3f GHz  gapSize=%.3f GHz  gapRat=%.4f  ', ...
         'fitness=%.6g  [%s]\n'], ...
    midGap*1e-9, gapSize*1e-9, gapRat, fitness, status);

logRow(OPT.itrPath, nEval, P, midGap, gapSize, gapRat, fitness, status);
end

%% ========================================================================
function bds = solveCrossBands(P, OPT, nEval)
%SOLVECROSSBANDS  Backend dispatcher: real solve, surrogate, or stub.
%
%   Validating the name a SECOND time here is deliberate. The check at
%   configuration time catches a typo before any path is derived from it; this
%   one catches an OPT that was hand-edited or reloaded from an old results
%   file afterwards.
assertKnownBackend(OPT.solverBackend);

% Deterministic pseudo-failures, to exercise the catch path in the objective
% without waiting for a real solver to misbehave. Keyed on the evaluation
% index so a repeated design fails reproducibly.
if OPT.isDryRun && OPT.dryRunFailEvery > 0 && mod(nEval, OPT.dryRunFailEvery) == 0
    error('cross_optimize:injectedFailure', ...
        'injected dry-run failure (dryRunFailEvery = %d)', OPT.dryRunFailEvery);
end

if OPT.isDryRun && OPT.dryRunDelay > 0
    pause(OPT.dryRunDelay);
end

switch OPT.solverBackend
    case 'comsol'
        bds = solveBands(P);

    case 'surrogate'
        bds = surrogateCrossBands(P, OPT);

    case 'stub'
        % No bands at all. Shape-compatible, every frequency NaN, so the
        % objective's no_gap / failure handling is exercised end to end.
        bds.sym  = struct('F', [], 'k_norm', []);
        bds.asym = struct('F', [], 'k_norm', []);
        bds.full = struct('F', [], 'midGap', [], 'gapSize', []);
        bds.isSynthetic = true;
end
end

%% ------------------------------------------------------------------------
function assertKnownBackend(backend)
%ASSERTKNOWNBACKEND  Reject an unrecognised OPT.solverBackend, loudly.
%
%   An error rather than a warning because there is no safe default. Falling
%   back to 'comsol' would spend hours of solve time on a typo; falling back
%   to 'surrogate' would hand back a folder of fabricated results that look
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
    error('cross_optimize:unknownBackend', ...
        ['Unknown OPT.solverBackend. Expected one of: ''comsol'', ', ...
         '''surrogate'', ''stub''.\nGot: %s\n', ...
         '  ''comsol''    - the real COMSOL solve (the only real physics)\n', ...
         '  ''surrogate'' - analytic fake bands, for debugging the loop\n', ...
         '  ''stub''      - no bands at all, every evaluation returns NaN'], got);
end
end

%% ========================================================================
function f = calFitness(midGap, gapRat, targetFreq, sigma)
%CALFITNESS  Objective for the mechanical bandgap. LOWER is better.
%
%   f = -( gapRat * exp(-((targetFreq - midGap)/sigma)^2) )
%
%   Negated because fminsearch minimizes, matching the sign convention in
%   boomerang_optimize_sweep_diamond.m. The frequency term is Gaussian, so
%   sigma is the 1/e half-width: landing sigma away from the target costs a
%   factor 1/e, 2*sigma costs 1/e^4.
freqPenalty = exp(-((targetFreq - midGap)/sigma)^2);
f = -(gapRat * freqPenalty);
end

%% ========================================================================
function [ok, why] = crossFabCheck(P, minFeature)
%CROSSFABCHECK  Scalar wrapper over the shared, vectorized isCrossFabricable.
%
%   Kept as a thin adapter rather than a second implementation: bayesopt_cross
%   needs the vectorized form for its XConstraintFcn, and two copies of a
%   lithography rule set are exactly how the two scripts would come to disagree
%   about which designs are buildable.
[ok, ~, whyCell] = isCrossFabricable(P.a, P.h, P.w, P.th, minFeature, ...
                                     P.r1, P.r2);
why = whyCell{1};
end

%% ========================================================================
function logRow(itrPath, nEval, P, midGap, gapSize, gapRat, fitness, status)
%LOGROW  Append one evaluation. Opened and closed per row so the log survives
%   a mid-run crash -- the whole point of keeping it separate from the .mat.
try
    itr = fopen(itrPath, 'at+');
    if itr < 0; return; end
    fprintf(itr, '%d\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.6g\t%s\r\n', ...
        nEval, P.a, P.h, P.w, P.th, midGap, gapSize, gapRat, fitness, status);
    fclose(itr);
catch
    % Logging must never abort a study that is already COMSOL-hours deep.
end
end
