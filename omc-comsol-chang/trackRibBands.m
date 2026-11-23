function res = trackRibBands(designs, opts)
%TRACKRIBBANDS  Follow the two gap-edge bands of the rib cell across designs.
%
%   res = trackRibBands(designs)
%   res = trackRibBands(designs, opts)
%
%   Standalone tool (not part of any optimizer). Given an ORDERED list of
%   rib unit-cell designs, e.g. a sweep of t, s, L or a, it identifies the
%   band below and the band above the gap at the X point (k = pi/a) in the
%   first design and follows the same two modes from design to design by
%   their mode shape, reporting each mode's frequency, the gap between them,
%   the match confidence and a status per design.
%
%   designs   struct array or cell array, in sweep order. Each element is
%       - a P struct (as passed to solveBands): solved here, with mode export
%       - a ds struct returned by solveBands with P.exportModes set
%       - the path of a saved <fileBase>_bds.mat of such a solve
%     The three forms can be mixed in a cell array. Solved results must
%     carry the modes (P.exportModes = struct('k','all') recommended) and
%     their P (runBands stores it in the sector).
%
%   opts (all optional)
%     dryRun      false. true: no COMSOL; every solve (designs given as P,
%                 and all bridge solves) calls surrogateRibBandsTracked
%     rootDir     output folder. Default: a new time-stamped folder
%                 test/ribBandTrack/run_<yyyymmdd_HHMMSS>
%     sector      field of ds to track in. Default: 'sym' if present, else
%                 'symy_symz' (y-even, z-even)
%     initIndex   [lo up] sorted X-point columns of the two bands in the
%                 first design; or
%     initFreq    [fLo fUp] (Hz): the X-point modes nearest these; or
%     targetFreq  (Hz) the two bands bounding the band gap nearest this
%                 frequency. Default: P.freq of the first design
%     vars        design variables interpolated for bridge solves, a subset
%                 of {'a','t','s','L','th','d'}. Default: those that differ
%                 between the designs
%     maxBridge   2. Levels of bridging; 0 disables it (see below)
%     refineK     false. true: re-solve a design with 2*kpts when a
%                 k-connection overlap is below kConfMin (0.5)
%     exportK     'all'. k-points whose fields are exported in solves made
%                 here ('X' is enough for the identification itself)
%     grid        ribModeGrid point counts (struct; default grid if empty)
%     useParity   true. x-mirror parity veto at X. Valid for rxtal = 0 (all
%                 rib scripts) and 45; set false for other orientations
%     plot        false. true: summary figure (plotRibBandTrack), saved as
%                 track_summary.png
%     surrogate   options of surrogateRibBandsTracked for the dry run
%     macHi, macLo, marginMin, fTol, ...   thresholds of matchBandModes
%
%   METHOD. The reference is the last TRUSTED design (status init, confident
%   or bridged). Each new design's X-point modes are matched to it with
%   matchBandModes. If the match of either band is 'ambiguous', an
%   intermediate design is solved at the midpoint of the two designs and the
%   bands are tracked reference -> midpoint -> design, recursively up to
%   opts.maxBridge levels (at most 2^maxBridge - 1 extra solves per design).
%   If it is still ambiguous the design is marked 'ambiguous' (the best
%   guess is reported) and the run continues from the same reference.
%   Tracking through a bridged avoided crossing follows the continuous
%   eigen-branch, so the mode's character can change; macGold_* (overlap
%   with the first design's modes) shows that.
%
%   STATUS per design: init | confident | bridged | ambiguous | lost |
%   failed (the solve or the mode data failed; the message is in 'error').
%
%   res.table   one row per design: design, the design variables (nm / deg),
%               source, status, status_lo, status_up, f_lo, f_up (Hz, at X),
%               gapX (f_up - f_lo), bandLoMax, bandUpMin, gapBand (extremes
%               over k of the two frequency levels the modes sit in),
%               nBetween (other bands between them at X), idx_*, mac_*,
%               margin_*, macGold_*, kConf, step (normalized distance to the
%               reference), nBridge, nRefine, truth_* (dry run), folder, error
%   res.designs per-design details (P, x, v, cand, S, F, R, ...)
%   res.bridges one entry per bridge solve
%   res.nBridgeSolves, res.nRefineSolves, res.nSolves
%   res.vars, res.opts, res.rootDir, res.csvFile, res.matFile
%   The table is also written to <rootDir>/track_table.csv and res to
%   <rootDir>/track_results.mat, both rewritten after every design.
%
%   Each solve made here gets its own folder under rootDir (design_###,
%   bridge_####, ..._k2), so solveBands never reuses a cached result.
%
%   Example (dry run):
%     P = struct('t',100e-9,'L',370e-9);
%     Ps = arrayfun(@(s) setfield(P,'s',s), (103:10:193)*1e-9);
%     res = trackRibBands(Ps, struct('dryRun',true,'targetFreq',11.2e9,'plot',true));

if nargin < 2 || isempty(opts), opts = struct(); end
opts = fillDefaults(opts);
if ischar(designs) || isstring(designs)
    D = cellstr(designs);
elseif isstruct(designs)
    D = num2cell(designs(:)');
elseif iscell(designs)
    D = designs(:)';
else
    error('trackRibBands:badDesigns', ...
        'designs must be a struct array or a cell array of P structs, ds structs or file names.');
end
N = numel(D);
if N == 0
    error('trackRibBands:badDesigns', 'designs is empty.');
end

%% read the designs (nothing is solved yet)
nodes = repmat(emptyNode(), 1, N);
for i = 1:N
    nodes(i) = readDesign(D{i}, opts.sector, i);
end
allVars = {'a','t','s','L','th','d'};
X = zeros(N, numel(allVars));
for i = 1:N
    X(i,:) = designVector(nodes(i).P, allVars);
end
if isempty(opts.vars)
    span = max(X,[],1) - min(X,[],1);
    opts.vars = allVars(span > 1e-9*max(abs(X),[],1));
else
    opts.vars = cellstr(opts.vars);
    if ~all(ismember(opts.vars, allVars))
        error('trackRibBands:badVars', ...
            'opts.vars must be a subset of {%s}.', strjoin(allVars, ', '));
    end
end
[~, col] = ismember(opts.vars, allVars);
X = X(:, col);
xLo = min(X,[],1);
xSpan = max(X,[],1) - xLo;
xSpan(xSpan <= 0) = 1;
toV = @(x) (x - xLo)./xSpan;        % normalized design coordinates, [0, 1]

%% output folder: unique per run
if isempty(opts.rootDir)
    base = fullfile('.','test','ribBandTrack', ...
        ['run_', char(datetime('now','Format','yyyyMMdd_HHmmss'))]);
    opts.rootDir = base;  q = 1;
    while exist(opts.rootDir, 'dir')
        q = q + 1;
        opts.rootDir = sprintf('%s_%d', base, q);
    end
end
if ~exist(opts.rootDir, 'dir'), mkdir(opts.rootDir); end

%% track
ref = [];  gold = [];  rows = [];
bridges = repmat(emptyBridge(), 1, 0);
nSolve = 0;  nBridgeSolve = 0;  nRefineSolve = 0;  bridgeCount = 0;  cur = 0;
res = struct();
fprintf('trackRibBands: %d designs, variables {%s}, output %s\n', ...
    N, strjoin(opts.vars, ', '), opts.rootDir);
for i = 1:N
    cur = i;
    nd = nodes(i);
    nd.x = X(i,:);  nd.v = toV(nd.x);
    nB = 0;  step = NaN;  mg = [NaN NaN];  kConf = NaN;
    g = gapInfo([], [NaN NaN], [NaN NaN]);
    solved = false;
    try
        nd = prepare(nd, fullfile(opts.rootDir, sprintf('design_%03d', i)));
        solved = true;
    catch ME
        R = blankResult('failed');
        nd.status = 'failed';  nd.err = ME.message;
        warning('trackRibBands:designFailed', 'Design %d failed: %s', i, ME.message);
    end
    if solved
        if isempty(ref)
            nd.tIdx = pickInitialBands(nd, opts);
            R = blankResult('init');
            R.idx = nd.tIdx;  R.f = nd.cand.f(nd.tIdx);
            R.mac(:) = 1;  R.mac2(:) = 0;  R.margin(:) = 1;  R.ok = true;
            nd.status = 'init';
            gold = struct('U', nd.cand.U(:,nd.tIdx), 'w', nd.cand.w);
            ref = nd;  step = 0;
        else
            step = norm(nd.v - ref.v);
            [R, nB] = resolveStep(ref, nd, 0);
            nd.tIdx = R.idx;
            nd.status = overallStatus(R);
            if R.ok, ref = nd; end
        end
        for q = 1:2
            if isfinite(R.idx(q))
                mg(q) = bandOverlap(gold.U(:,q), nd.cand.U(:,R.idx(q)), ...
                    {gold.w, nd.cand.w});
            end
        end
        g = gapInfo(nd.S, R.idx, R.f);
        kConf = min(nd.S.conf, [], 'omitnan');
        if isempty(kConf), kConf = NaN; end
    end
    nd.R = R;
    nodes(i) = nd;
    rows = [rows; makeRow(i, nd, opts.vars, R, g, mg, kConf, step, nB)]; %#ok<AGROW>
    fprintf('  design %3d  %-9s  f_lo = %8.4f GHz  f_up = %8.4f GHz  mac = %.2f / %.2f  bridges = %d\n', ...
        i, nd.status, R.f(1)*1e-9, R.f(2)*1e-9, R.mac(1), R.mac(2), nB);
    writeOutputs();
end

st = {rows.status};
fprintf(['trackRibBands done: %d designs (%d ambiguous, %d lost, %d failed), ' ...
    '%d bridge solves, %d refine solves.\n  table: %s\n'], N, ...
    nnz(strcmp(st,'ambiguous')), nnz(strcmp(st,'lost')), nnz(strcmp(st,'failed')), ...
    nBridgeSolve, nRefineSolve, res.csvFile);
if opts.plot
    fig = plotRibBandTrack(res);
    saveas(fig, fullfile(opts.rootDir, 'track_summary.png'));
end

%% ===================== NESTED FUNCTIONS =====================
    function nd = prepare(nd, folder)
        % Make sure the node has sampled modes: solve it if it came as a P
        % struct, then build the tracker input. Optionally refine in k.
        if isempty(nd.ds)
            nd.ds = solveOne(nd.P, folder);
            if ~opts.dryRun, nd.folder = folder; end
        end
        sec = pickSector(nd.ds, opts.sector);
        [nd.cand, nd.S] = bandTrackCandidate(sec);
        if opts.refineK && any(nd.S.conf < opts.kConfMin)
            Pk = nd.P;
            Pk.kpts = 2*(size(sec.F,1) - 1);
            nRefineSolve = nRefineSolve + 1;
            nd.nRefine = 1;
            nd.ds = solveOne(Pk, [folder, '_k2']);
            sec = pickSector(nd.ds, opts.sector);
            [nd.cand, nd.S] = bandTrackCandidate(sec);
        end
        nd.cand.truth = nan(1, numel(nd.cand.f));
        if isfield(sec, 'truth')                    % dry-run ground truth
            nd.cand.truth = sec.truth(end, nd.cand.order);
        end
        nd.F = sec.F;  nd.k_norm = sec.k_norm;
        nd.ds = [];                                 % keep res small
    end

    function ds = solveOne(P, folder)
        % One band-structure solve with mode export, in its own folder.
        nSolve = nSolve + 1;
        if opts.dryRun
            ds = surrogateRibBandsTracked(P, opts.surrogate);
            return
        end
        P.datLoc = [folder, filesep];
        if isfield(P,'fileBase'), P = rmfield(P,'fileBase'); end
        P.exportModes = struct('k', opts.exportK, 'n', opts.grid);
        P.savedat = 1;
        ds = solveBands(P);
        removeStrayBackslashDir(P.datLoc);
        if ~hasSector(ds, opts.sector)
            % solveBands found a cached result and returned empty sectors
            f = dir(fullfile(folder, '*_bds.mat'));
            if isempty(f)
                error('trackRibBands:noResult', ...
                    'solveBands returned no bands and %s holds no _bds.mat.', folder);
            end
            Ld = load(fullfile(folder, f(1).name), 'ds');
            ds = Ld.ds;
        end
    end

    function [R, nB] = resolveStep(rf, nd, depth)
        % Match nd to the trusted node rf; bridge through the midpoint
        % design while the match is ambiguous and depth < maxBridge.
        R = matchBandModes(rf.cand, rf.tIdx, nd.cand, opts);
        R.bridged = false;
        nB = 0;
        if ~R.anyAmbiguous || depth >= opts.maxBridge || norm(nd.v - rf.v) < 1e-9
            return
        end
        bridgeCount = bridgeCount + 1;
        id = bridgeCount;
        mid = emptyNode();
        mid.source = 'bridge';
        entry = emptyBridge();
        entry.id = id;  entry.design = cur;  entry.depth = depth + 1;
        try
            mid.P = midDesign(rf.P, nd.P, opts.vars);
            mid.P.plotgeom = 0;  mid.P.savebndplot = 0;  mid.P.saveplots = 0;
            checkRibDesign(mid.P);
            mid.x = designVector(mid.P, opts.vars);
            mid.v = toV(mid.x);
            entry.x = mid.x;  entry.v = mid.v;
            nBridgeSolve = nBridgeSolve + 1;
            nB = 1;
            mid = prepare(mid, fullfile(opts.rootDir, sprintf('bridge_%04d', id)));
        catch ME
            warning('trackRibBands:bridgeFailed', ...
                'Bridge solve %d (design %d) failed: %s', id, cur, ME.message);
            entry.status = 'failed';  entry.err = ME.message;
            bridges(id) = entry;
            return
        end
        [Rm, n1] = resolveStep(rf, mid, depth + 1);
        nB = nB + n1;
        entry.status = overallStatus(Rm);
        entry.idx = Rm.idx;  entry.f = Rm.f;  entry.mac = Rm.mac;
        entry.folder = mid.folder;
        bridges(id) = entry;
        if ~Rm.ok, return; end
        mid.tIdx = Rm.idx;
        [R, n2] = resolveStep(mid, nd, depth + 1);
        nB = nB + n2;
        R.bridged = true;
    end

    function writeOutputs()
        res.table = struct2table(rows, 'AsArray', true);
        res.designs = nodes;
        res.bridges = bridges;
        res.nBridgeSolves = nBridgeSolve;
        res.nRefineSolves = nRefineSolve;
        res.nSolves = nSolve;
        res.vars = opts.vars;
        res.opts = opts;
        res.rootDir = opts.rootDir;
        res.csvFile = fullfile(opts.rootDir, 'track_table.csv');
        res.matFile = fullfile(opts.rootDir, 'track_results.mat');
        writetable(res.table, res.csvFile);
        save(res.matFile, 'res');
    end
end

%% ===================== LOCAL HELPERS =====================
function opts = fillDefaults(opts)
def = struct('dryRun',false, 'rootDir','', 'sector','', 'targetFreq',[], ...
    'initIndex',[], 'initFreq',[], 'vars',{{}}, 'maxBridge',2, ...
    'refineK',false, 'kConfMin',0.5, 'exportK','all', 'grid',struct(), ...
    'useParity',true, 'plot',false, 'surrogate',struct());
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts, fn{k}), opts.(fn{k}) = def.(fn{k}); end
end
if ~isfield(opts.surrogate, 'n'), opts.surrogate.n = opts.grid; end
end

% -------------------------------------------------------------------------

function nd = emptyNode()
nd = struct('P',[], 'ds',[], 'cand',[], 'S',[], 'F',[], 'k_norm',[], ...
    'x',[], 'v',[], 'source','', 'folder','', 'tIdx',[NaN NaN], ...
    'nRefine',0, 'R',[], 'status','', 'err','');
end

function b = emptyBridge()
b = struct('id',NaN, 'design',NaN, 'depth',NaN, 'x',[], 'v',[], ...
    'status','', 'idx',[NaN NaN], 'f',[NaN NaN], 'mac',[NaN NaN], ...
    'folder','', 'err','');
end

function R = blankResult(status)
R = struct('idx',[NaN NaN], 'f',[NaN NaN], 'mac',[NaN NaN], ...
    'mac2',[NaN NaN], 'margin',[NaN NaN], 'ok',false, 'anyAmbiguous',false, ...
    'bridged',false);
R.status = {status, status};
end

% -------------------------------------------------------------------------

function nd = readDesign(d, sector, i)
% Classify one input: file name, solved ds struct, or P struct.
nd = emptyNode();
if ischar(d) || (isstring(d) && isscalar(d))
    file = char(d);
    if ~isfile(file)
        error('trackRibBands:noFile', 'Design %d: file not found: %s', i, file);
    end
    Ld = load(file);
    if ~isfield(Ld, 'ds')
        error('trackRibBands:noDs', 'Design %d: %s holds no variable ds.', i, file);
    end
    nd.ds = Ld.ds;  nd.source = 'file';  nd.folder = file;
elseif isstruct(d) && isscalar(d) && any(isfield(d, {'sym','symy_symz','full'}))
    nd.ds = d;  nd.source = 'given';
elseif isstruct(d) && isscalar(d)
    nd.P = d;  nd.source = 'solve';
else
    error('trackRibBands:badDesigns', ...
        'Design %d is not a P struct, a ds struct or a file name.', i);
end
if ~isempty(nd.ds)
    sec = pickSector(nd.ds, sector);
    if isfield(sec, 'P')
        nd.P = sec.P;
    elseif isfield(nd.ds, 'P')
        nd.P = nd.ds.P;
    else
        error('trackRibBands:noP', ...
            'Design %d: the solved result does not carry its P struct.', i);
    end
end
Pv = ribDefaultParams(nd.P);
if ~strcmp(Pv.celltype, 'rib')
    error('trackRibBands:celltype', ...
        'Design %d has celltype ''%s''; only ''rib'' is supported.', i, Pv.celltype);
end
end

% -------------------------------------------------------------------------

function sec = pickSector(ds, name)
% The symmetry sector to track in.
if isempty(name)
    if isfield(ds,'sym') && ~isempty(ds.sym)
        name = 'sym';
    elseif isfield(ds,'symy_symz') && ~isempty(ds.symy_symz)
        name = 'symy_symz';
    else
        error('trackRibBands:noSector', ...
            ['The result has neither ds.sym nor ds.symy_symz (a cached ' ...
             'solveBands call returns empty sectors; pass the _bds.mat file).']);
    end
end
if ~isfield(ds, name) || isempty(ds.(name))
    error('trackRibBands:noSector', 'The result has no sector ''%s''.', name);
end
sec = ds.(name);
end

function tf = hasSector(ds, name)
tf = true;
try
    pickSector(ds, name);
catch
    tf = false;
end
end

% -------------------------------------------------------------------------

function x = designVector(P, vars)
% Effective values of the design variables (L through ribLength).
Pv = ribDefaultParams(P);
x = zeros(1, numel(vars));
for k = 1:numel(vars)
    if strcmp(vars{k}, 'L')
        x(k) = ribLength(Pv, Pv.a);
    else
        x(k) = Pv.(vars{k});
    end
end
end

% -------------------------------------------------------------------------

function Pm = midDesign(Pa, Pb, vars)
% The design half-way between two designs. All other fields come from Pa.
% (The midpoint is the same in raw and in normalized coordinates.)
Pm = Pa;
A = ribDefaultParams(Pa);  B = ribDefaultParams(Pb);
for k = 1:numel(vars)
    f = vars{k};
    if strcmp(f, 'L'), continue; end
    Pm.(f) = (A.(f) + B.(f))/2;
end
[La, srcA] = ribLength(A, A.a);
[Lb, srcB] = ribLength(B, B.a);
if ~(strcmp(srcA,'rule') && strcmp(srcB,'rule'))
    % Explicit half-length: drop the width and any builder-written record so
    % ribLength sees exactly this value.
    Pm.L = (La + Lb)/2;
    drop = {'w','ribDerived'};
    for k = 1:numel(drop)
        if isfield(Pm, drop{k}), Pm = rmfield(Pm, drop{k}); end
    end
end
if isfield(Pa,'freq') && isfield(Pb,'freq')
    Pm.freq = (Pa.freq + Pb.freq)/2;
end
if isfield(Pm,'fileBase'), Pm = rmfield(Pm,'fileBase'); end
end

% -------------------------------------------------------------------------

function checkRibDesign(P)
% The checks the builder runs; errors if the design cannot be built.
Pv = ribDefaultParams(P);
ribLength(Pv, Pv.a);
validateRibParams(Pv, 'unitcell');
end

% -------------------------------------------------------------------------

function tIdx = pickInitialBands(nd, opts)
% Sorted X-point columns [lower upper] of the two bands to follow.
f = nd.cand.f;  S = nd.S;  nb = numel(f);
if ~isempty(opts.initIndex)
    tIdx = sort(reshape(double(opts.initIndex), 1, []));
    if numel(tIdx) ~= 2 || tIdx(1) == tIdx(2) || any(tIdx < 1 | tIdx > nb | tIdx ~= round(tIdx))
        error('trackRibBands:badInitIndex', ...
            'opts.initIndex must be two different sorted X-point columns in 1..%d.', nb);
    end
    return
end
if ~isempty(opts.initFreq)
    if numel(opts.initFreq) ~= 2
        error('trackRibBands:badInitFreq', 'opts.initFreq must be [fLower fUpper] in Hz.');
    end
    [~, j1] = min(abs(f - min(opts.initFreq)));
    [~, j2] = min(abs(f - max(opts.initFreq)));
    if j1 == j2
        error('trackRibBands:badInitFreq', ...
            'Both opts.initFreq values select the same X-point mode (%.4f GHz).', f(j1)*1e-9);
    end
    tIdx = [j1 j2];
    return
end
target = opts.targetFreq;
if isempty(target) && isfield(nd.P,'freq') && nd.P.freq > 0
    target = nd.P.freq;
end
if isempty(target)
    error('trackRibBands:noTarget', ...
        'Set opts.targetFreq, opts.initFreq or opts.initIndex to choose the two bands.');
end
% The band gap (over all k) nearest the target, from the aligned levels.
best = Inf;  levLo = NaN;
for n = 1:size(S.Flev, 2) - 1
    if ~(S.complete(n) && S.complete(n+1)), continue; end
    lo = max(S.Flev(:,n));  hi = min(S.Flev(:,n+1));
    if hi <= lo, continue; end
    d = abs((lo + hi)/2 - target);
    if d < best, best = d;  levLo = n; end
end
j1 = [];  j2 = [];
if ~isnan(levLo)
    j1 = find(S.levAnchor == levLo, 1);
    j2 = find(S.levAnchor == levLo + 1, 1);
end
if isempty(j1) || isempty(j2)
    % No band gap in the window: the two X-point modes around the target.
    j1 = find(f < target, 1, 'last');
    j2 = j1 + 1;
    if isempty(j1) || j2 > nb
        error('trackRibBands:noInitialBands', ...
            ['The first design has no band gap and opts.targetFreq (%.4g GHz) ' ...
             'is outside its X-point frequencies. Set opts.initIndex or opts.initFreq.'], ...
            target*1e-9);
    end
end
tIdx = [j1 j2];
end

% -------------------------------------------------------------------------

function g = gapInfo(S, idx, f)
% Gap between the two followed modes at X, and between the band extremes of
% the frequency levels they sit in.
g = struct('gapX',f(2) - f(1), 'bandLoMax',NaN, 'bandUpMin',NaN, ...
    'gapBand',NaN, 'nBetween',NaN);
if isempty(S) || any(~isfinite(idx)), return; end
lev = S.levAnchor(idx);
g.nBetween = abs(lev(2) - lev(1)) - 1;
if S.complete(lev(1)), g.bandLoMax = max(S.Flev(:,lev(1))); end
if S.complete(lev(2)), g.bandUpMin = min(S.Flev(:,lev(2))); end
g.gapBand = g.bandUpMin - g.bandLoMax;
end

% -------------------------------------------------------------------------

function st = overallStatus(R)
if R.ok
    if R.bridged, st = 'bridged'; else, st = 'confident'; end
elseif any(strcmp(R.status, 'lost'))
    st = 'lost';
else
    st = 'ambiguous';
end
end

% -------------------------------------------------------------------------

function row = makeRow(i, nd, vars, R, g, mg, kConf, step, nB)
row = struct('design', i);
for k = 1:numel(vars)
    if strcmp(vars{k}, 'd')
        row.d_deg = nd.x(k)*180/pi;
    else
        row.([vars{k}, '_nm']) = nd.x(k)*1e9;
    end
end
tr = [NaN NaN];
for q = 1:2
    if isfinite(R.idx(q)) && ~isempty(nd.cand)
        tr(q) = nd.cand.truth(R.idx(q));
    end
end
row.source = nd.source;
row.status = nd.status;
row.status_lo = R.status{1};   row.status_up = R.status{2};
row.f_lo = R.f(1);             row.f_up = R.f(2);
row.gapX = g.gapX;
row.bandLoMax = g.bandLoMax;   row.bandUpMin = g.bandUpMin;
row.gapBand = g.gapBand;       row.nBetween = g.nBetween;
row.idx_lo = R.idx(1);         row.idx_up = R.idx(2);
row.mac_lo = R.mac(1);         row.mac_up = R.mac(2);
row.margin_lo = R.margin(1);   row.margin_up = R.margin(2);
row.macGold_lo = mg(1);        row.macGold_up = mg(2);
row.kConf = kConf;             row.step = step;
row.nBridge = nB;              row.nRefine = nd.nRefine;
row.truth_lo = tr(1);          row.truth_up = tr(2);
row.folder = nd.folder;        row.error = nd.err;
end

% -------------------------------------------------------------------------

function removeStrayBackslashDir(datLoc)
% solveBands appends a literal '\' to datLoc. On macOS/Linux that creates an
% empty sub-folder named '\' inside the solve folder; remove it. (Same helper
% as in optimize_rib_unitCell.m.)
if ispc, return; end
stray = [datLoc, '\'];
if exist(stray, 'dir')
    try
        rmdir(stray);   % fails, harmlessly, if it is not empty
    catch
    end
end
end
