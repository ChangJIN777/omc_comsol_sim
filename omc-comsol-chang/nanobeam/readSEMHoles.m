function S = readSEMHoles(csvFile, opts)
%READSEMHOLES  Full-beam SEM hole table -> centre-first half-beam hole arrays.
%   S = readSEMHoles(csvFile) reads a CSV with header hole,a_nm,hx_nm,hy_nm
%   (values in nm, one row per hole, "hole" = position index along the beam:
%   any contiguous integer range, any row order). It finds the cavity centre
%   from the mirror symmetry of hx and hy, and returns the arrays that
%   CreateNanobeamGeom reads when P.useManualGeom = 1 (and, for an end taper,
%   P.useManualWvgMir = 1).
%
%   a_nm convention (option AConvention):
%     'next' (default): a_k = centre-to-centre distance from hole k to hole
%                       k+1. The last row's a is not used and may be NaN.
%     'prev':           a_k = distance from hole k-1 to hole k. The first
%                       row's a is not used and may be NaN.
%   The table is first turned into hole positions. The output a_hole then
%   follows CreateNanobeamGeom.m:88 (centre-first, a_hole(i) = spacing from
%   hole i-1 to hole i). a_hole(1) is NaN for holeatctr = 1, and the spacing
%   of the two central holes for holeatctr = 0.
%
%   Name-value options:
%     CenterHole  [] = auto | file hole index k (hole at the centre) or k+0.5
%                 (dielectric centre between holes k and k+1)
%     HoleAtCtr   [] = auto | 1 | 0. Restricts the auto search to integer or
%                 half-integer centres.
%     AConvention 'next' | 'prev'
%     Half        'average' (default): mean of the two halves over the
%                     region measured on both sides (symmetric model)
%                 'right' | 'left': one half as built (symmetric model)
%                 'both': right half in S, left half in S.L, for P.asymCav = 1
%                     with P.PL = left half (CreateNanobeamGeom_asym)
%     EndTaper    what to do with the end taper. The end taper is the holes
%                 measured on one side only, past the two-sided region.
%                 'drop' (default): discard them; pad with Mirror to NHoles
%                 'holes':  keep them as ordinary a_hole/hx_hole/hy_hole entries
%                 'wvgmir': return them as a_wvgmir/hx_wvgmir/hy_wvgmir
%                           (CreateNanobeamGeom.m:238-250, mirror side first)
%     EndTaperStart  [] | file hole index of the first end-taper hole. Holes
%                 between the two-sided region and it stay core holes.
%     Mirror      [] | [a hx hy] in m, used for padding. [] = the placeholder
%                 rows of the file (>= 3 identical rows at one end).
%     NHoles      [] | half-beam hole count, including the centre hole.
%                 Reached by padding with Mirror cells or by dropping outer
%                 cells. [] = the longest half in the file. A half that keeps
%                 its end taper ('holes'/'wvgmir') is never padded or cut.
%     MirrorTol   0.02. Relative tolerance used to estimate ndef.
%     Verbose     true. Prints the centre search, pair mismatch and segments.
%
%   Output S (lengths in m, arrays are columns):
%     a_hole, hx_hole, hy_hole, nholes, ndef     main half (right if 'both')
%     wvgmir, a_wvgmir, hx_wvgmir, hy_wvgmir     wvgmir = 0, [] when unused
%     L            left half in the same format ('both' only, else [])
%     centerHole, holeatctr, aConvention, half, endTaper
%     mirror       [a hx hy] used for padding ([] if none)
%     placeholder  file hole indices of the placeholder rows
%     pairs        table: mirrored pairs, left/right values
%     centreScan   table: candidate centres and their mismatch
%     aSymmetry    table: spacing mismatch for both a conventions
%     file         the file sorted by hole, plus x_m (position, centre = 0)
%
%   Example:
%     S = readSEMHoles('SEM_data/omc_v2_dose430_a0MN3w1.csv','Half','average');
%     P.useManualGeom = 1;  P.holeatctr = S.holeatctr;  P.nholes = S.nholes;
%     P.a_hole = S.a_hole;  P.hx_hole = S.hx_hole;  P.hy_hole = S.hy_hole;

arguments
    csvFile {mustBeTextScalar}
    opts.CenterHole double = []
    opts.HoleAtCtr double = []
    opts.AConvention {mustBeMember(opts.AConvention,{'next','prev'})} = 'next'
    opts.Half {mustBeMember(opts.Half,{'average','right','left','both'})} = 'average'
    opts.EndTaper {mustBeMember(opts.EndTaper,{'drop','holes','wvgmir'})} = 'drop'
    opts.EndTaperStart double = []
    opts.Mirror double = []
    opts.NHoles double = []
    opts.MirrorTol (1,1) double {mustBePositive} = 0.02
    opts.Verbose (1,1) logical = true
end
csvFile = char(csvFile);

%% read and validate the table
if ~isfile(csvFile)
    error('readSEMHoles:noFile','SEM hole table not found: %s',csvFile);
end
T = readtable(csvFile,'Delimiter',',','VariableNamingRule','preserve');
need = {'hole','a_nm','hx_nm','hy_nm'};
miss = setdiff(need,T.Properties.VariableNames,'stable');
if ~isempty(miss)
    error('readSEMHoles:columns',['%s: missing column(s) %s; the header must be ',...
        'hole,a_nm,hx_nm,hy_nm'],csvFile,strjoin(miss,', '));
end
for c = need
    if ~isnumeric(T.(c{1}))
        error('readSEMHoles:nonNumeric','%s: column %s has non-numeric entries',...
            csvFile,c{1});
    end
end
hole = T.hole(:);
if any(~isfinite(hole)) || any(hole ~= round(hole))
    error('readSEMHoles:holeIndex','%s: every "hole" entry must be a finite integer',...
        csvFile);
end
[hole,ord] = sort(hole);
T = T(ord,:);
if any(diff(hole) == 0)
    error('readSEMHoles:holeIndex','%s: hole index repeated: %s',csvFile,...
        mat2str(unique(hole(diff(hole) == 0))'));
end
if any(diff(hole) ~= 1)
    error('readSEMHoles:holeIndex','%s: hole indices must be contiguous; missing %s',...
        csvFile,mat2str(setdiff(hole(1):hole(end),hole)'));
end
N = numel(hole);
if N < 7
    error('readSEMHoles:tooShort','%s: %d holes; need at least 7',csvFile,N);
end
a  = T.a_nm(:)';  hx = T.hx_nm(:)';  hy = T.hy_nm(:)';
aUsed = true(1,N);
if strcmp(opts.AConvention,'next'), aUsed(N) = false; else, aUsed(1) = false; end
bad = ~isfinite(hx) | hx <= 0 | ~isfinite(hy) | hy <= 0 | ...
      (aUsed & (~isfinite(a) | a <= 0));
if any(bad)
    error('readSEMHoles:badValue',['%s: a/hx/hy must be finite and > 0 (only the ',...
        'unused end-row a may be NaN); check hole(s) %s'],csvFile,mat2str(hole(bad)'));
end

%% hole positions [nm] along the beam
x = holePositions(a,opts.AConvention);

%% placeholder rows: >= 3 identical rows at either end of the table
isPh = false(1,N);
rowsEq = @(i,j) hx(i) == hx(j) && hy(i) == hy(j) && ...
         (a(i) == a(j) || ~aUsed(i) || ~aUsed(j));   % unused end-row a is ignored
k = 1;  while k < N && rowsEq(k+1,1), k = k+1; end
if k >= 3, isPh(1:k) = true; end
k = N;  while k > 1 && rowsEq(k-1,N), k = k-1; end
if N-k+1 >= 3, isPh(k:N) = true; end

%% centre: minimise the left/right mismatch of hx and hy (measured rows only)
if isempty(opts.CenterHole)
    cand = (hole(1)+1):0.5:(hole(N)-1);
    if ~isempty(opts.HoleAtCtr)
        cand = cand((mod(cand,1) == 0) == logical(opts.HoleAtCtr));
    end
else
    cand = opts.CenterHole;
    if ~isscalar(cand) || mod(2*cand,1) ~= 0 || cand <= hole(1) || cand >= hole(N)
        error('readSEMHoles:center',...
            'CenterHole = %s must be an inner hole index k or k+0.5',mat2str(cand));
    end
    if ~isempty(opts.HoleAtCtr) && (mod(cand,1) == 0) ~= logical(opts.HoleAtCtr)
        error('readSEMHoles:center','CenterHole = %g contradicts HoleAtCtr = %d',...
            cand,opts.HoleAtCtr);
    end
end
hxRef = mean(hx(~isPh));  hyRef = mean(hy(~isPh));
scan = nan(numel(cand),5);                  % [centre nPairs rmsHx rmsHy score]
for ic = 1:numel(cand)
    [iL,iR] = pairRows(cand(ic)-hole(1)+1,N,isPh);
    scan(ic,1:2) = [cand(ic), numel(iL)];
    if numel(iL) >= 3 || ~isempty(opts.CenterHole)
        dhx = hx(iL) - hx(iR);  dhy = hy(iL) - hy(iR);
        scan(ic,3:5) = [rms(dhx), rms(dhy), sqrt(mean((dhx/hxRef).^2 + (dhy/hyRef).^2))];
    end
end
scan = scan(isfinite(scan(:,5)),:);
if isempty(scan)
    error('readSEMHoles:center',...
        '%s: no centre candidate has >= 3 measured mirror pairs',csvFile);
end
scan = sortrows(scan,5);
centreScan = array2table(scan,'VariableNames',...
    {'centre','nPairs','rmsHx_nm','rmsHy_nm','score'});
if size(scan,1) > 1 && scan(2,5) < 1.5*scan(1,5)
    warning('readSEMHoles:ambiguousCentre',['centre %g is only marginally better ',...
        'than %g (score %.4f vs %.4f); set CenterHole explicitly'],...
        scan(1,1),scan(2,1),scan(1,5),scan(2,5));
end
cH = scan(1,1);                             % file hole index of the centre
holeatctr = double(mod(cH,1) == 0);
cRow = cH - hole(1) + 1;                    % row (half-integer: between rows)

%% rows of each half, centre-first; region measured on both sides
if holeatctr
    rowsR = cRow:N;          rowsL = cRow:-1:1;
else
    rowsR = cRow+0.5:N;      rowsL = cRow-0.5:-1:1;
end
measR = firstPh(isPh(rowsR));               % measured entries per half
measL = firstPh(isPh(rowsL));
nBoth = min(measR,measL);                   % entries measured on both halves
if holeatctr
    a0 = NaN;
else
    a0 = x(rowsR(1)) - x(rowsL(1));         % spacing of the two central holes
end
spacR = [a0, abs(diff(x(rowsR)))];          % a_hole convention, nm
spacL = [a0, abs(diff(x(rowsL)))];

% end taper = measured entries past the two-sided region (one side only)
tapR = nBoth+1:measR;
tapL = nBoth+1:measL;
if ~isempty(opts.EndTaperStart)
    iR = find(hole(rowsR) == opts.EndTaperStart);
    iL = find(hole(rowsL) == opts.EndTaperStart);
    if ~isempty(iR) && iR > nBoth && iR <= measR
        tapR = iR:measR;   tapL = [];
    elseif ~isempty(iL) && iL > nBoth && iL <= measL
        tapL = iL:measL;   tapR = [];
    else
        error('readSEMHoles:endTaperStart',['EndTaperStart = %g is not a measured ',...
            'hole past the two-sided region (centre %g, %d two-sided entries)'],...
            opts.EndTaperStart,cH,nBoth);
    end
end

%% mirror cell used for padding [m]
mirror = opts.Mirror(:)';
if isempty(mirror) && any(isPh)
    iPh = find(isPh);
    if isPh(1) && isPh(N) && ~rowsEq(1,N)
        warning('readSEMHoles:twoMirrors',...
            'placeholder rows differ at the two ends; padding with hole %d',hole(1));
    end
    mirror = [a(iPh(1)), hx(iPh(1)), hy(iPh(1))]*1e-9;
end
if ~isempty(mirror) && (numel(mirror) ~= 3 || any(~isfinite(mirror)) || any(mirror <= 0))
    error('readSEMHoles:mirror','Mirror must be [a hx hy] in m, finite and > 0');
end

%% halves, centre-first, in nm
hR = splitHalf(rowsR,spacR,tapR);
hL = splitHalf(rowsL,spacL,tapL);
if isempty(opts.NHoles)
    nTarget = max(numel(hR.hx),numel(hL.hx));   % longest core in the file
else
    nTarget = opts.NHoles;
end
switch opts.Half
    case {'right','both'}
        main = hR;
    case 'left'
        main = hL;
    case 'average'
        if ~isempty(hR.ext) && ~isempty(hL.ext)
            error('readSEMHoles:average',['both halves have one-sided measured ',...
                'holes outside the end taper; use Half = ''right''/''left''/''both''']);
        end
        n = nBoth;
        main.a  = [(hR.a(1:n)  + hL.a(1:n))/2,  hR.extA,  hL.extA];
        main.hx = [(hR.hx(1:n) + hL.hx(1:n))/2, hR.extHx, hL.extHx];
        main.hy = [(hR.hy(1:n) + hL.hy(1:n))/2, hR.extHy, hL.extHy];
        main.a(1) = a0;
        if ~isempty(hR.tap), main.tap = hR.tap; else, main.tap = hL.tap; end
end
if strcmp(opts.Half,'both') && strcmp(opts.EndTaper,'drop') && ~isempty([tapR,tapL])
    warning('readSEMHoles:bothDrop',['Half = ''both'' with EndTaper = ''drop'' ',...
        'replaces the measured end taper by mirror cells; use ''holes'' or ''wvgmir'' ',...
        'for the as-built beam']);
end

S = finishHalf(main);
if strcmp(opts.Half,'both')
    S.L = finishHalf(hL);
else
    S.L = [];
end
S.centerHole = cH;
S.holeatctr = holeatctr;
S.aConvention = opts.AConvention;
S.half = opts.Half;
S.endTaper = opts.EndTaper;
S.mirror = mirror;
S.placeholder = hole(isPh)';
S.centreScan = centreScan;

% mirrored pairs measured on both sides (entry 1 of holeatctr = 1 is the centre)
ip = 1+holeatctr:nBoth;
S.pairs = table((1:numel(ip))',hole(rowsL(ip)),hole(rowsR(ip)),...
    hx(rowsL(ip))',hx(rowsR(ip))',hy(rowsL(ip))',hy(rowsR(ip))',...
    spacL(ip)',spacR(ip)','VariableNames',{'j','holeL','holeR','hxL_nm','hxR_nm',...
    'hyL_nm','hyR_nm','spacingL_nm','spacingR_nm'});

% a-convention check: symmetry of the spacings under both readings
conv = {'next','prev'};
aSym = nan(1,2);
for kc = 1:2
    xx = holePositions(a,conv{kc});
    aSym(kc) = rms(abs(diff(xx(rowsR(1:nBoth)))) - abs(diff(xx(rowsL(1:nBoth)))));
end
S.aSymmetry = table(conv',aSym','VariableNames',{'aConvention','rmsSpacingMismatch_nm'});
kUse = find(strcmp(conv,opts.AConvention));
if aSym(3-kUse) < 0.5*aSym(kUse)
    warning('readSEMHoles:aConvention',['spacings are more symmetric with ',...
        'AConvention = ''%s'' (rms %.2f nm vs %.2f nm); check what a_nm means'],...
        conv{3-kUse},aSym(3-kUse),aSym(kUse));
end

xCtr = interp1(1:N,x,cRow);
S.file = T;
S.file.x_m = (x(:) - xCtr)*1e-9;

if opts.Verbose
    printSummary(S,csvFile,measL,measR);
end

%% ---------------------------------------------------------------- nested
    function h = splitHalf(rows,spac,tap)
        % one half, nm: core = all but the end taper; ext = measured one-sided
        % core entries (past the two-sided region, before the taper)
        core = setdiff(1:numel(rows),tap,'stable');
        h.a  = spac(core);   h.hx = hx(rows(core));   h.hy = hy(rows(core));
        ext = setdiff(nBoth+1:firstPh(isPh(rows)),tap);
        h.ext = ext;
        h.extA = spac(ext);  h.extHx = hx(rows(ext));  h.extHy = hy(rows(ext));
        if isempty(tap)
            h.tap = [];
        else
            h.tap = struct('a',spac(tap),'hx',hx(rows(tap)),'hy',hy(rows(tap)));
        end
    end

    function out = finishHalf(h)
        % end-taper handling, padding/truncation, SI units, ndef estimate
        aH = h.a;  hxH = h.hx;  hyH = h.hy;  tp = h.tap;
        if ~isempty(tp) && strcmp(opts.EndTaper,'holes')
            aH = [aH, tp.a];  hxH = [hxH, tp.hx];  hyH = [hyH, tp.hy];
        end
        if isempty(tp) || strcmp(opts.EndTaper,'drop')
            n0 = numel(hxH);
            if nTarget > n0
                if isempty(mirror)
                    error('readSEMHoles:noMirror',['NHoles = %d needs %d padding ',...
                        'cells but the file has no placeholder rows; pass Mirror = ',...
                        '[a hx hy] in m'],nTarget,nTarget-n0);
                end
                mm = mirror*1e9;
                aH  = [aH,  mm(1)*ones(1,nTarget-n0)];
                hxH = [hxH, mm(2)*ones(1,nTarget-n0)];
                hyH = [hyH, mm(3)*ones(1,nTarget-n0)];
            elseif nTarget < n0
                if nTarget < nBoth
                    warning('readSEMHoles:truncate',['NHoles = %d cuts into the ',...
                        'two-sided measured region (%d entries)'],nTarget,nBoth);
                end
                aH = aH(1:nTarget);  hxH = hxH(1:nTarget);  hyH = hyH(1:nTarget);
            end
        end
        out.a_hole  = aH(:)*1e-9;
        out.hx_hole = hxH(:)*1e-9;
        out.hy_hole = hyH(:)*1e-9;
        out.nholes  = numel(hxH);
        if ~isempty(tp) && strcmp(opts.EndTaper,'wvgmir')
            out.wvgmir    = numel(tp.hx);
            out.a_wvgmir  = tp.a(:)*1e-9;
            out.hx_wvgmir = tp.hx(:)*1e-9;
            out.hy_wvgmir = tp.hy(:)*1e-9;
        else
            out.wvgmir = 0;  out.a_wvgmir = [];  out.hx_wvgmir = [];  out.hy_wvgmir = [];
        end
        % ndef = last entry of the two-sided region still outside MirrorTol
        if isempty(mirror)
            out.ndef = nBoth;
        else
            mm = mirror*1e9;
            nC = min(numel(hxH),nBoth);
            dev = max(abs([aH(1:nC)'/mm(1), hxH(1:nC)'/mm(2), hyH(1:nC)'/mm(3)] - 1),...
                [],2,'omitnan');
            last = find(dev > opts.MirrorTol,1,'last');
            if isempty(last), last = 0; end
            out.ndef = last;
        end
    end
end

% =========================================================================
function x = holePositions(a,conv)
% hole positions [nm] from the a column
if strcmp(conv,'next')
    x = [0, cumsum(a(1:end-1))];
else
    x = [0, cumsum(a(2:end))];
end
end

function m = firstPh(ph)
% number of leading non-placeholder entries
m = find(ph,1) - 1;
if isempty(m), m = numel(ph); end
end

function [iL,iR] = pairRows(cRow,N,isPh)
% row pairs mirrored about cRow (integer or half-integer), measured rows only,
% contiguous from the centre
j = (1 - 0.5*(mod(cRow,1) ~= 0)):1:min(cRow-1,N-cRow);
iL = cRow - j;  iR = cRow + j;
last = firstPh(isPh(iL) | isPh(iR));
iL = iL(1:last);  iR = iR(1:last);
end

function printSummary(S,csvFile,measL,measR)
[~,nm,ext] = fileparts(csvFile);
fprintf('\nreadSEMHoles: %s%s\n',nm,ext);
fprintf('  centre search, best first (rms left-right mismatch of hx, hy):\n');
cs = S.centreScan(1:min(4,height(S.centreScan)),:);
for k = 1:height(cs)
    fprintf('    centre %5.1f  %2d pairs  hx %6.2f nm  hy %6.2f nm  score %.4f\n',...
        cs.centre(k),cs.nPairs(k),cs.rmsHx_nm(k),cs.rmsHy_nm(k),cs.score(k));
end
fprintf('  -> centre = hole %g, holeatctr = %d\n',S.centerHole,S.holeatctr);
fprintf(['  a_nm read as ''%s''; spacing mismatch rms: ''next'' %.2f nm, ',...
    '''prev'' %.2f nm\n'],S.aConvention,S.aSymmetry.rmsSpacingMismatch_nm);
if ~isempty(S.placeholder)
    iPh = S.file.hole == S.placeholder(1);
    fprintf('  placeholder rows: holes %d-%d (a/hx/hy = %s nm)\n',S.placeholder(1),...
        S.placeholder(end),mat2str([S.file.a_nm(iPh),S.file.hx_nm(iPh),S.file.hy_nm(iPh)]));
end
if ~isempty(S.mirror)
    fprintf('  padding mirror cell a/hx/hy = %s nm\n',mat2str(round(S.mirror*1e9,3)));
end
fprintf('  measured entries from the centre: left %d, right %d\n',measL,measR);
p = S.pairs;
fprintf('    j   L   R |  hxL   hxR |  hyL   hyR |  sL    sR    [nm]\n');
for k = 1:height(p)
    fprintf('  %3d %3d %3d | %5.1f %5.1f | %5.1f %5.1f | %5.1f %5.1f\n',p.j(k),...
        p.holeL(k),p.holeR(k),p.hxL_nm(k),p.hxR_nm(k),p.hyL_nm(k),p.hyR_nm(k),...
        p.spacingL_nm(k),p.spacingR_nm(k));
end
d = [p.hxL_nm - p.hxR_nm, p.hyL_nm - p.hyR_nm, p.spacingL_nm - p.spacingR_nm];
fprintf('  L-R mismatch rms/max: hx %.2f/%.2f  hy %.2f/%.2f  a %.2f/%.2f nm\n',...
    rms(d(:,1)),max(abs(d(:,1))),rms(d(:,2)),max(abs(d(:,2))),...
    rms(d(:,3),'omitnan'),max(abs(d(:,3))));
fprintf('  half ''%s'': nholes %d, ndef %d, wvgmir %d (EndTaper ''%s'')\n',...
    S.half,S.nholes,S.ndef,S.wvgmir,S.endTaper);
if ~isempty(S.L)
    fprintf('  left half: nholes %d, ndef %d, wvgmir %d\n',S.L.nholes,S.L.ndef,S.L.wvgmir);
end
end
