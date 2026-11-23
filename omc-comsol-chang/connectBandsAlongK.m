function S = connectBandsAlongK(F, U, w, anchorRow, opts)
%CONNECTBANDSALONGK  Align and connect band levels across k-points.
%
%   S = connectBandsAlongK(F)                      frequency continuity only
%   S = connectBandsAlongK(F, U, w)                mode-shape overlap
%   S = connectBandsAlongK(F, U, w, anchorRow, opts)
%
%   F   nk x nb frequencies as runBands returns them. The eigensolver keeps
%       the nb eigenvalues nearest its shift, so column j need not be the
%       same absolute band at every k: the window can slide.
%   U   3*nPts x nb x nk sampled modes at EVERY row of F (columns as F), or
%       [] to fall back to frequency continuity.
%   w   grid weights (ribModeGrid's G.w)
%   anchorRow   row whose columns define the labels (default nk = X)
%   opts.maxOff window slide searched in the frequency fallback (default 2)
%   opts.macMin smallest overlap accepted as a connection (default 0.5)
%
%   S.Fs, S.order   F sorted ascending along each row, and the permutation
%   S.off       nk x 1 window offset: sorted column j at row r is absolute
%               level j + S.off(r) (the lowest level seen anywhere is 1)
%   S.Flev      nk x nLev levels aligned by S.off, NaN where a level is
%               outside the solver window at that k
%   S.complete  1 x nLev, true if the level was found at every k
%   S.levAnchor 1 x nb level of each sorted column at the anchor row
%   S.branch    nk x nb id of the mode-shape branch each sorted column
%               belongs to (ids = sorted anchor columns; new ids for bands
%               with no connection). Follows the CHARACTER through
%               crossings. [] without U.
%   S.conf      nk x 1 smallest connected overlap from the previous row (NaN
%               at the anchor and without U)
%
%   Gaps are a property of the aligned LEVELS (S.Flev), not of the branches.

if nargin < 2, U = []; end
if nargin < 3, w = []; end
[nk, nb] = size(F);
if nargin < 4 || isempty(anchorRow), anchorRow = nk; end
if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts,'maxOff'), opts.maxOff = 2; end
if ~isfield(opts,'macMin'), opts.macMin = 0.5; end

[Fs, order] = sort(F, 2);
hasU = ~isempty(U);
off  = zeros(nk, 1);
conf = nan(nk, 1);
branch = nan(nk, nb);
branch(anchorRow,:) = 1:nb;
nextId = nb + 1;

% [to; from] row pairs: from the anchor down to row 1, then up to row nk
steps = [anchorRow-1:-1:1, anchorRow+1:nk
         anchorRow:-1:2,   anchorRow:nk-1];
for c = 1:size(steps, 2)
    r1 = steps(1,c);  r0 = steps(2,c);
    o = [];
    if hasU
        M = bandOverlap(U(:,order(r0,:),r0), U(:,order(r1,:),r1), w);
        % matchpairs pairs (i,j) only when cost < 2*costUnmatched
        pairs = matchpairs(1 - M, (1 - opts.macMin)/2);
        if ~isempty(pairs)
            conf(r1) = min(M(sub2ind(size(M), pairs(:,1), pairs(:,2))));
            o = round(median(pairs(:,2) - pairs(:,1)));
            branch(r1, pairs(:,2)) = branch(r0, pairs(:,1));
        end
        isNew = isnan(branch(r1,:));
        branch(r1, isNew) = nextId + (0:nnz(isNew)-1);
        nextId = nextId + nnz(isNew);
    end
    if isempty(o)
        o = freqOffset(Fs(r0,:), Fs(r1,:), opts.maxOff);
    end
    % column j at r0 continues as column j + o at r1
    off(r1) = off(r0) - o;
end
off = off - min(off);

Flev = nan(nk, nb + max(off));
for r = 1:nk
    Flev(r, off(r) + (1:nb)) = Fs(r,:);
end
if ~hasU, branch = []; end
S = struct('Fs',Fs, 'order',order, 'off',off, 'Flev',Flev, ...
    'complete',all(isfinite(Flev), 1), 'levAnchor',off(anchorRow) + (1:nb), ...
    'branch',branch, 'conf',conf, 'anchorRow',anchorRow);
end

% -------------------------------------------------------------------------

function o = freqOffset(fa, fb, maxOff)
% Window slide o that makes fb(j+o) the smoothest continuation of fa(j).
% Zero is tried first and wins ties.
nb = numel(fa);
best = Inf;  o = 0;
for cand = [0, reshape([1:maxOff; -(1:maxOff)], 1, [])]
    ja = max(1, 1-cand):min(nb, nb-cand);
    if isempty(ja), continue; end
    d = mean(abs(fb(ja + cand) - fa(ja)));
    if d < best*(1 - 1e-9)
        best = d;  o = cand;
    end
end
end
