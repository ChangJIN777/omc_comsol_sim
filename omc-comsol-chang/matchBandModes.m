function R = matchBandModes(ref, tIdx, cand, opts)
%MATCHBANDMODES  Find chosen modes of a reference design in a new design.
%
%   R = matchBandModes(ref, tIdx, cand)
%   R = matchBandModes(ref, tIdx, cand, opts)
%
%   ref, cand   bandTrackCandidate outputs (X-point modes sorted by
%               frequency) of the reference design and of the new design
%   tIdx        sorted columns of ref of the modes to follow (1 x nT)
%
%   METHOD. ALL reference modes are assigned to the candidate's modes at once
%   (matchpairs), so no candidate is claimed twice, with the cost
%     C = wO*(1 - MAC) + wF*min(1, |f - fRef|/(fTol*fRef))
%         + wP*mean|feature difference| + gateCost*[x parity flipped]
%   Each followed mode is then classified:
%     'confident'  assigned, MAC >= macHi, MAC - runner-up >= marginMin and
%                  |f - fRef|/fRef <= fTol
%     'lost'       no admissible candidate resembles it (best MAC < macLo;
%                  candidates with flipped x parity are not admissible)
%     'ambiguous'  everything in between: typically a hybridized mode near an
%                  avoided crossing, or too large a design step. idx is then
%                  the best guess.
%   Candidates within degTol of the match form a degenerate group; if the
%   single-mode MAC is low but the reference lies in the group's span, the
%   group overlap is used and R.f is the group mean.
%
%   opts (defaults)  macHi 0.80, macLo 0.50, marginMin 0.30, fTol 0.08,
%       degTol 2e-3, wO 0.6, wF 0.2, wP 0.2, costUnmatched 0.4,
%       useParity true (x-mirror parity veto at X; valid when the crystal
%       orientation keeps the x mirror, rxtal = 0 or 45; set false otherwise),
%       parityMin 0.8, gateCost 1e3
%
%   R (fields are 1 x nT)  idx (sorted candidate column, NaN if lost), f,
%       mac, mac2 (runner-up), margin, groupMac, degenerate, fErr, dFeat,
%       cost, status (cellstr); and R.ok (all confident), R.anyAmbiguous,
%       R.M (the MAC matrix), R.nValid (grid points used).

if nargin < 4 || isempty(opts), opts = struct(); end
def = struct('macHi',0.80, 'macLo',0.50, 'marginMin',0.30, 'fTol',0.08, ...
    'degTol',2e-3, 'wO',0.6, 'wF',0.2, 'wP',0.2, 'costUnmatched',0.4, ...
    'useParity',true, 'parityMin',0.8, 'gateCost',1e3);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts, fn{k}) || isempty(opts.(fn{k})), opts.(fn{k}) = def.(fn{k}); end
end
tIdx = tIdx(:)';  nT = numel(tIdx);
if any(~isfinite(tIdx)) || any(tIdx < 1 | tIdx > numel(ref.f) | tIdx ~= round(tIdx))
    error('matchBandModes:badTargets', ...
        'tIdx must hold sorted column indices of the reference (1..%d).', numel(ref.f));
end

R = struct('idx',nan(1,nT), 'f',nan(1,nT), 'mac',nan(1,nT), 'mac2',nan(1,nT), ...
    'margin',nan(1,nT), 'groupMac',nan(1,nT), 'degenerate',false(1,nT), ...
    'fErr',nan(1,nT), 'dFeat',nan(1,nT), 'cost',nan(1,nT));
R.status = repmat({'lost'}, 1, nT);

%% cost matrix: rows = reference modes, columns = candidate modes
wPair = {ref.w, cand.w};
[M, R.nValid] = bandOverlap(ref.U, cand.U, wPair);
fRef = ref.f(:);
dF = min(1, abs(cand.f(:)' - fRef)./(opts.fTol*fRef));
dP = mean(abs(permute(ref.fp,[2 3 1]) - permute(cand.fp,[3 2 1])), 3, 'omitnan');
dP(isnan(dP)) = 0;
gate = false(size(M));
if opts.useParity
    pr = ref.px(:);  pc = cand.px(:)';
    gate = (abs(pr) > opts.parityMin) & (abs(pc) > opts.parityMin) & ...
           (sign(pr) ~= sign(pc));
end
C = opts.wO*(1 - M) + opts.wF*dF + opts.wP*dP + opts.gateCost*gate;

%% global assignment, then classify each followed mode
pairs = matchpairs(C, opts.costUnmatched);
for q = 1:nT
    t = tIdx(q);
    row = M(t,:);
    rowOK = row;  rowOK(gate(t,:)) = 0;         % parity-flipped: not admissible
    [macBest, jBest] = max(rowOK);
    j = pairs(pairs(:,1) == t, 2);
    paired = ~isempty(j);
    if ~paired, j = jBest; end

    R.idx(q) = j;  R.f(q) = cand.f(j);  R.mac(q) = rowOK(j);
    R.cost(q) = C(t,j);  R.dFeat(q) = dP(t,j);
    macEff = R.mac(q);  excl = j;
    grp = find(abs(cand.f - cand.f(j)) <= opts.degTol*cand.f(j));
    if numel(grp) > 1 && ~any(gate(t,grp))
        R.groupMac(q) = bandOverlap(ref.U(:,t), cand.U(:,grp), wPair, 'subspace');
        if R.mac(q) < opts.macHi && R.groupMac(q) >= opts.macHi
            R.degenerate(q) = true;
            macEff = R.groupMac(q);  excl = grp;  R.f(q) = mean(cand.f(grp));
        end
    end
    others = rowOK;  others(excl) = [];
    R.mac2(q) = max([others, 0]);
    R.margin(q) = macEff - R.mac2(q);
    R.fErr(q) = abs(R.f(q) - fRef(t))/fRef(t);

    if max(macBest, macEff) < opts.macLo
        R.idx(q) = NaN;  R.f(q) = NaN;  R.mac(q) = macBest;   % 'lost'
    elseif paired && macEff >= opts.macHi && R.margin(q) >= opts.marginMin ...
            && R.fErr(q) <= opts.fTol
        R.status{q} = 'confident';
    else
        R.status{q} = 'ambiguous';
    end
end
R.ok = all(strcmp(R.status, 'confident'));
R.anyAmbiguous = any(strcmp(R.status, 'ambiguous'));
R.M = M;
end
