function [M, nValid] = bandOverlap(UA, UB, w, mode)
%BANDOVERLAP  Modal assurance criterion between two sets of sampled modes.
%
%   M = bandOverlap(UA, UB)             unweighted
%   M = bandOverlap(UA, UB, w)          w: nPts weights, or {wA, wB} (two
%                                       designs; their geometric mean is used)
%   S = bandOverlap(UA, UB, w, 'subspace')
%   [M, nValid] = ...
%
%   UA is 3*nPts x nA and UB is 3*nPts x nB, rows [ux; uy; uz] on the same
%   ribModeGrid (same point counts).
%     M(i,j) = |a_i^H W b_j|^2 / ((a_i^H W a_i)(b_j^H W b_j))   in [0, 1]
%   It is invariant to the sign and complex phase of each mode. A point where
%   any component of any mode is not finite (outside the solid in either
%   design) is left out; nValid is the number of points used.
%
%   'subspace' returns the scalar ||QA^H W QB||_F^2 / min(nA, nB), the
%   overlap of the two spans. Use it for degenerate groups, whose individual
%   members are only defined up to a rotation.

UA = double(UA);  UB = double(UB);
if size(UA,1) ~= size(UB,1) || mod(size(UA,1), 3) ~= 0
    error('bandOverlap:sizeMismatch', ...
        'UA and UB must both have 3*nPts rows on the same grid (%d vs %d).', ...
        size(UA,1), size(UB,1));
end
nPts = size(UA,1)/3;
if nargin < 3 || isempty(w), w = ones(nPts,1); end
if iscell(w), w = sqrt(w{1}(:).*w{2}(:)); end
ok = all(isfinite(UA),2) & all(isfinite(UB),2);
ok = repmat(all(reshape(ok, nPts, 3), 2), 3, 1);
nValid = nnz(ok)/3;
sw = sqrt(repmat(w(:), 3, 1));
A = UA(ok,:).*sw(ok);
B = UB(ok,:).*sw(ok);
if nargin > 3 && strcmpi(mode, 'subspace')
    [QA, ~] = qr(A, 0);  [QB, ~] = qr(B, 0);
    M = norm(QA'*QB, 'fro')^2/min(size(A,2), size(B,2));
    return
end
M = abs(A'*B).^2 ./ max(sum(abs(A).^2,1)'*sum(abs(B).^2,1), realmin);
end
