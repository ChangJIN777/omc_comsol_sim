function FP = modeFingerprints(U, G, divU)
%MODEFINGERPRINTS  Cheap scalar descriptors of sampled modes (pure MATLAB).
%
%   FP = modeFingerprints(U, G)
%   FP = modeFingerprints(U, G, divU)
%
%   U is 3*nPts x nb (one k-point, rows [ux; uy; uz]), G the ribModeGrid,
%   divU the optional nPts x nb volumetric strain. All outputs are 1 x nb
%   (or k x nb) and independent of the amplitude and phase of each mode:
%     FP.pol   3 x nb  fractions of |u|^2 in ux, uy, uz (sum 1)
%     FP.reg   3 x nb  fractions of |u|^2 in core, rib, arms (kinetic-energy
%                      fractions; uniform density)
%     FP.px, FP.py, FP.pz   mirror parities <u|M u>/<u|u> in [-1, 1]. px is
%                      meaningful at Gamma and X only (the x mirror maps k to
%                      -k). py / pz are NaN when only half the cell is built:
%                      the boundary condition fixes that parity.
%     FP.tip   coherent outward rib-tip displacement / RMS amplitude, mapped
%              to [0,1) by x/(1+x). A moving-boundary proxy for "breathing".
%     FP.dil   |mean volumetric strain| / RMS volumetric strain, in [0,1]. A
%              photoelastic proxy; NaN without divU.
%     FP.vec   11 x nb, rows FP.names. Proxies only: none of them is g_OM.

if nargin < 3, divU = []; end
U = double(U);
nPts = numel(G.w);  nb = size(U,2);
ok = all(reshape(all(isfinite(U),2), nPts, 3), 2);     % usable in every band
w  = G.w(:).*ok;
U(~repmat(ok,3,1), :) = 0;
ix = 1:nPts;  iy = nPts + ix;  iz = 2*nPts + ix;
E  = abs(U).^2;
Ep = E(ix,:) + E(iy,:) + E(iz,:);
Et = max(w'*Ep, realmin);

FP.pol = [w'*E(ix,:); w'*E(iy,:); w'*E(iz,:)]./Et;
FP.reg = zeros(3, nb);
for r = 1:3
    FP.reg(r,:) = ((w.*(G.region(:) == r))'*Ep)./Et;
end
FP.px = parity(U, w, ok, G.mirrorX, [-1  1  1], nPts, Et);
FP.py = parity(U, w, ok, G.mirrorY, [ 1 -1  1], nPts, Et);
FP.pz = parity(U, w, ok, G.mirrorZ, [ 1  1 -1], nPts, Et);

qy = G.q(2,:)';
tipPts = G.region(:) == 2 & abs(abs(qy) - max(abs(qy))) < 1e-9;
wt  = w.*tipPts;
tip = abs((wt.*sign(qy))'*U(iy,:))./max(sum(wt), realmin) ...
      ./sqrt(Et/max(sum(w), realmin));
FP.tip = tip./(1 + tip);

FP.dil = nan(1, nb);
if ~isempty(divU)
    D = double(divU);
    D(~ok,:) = 0;  D(~isfinite(D)) = 0;
    FP.dil = abs(w'*D)./sqrt(max(sum(w)*(w'*abs(D).^2), realmin));
end

FP.names = {'polx','poly','polz','regCore','regRib','regArm', ...
            'px','py','pz','tip','dil'};
FP.vec = [FP.pol; FP.reg; FP.px; FP.py; FP.pz; FP.tip; FP.dil];
end

% -------------------------------------------------------------------------

function p = parity(U, w, ok, perm, sg, nPts, Et)
% <u|M u>/<u|u> for the mirror with point permutation perm and component
% signs sg. NaN when the mirror image is not part of the built cell.
p = nan(1, size(U,2));
if isempty(perm), return; end
perm = perm(:);
wp = w.*ok(perm);
num = zeros(1, size(U,2));
for c = 1:3
    Uc = U((c-1)*nPts + (1:nPts), :);
    num = num + sg(c)*real(sum(wp.*conj(Uc).*Uc(perm,:), 1));
end
p = num./Et;
end
