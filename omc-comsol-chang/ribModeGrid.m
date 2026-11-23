function G = ribModeGrid(P, n)
%RIBMODEGRID  Design-independent sampling grid for the 'rib' unit cell.
%
%   G = ribModeGrid(P)
%   G = ribModeGrid(P, n)
%
%   The cell is split into three regions, each sampled at cell midpoints in
%   normalized coordinates, so point p is the SAME material location in every
%   design and mode shapes of different designs can be compared index-wise:
%     1 core  |x| <= t/2, |y| <= s/2        x = t/2*xi,  y = s/2*etaS
%     2 rib   |x| <= t/2, s/2 < |y| <= L    x = t/2*xi,  y = s/2 + (L-s/2)*eta
%     3 arms  t/2 < |x| <= a/2, |y| <= s/2  x = +-(t/2 + (a/2-t/2)*alpha)
%   and z = th/2*zeta. Only the part of the cell the builder builds is
%   sampled: y >= 0 when P.mbeveny ~= 0, z >= 0 when P.mbevenz ~= 0 (the rule
%   of buildRibUnitCell_LN for mechanical runs).
%
%   n   point counts (all optional): xi (across t, default 4), eta (along the
%       rib, 6), etaS (across the spine half-width, 2), arm (along each arm,
%       3), zeta (per z half, 2). Compared designs must use the same n.
%
%   G.xyz      3 x nPts coordinates (m) for mphinterp
%   G.q        3 x nPts normalized coordinates (design independent):
%              qx = xi in (-1,1), or +-(1+alpha) on the arms; qy = +-etaS, or
%              +-(1+eta) on the rib; qz = +-zeta
%   G.region   1 x nPts region index (1 core, 2 rib, 3 arms)
%   G.w        1 x nPts volume weights, sum 1 (uniform density: mass weights)
%   G.mirrorX  1 x nPts index of the point at (-x, y, z)
%   G.mirrorY, G.mirrorZ   same for y and z; [] when only half is built
%   G.symY, G.symZ, G.n, G.dims
%
%   Points near the rib-tip fillets (r2) or, with a sidewall angle, near the
%   leaning walls can lie outside the solid. mphinterp returns NaN there and
%   bandOverlap / modeFingerprints skip those points.

if nargin < 2 || isempty(n), n = struct(); end
def = struct('xi',4, 'eta',6, 'etaS',2, 'arm',3, 'zeta',2);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(n, fn{k}), n.(fn{k}) = def.(fn{k}); end
end

P = ribDefaultParams(P);
a = P.a;  t = P.t;  s = P.s;  th = P.th;
L = ribLength(P, a);
symY = isfield(P,'mbeveny') && P.mbeveny ~= 0;
symZ = isfield(P,'mbevenz') && P.mbevenz ~= 0;

mid  = @(m) ((1:m) - 0.5)/m;            % cell midpoints in (0,1)
xi   = 2*mid(n.xi) - 1;                 % (-1,1) across the rib width
eta  = mid(n.eta);  etaS = mid(n.etaS);  al = mid(n.arm);  ze = mid(n.zeta);

% In-plane points for y > 0. Columns: [qx qy x y region areaWeight]
[XI, ES] = ndgrid(xi, etaS);                                    % 1 core
m = numel(XI);
core = [XI(:), ES(:), t/2*XI(:), s/2*ES(:), ones(m,1), repmat(t*(s/2)/m, m, 1)];
[XI, ER] = ndgrid(xi, eta);                                     % 2 rib
m = numel(XI);
rib = [XI(:), 1+ER(:), t/2*XI(:), s/2+(L-s/2)*ER(:), 2*ones(m,1), ...
       repmat(t*(L-s/2)/m, m, 1)];
[AL, ES] = ndgrid(al, etaS);                                    % 3 arms
m = numel(AL);
xa = t/2 + (a/2 - t/2)*AL(:);
wa = repmat((a/2 - t/2)*(s/2)/m, m, 1);
arm = [ 1+AL(:), ES(:),  xa, s/2*ES(:), 3*ones(m,1), wa
       -1-AL(:), ES(:), -xa, s/2*ES(:), 3*ones(m,1), wa];
pl = [core; rib; arm];
if ~symY
    plm = pl;  plm(:,[2 4]) = -plm(:,[2 4]);
    pl = [pl; plm];
end
zq = ze;
if ~symZ, zq = [-fliplr(ze), ze]; end

np = size(pl,1);  nz = numel(zq);
G.q      = [repmat(pl(:,1:2)', 1, nz); kron(zq, ones(1,np))];
G.xyz    = [repmat(pl(:,3:4)', 1, nz); kron(th/2*zq, ones(1,np))];
G.region = repmat(pl(:,5)', 1, nz);
w = repmat(pl(:,6)', 1, nz);
G.w = w/sum(w);
G.mirrorX = mirrorPerm(G.q, 1);
G.mirrorY = [];
G.mirrorZ = [];
if ~symY, G.mirrorY = mirrorPerm(G.q, 2); end
if ~symZ, G.mirrorZ = mirrorPerm(G.q, 3); end
G.symY = symY;  G.symZ = symZ;  G.n = n;
G.dims = struct('a',a, 't',t, 's',s, 'L',L, 'th',th);
end

% -------------------------------------------------------------------------

function p = mirrorPerm(q, ax)
% Index of each point's mirror image across the plane normal to axis ax.
qm = q;  qm(ax,:) = -qm(ax,:);
[tf, p] = ismembertol(qm', q', 1e-9, 'ByRows', true);
if all(tf), p = p'; else, p = []; end
end
