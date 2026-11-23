function ds = surrogateRibBandsTracked(P, opt)
%SURROGATERIBBANDSTRACKED  Synthetic rib bands WITH mode shapes (no COMSOL).
%
%   ds = surrogateRibBandsTracked(P)
%   ds = surrogateRibBandsTracked(P, opt)
%
%   NOT physics. A stand-in for solveBands(P) with P.exportModes, used by
%   trackRibBands' dry run and by BandTrackingTest. It returns ds.sym (F,
%   k_norm, modes, midGap, gapSize, P, truth, detune), ds.full and ds.P in
%   the real format.
%
%   Each band is a fixed "character" c with its own mode shape:
%     1  valence band (x-even), top at X at fc - gap/2
%     2  conduction band (x-odd), bottom at X at fc + gap/2
%     3  INTRUDER (x-even): nearly flat; its frequency falls linearly with
%        the spine width s, from 1.45*fc at s = 50 nm to 0.85*fc at 200 nm.
%        With the default t, L, a it CROSSES the conduction band at X near
%        s = 140.5 nm (opposite parity: exact crossing) and meets the valence
%        band near s = 174.3 nm, where the two ANTICROSS and hybridize
%        (same parity, coupling opt.coupling).
%     4.. filler bands below and above (alternating parity)
%   Rows of ds.sym.F are sorted ascending, as COMSOL returns them.
%   ds.sym.truth(k, j) is the (dominant) character in column j at row k, and
%   ds.sym.detune(k) the uncoupled intruder-minus-valence frequency (Hz).
%
%   opt.n         ribModeGrid point counts
%   opt.coupling  coupling (Hz) between characters 1 and 3. Default 0.1e9.
%                 0: they cross exactly.
%
%   Mode shapes are polynomials in the grid's normalized coordinates with a
%   definite x parity, orthonormal, times an arbitrary complex phase per
%   mode and design. Characters 1 and 3 are the SAME scalar profile in uy
%   and in uz respectively, so the overlap between two hybrids of them is
%   exactly cos^2 of the change of the mixing angle, for any design.

if nargin < 2 || isempty(opt), opt = struct(); end
if ~isfield(opt,'n'),        opt.n = struct(); end
if ~isfield(opt,'coupling'), opt.coupling = 0.1e9; end

P = ribDefaultParams(P);
if ~isfield(P,'mbeveny'), P.mbeveny = 1; end
if ~isfield(P,'mbevenz'), P.mbevenz = 1; end
if ~isfield(P,'kpts'),    P.kpts = 10; end
if ~isfield(P,'nbands'),  P.nbands = 10; end
nb = P.nbands;  nk = P.kpts + 1;
if nb < 4
    error('surrogateRibBandsTracked:tooFewBands', 'P.nbands must be at least 4.');
end
G = ribModeGrid(P, opt.n);
nPts = size(G.q, 2);
L  = ribLength(P, P.a);
kn = (0:P.kpts)'/P.kpts;

%% frequencies per character
fc  = 11.2e9*(300e-9/P.a)*sqrt(370e-9/L)*(P.t/100e-9)^0.2;
gap = 1.6e9*exp(-((P.t-130e-9)/120e-9)^2 - ((L-400e-9)/240e-9)^2);
fV = fc - gap/2;  fC = fc + gap/2;
nLow = floor((nb - 3)/2);  nHigh = nb - 3 - nLow;
f = zeros(nk, nb);
f(:,1) = fV*(0.55 + 0.45*sin(pi*kn/2));
f(:,2) = fC*(1 + 0.30*cos(pi*kn/2).^2);
f(:,3) = fc*(1.45 - 0.6*(P.s - 50e-9)/150e-9)*(1 + 0.02*cos(pi*kn));
for j = 1:nLow
    f(:,3+j) = fV*0.5*(1 - j/(nLow + 1))*(1 + 0.3*kn);
end
for j = 1:nHigh
    f(:,3+nLow+j) = fC*(1.45 + 0.15*j)*(1 + 0.05*kn);
end
detune = f(:,3) - f(:,1);

%% design-independent shapes with definite x parity
qx = G.q(1,:)';  qy = G.q(2,:)';  qz = G.q(3,:)';
Be = [ones(nPts,1), qy, qy.^2, qx.^2, qz, qx.^2.*qy];      % even in x
Bo = [qx, qx.*qy, qx.^3, qx.*qz, qx.*qy.^2, qx.^3.*qy];    % odd in x
par = ones(1, nb);  par(2:2:end) = -1;
Phi = zeros(3*nPts, nb);
for c = 1:nb
    rs = RandStream('mt19937ar', 'Seed', 100 + c);
    cf = rs.randn(6, 3);
    if par(c) > 0          % Mx u = +u: ux odd, uy and uz even
        Phi(:,c) = [Bo*cf(:,1); Be*cf(:,2); Be*cf(:,3)];
    else
        Phi(:,c) = [Be*cf(:,1); Bo*cf(:,2); Bo*cf(:,3)];
    end
end
% Characters 1 and 3: one even scalar profile, in uy and in uz. They are
% orthogonal at every grid point and have equal norms under any weights.
g13 = Be*[1; 0.6; -0.8; 0.5; 0.3; 0.4];
Phi(:,1) = [zeros(nPts,1); g13; zeros(nPts,1)];
Phi(:,3) = [zeros(2*nPts,1); g13];
[Phi, ~] = qr(Phi, 0);      % columns 1-3 are already mutually orthogonal

%% assemble, sorted by frequency at each k
U = complex(zeros(3*nPts, nb, nk, 'single'));
F = zeros(nk, nb);  truth = zeros(nk, nb);
for m = 1:nk
    fk = f(m,:);  Uk = Phi;  lab = 1:nb;
    if opt.coupling > 0
        [V, D] = eig([fk(1) opt.coupling; opt.coupling fk(3)]);
        fk([1 3]) = diag(D)';
        Uk(:,[1 3]) = Phi(:,[1 3])*V;
        [~, dom] = max(abs(V), [], 1);
        ids = [1 3];
        lab([1 3]) = ids(dom);
    end
    [F(m,:), ord] = sort(fk);
    ph = exp(1i*(1.7*(1:nb) + 0.9*m + 1e9*P.t + 7e8*P.s));
    U(:,:,m) = single(Uk(:,ord).*ph(ord));
    truth(m,:) = lab(ord);
end

sec = struct('F',F, 'k_norm',kn, 'kx_norm',kn, 'ky_norm',zeros(nk,1), ...
    'truth',truth, 'detune',detune);
sec.modes = struct('U',U, 'divU',[], 'grid',G, 'kIdx',1:nk, 'k_norm',kn', ...
    'F',F, 'a',P.a, 'sector',[P.mbeveny P.mbevenz], 'freqCheck',zeros(1,nk));
[sec.midGap, sec.gapSize] = findGaps(sec);
sec.P = P;
full.F = F;
[full.midGap, full.gapSize] = findGaps(full);
ds.sym = sec;  ds.full = full;  ds.P = P;
end
