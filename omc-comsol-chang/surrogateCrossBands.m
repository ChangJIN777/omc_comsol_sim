function ds = surrogateCrossBands(P, OPT)
%SURROGATECROSSBANDS  Analytic stand-in for a band structure. NOT PHYSICS.
%
% =========================================================================
%  THIS FUNCTION DOES NOT COMPUTE A BAND STRUCTURE. It evaluates a few
%  closed-form expressions chosen to LOOK like one. The frequencies are not
%  approximations of the cross cell's modes, they are not converged, and there
%  is no eigenproblem anywhere below. NOTHING THAT COMES OUT OF HERE MAY BE
%  USED TO JUDGE A DESIGN. Its only purpose is to make the optimization loop
%  run in milliseconds so the loop itself can be debugged.
% =========================================================================
%
% What it must nonetheless get right:
%
% 1. SHAPE. Returns the same ds a real 2-sector solveBands call does: .sym and
%    .asym each with F [nk x nbands] in Hz and k_norm [nk x 1] running 0->3
%    (the Gamma-X-M-Gamma circuit runBands_2D walks), and .full with F, midGap
%    and gapSize. Get this wrong and the dry run tests the surrogate instead
%    of the pipeline.
%
% 2. DETERMINISM. No rand, no clock, no counters. The result depends only on
%    (a,h,w,th) and on OPT.kpts/OPT.nbands. A surrogate that answered
%    differently on a repeated design would be exactly the class of bug a dry
%    run exists to find rather than introduce.
%
% 3. A NON-TRIVIAL INTERIOR OPTIMUM, so a convergence trace is a real test
%    rather than a slide into a bound. Four competing Gaussian factors:
%       arm length ratio  h/a   -> best near 0.75
%       arm width ratio   w/h   -> best near 0.32
%       slab aspect       th/a  -> best near 0.70
%       void filling      fill  -> best near 0.35
%    against one scaling law, f ~ 1/a, so the target-frequency penalty in
%    calFitness FIGHTS the gap-width factors instead of agreeing with them.
%    A floor is subtracted so part of the box has no gap at all and the
%    gapless branch of the objective is exercised too.
%
% 4. SELF-CONSISTENCY. midGap and gapSize are NOT written by hand: the
%    synthesized band matrix is passed through the REAL findGaps, exactly as
%    solveBands does. Hand-written values could disagree with the bands beside
%    them, and would leave findGaps itself untested.

nk     = 3*OPT.kpts + 1;          % Gamma-X-M-Gamma, plus the wrapped Gamma
k_norm = linspace(0, 3, nk).';
nb     = OPT.nbands;

% --- dimensionless design ratios ---
hOverA  = P.h / P.a;
wOverH  = P.w / P.h;
thOverA = P.th / P.a;

% Void area: two crossed bars of length h and width w, overlapping in a w x w
% square at the centre, over a square cell of side a.
fill = (2*P.h*P.w - P.w^2) / P.a^2;

% --- gap strength: product of Gaussians, one interior optimum each ---
gap = 0.35 ...
    * exp(-((hOverA  - 0.75)/0.18)^2) ...
    * exp(-((wOverH  - 0.32)/0.12)^2) ...
    * exp(-((thOverA - 0.70)/0.25)^2) ...
    * exp(-((fill    - 0.35)/0.15)^2);

% Floor, so a good fraction of the design box has NO complete gap.
gap = gap - 0.04;

% --- frequency scale: Bragg-like, f ~ v/(2a) ---
vEff   = 12000;                   % m/s, diamond shear-ish. Plausible, not fitted.
fBragg = vEff / (2*P.a);

% --- synthesize two symmetry sectors ---
% The sym sector carries the gap; the asym sector is offset so the COMPLETE
% gap (their union, via findGaps on .full) is narrower than the sym-only gap,
% as it would be in a real calculation.
% The asym sector is offset slightly and given a slightly narrower gap, so
% the COMPLETE gap (their union, via findGaps on .full) is narrower than the
% sym-only gap, as it would be in a real calculation.
%
% These two numbers are load-bearing and were WRONG on the first attempt
% (offset 0.05, reduction 0.06). A complete gap only survives where the two
% sector gaps OVERLAP, which needs roughly gap > offset + reduction/2; at
% those values that meant gap > 0.08, and almost nothing in the design box
% cleared it. A 40-evaluation bayesopt study then returned "no_gap" for every
% single point -- a landscape with no signal at all to optimize. Keep the sum
% of these two small relative to the attainable gap, and re-check the
% complete-gap fraction if you change them.
ds.sym  = makeSector(k_norm, nb, fBragg, max(gap, 0),        0.00);
ds.asym = makeSector(k_norm, nb, fBragg, max(gap - 0.02, 0), 0.02);

[ds.sym.midGap,  ds.sym.gapSize ] = findGaps(ds.sym);
[ds.asym.midGap, ds.asym.gapSize] = findGaps(ds.asym);

ds.full.F = [ds.sym.F, ds.asym.F];
[ds.full.midGap, ds.full.gapSize] = findGaps(ds.full);
ds.full.k_norm = k_norm;

% PROVENANCE: the in-file marker the objective checks on every evaluation.
ds.isSynthetic  = true;
ds.surrogateTag = sprintf(['SYNTHETIC band data from surrogateCrossBands - ', ...
    'NOT physics. a=%.1fnm h=%.1fnm w=%.1fnm th=%.1fnm'], ...
    P.a*1e9, P.h*1e9, P.w*1e9, P.th*1e9);
end

%% ------------------------------------------------------------------------
function sec = makeSector(k_norm, nb, fBragg, gapFrac, offset)
%MAKESECTOR  One synthetic symmetry sector with a controllable gap.
%
%   Bands are split into a group below fLow and a group above fHigh, so a gap
%   of width fHigh-fLow exists by construction and findGaps must find it. With
%   gapFrac = 0 the two groups meet and there is no gap -- which is how the
%   gapless branch gets exercised.
fCentre = fBragg * (1 + offset);
fLow    = fCentre * (1 - gapFrac/2);
fHigh   = fCentre * (1 + gapFrac/2);

nLow = max(1, floor(nb/2));
nUp  = nb - nLow;

% Dispersion ramp, 0 at Gamma rising to 1 at the end of the circuit.
ramp = (1 - cos(pi*k_norm/3))/2;

% CONSECUTIVE BANDS MUST OVERLAP IN FREQUENCY, or findGaps discovers spurious
% gaps BETWEEN them and the objective's max(gapSize) can pick one of those
% instead of the gap this function is trying to place. (First version here
% used non-overlapping bands and reported a 4.96 GHz mid-gap for a surrogate
% designed around 37.5 GHz.) Each band spans 1.6 spacings while consecutive
% starts are 1 spacing apart, so neighbours always overlap by 0.6.
span = 1.6;

F = zeros(numel(k_norm), nb);

% Below the gap: densely fill [0, fLow]; the top band reaches exactly fLow.
for m = 1:nLow
    F(:, m) = fLow * ((m-1) + span*ramp) / (nLow - 1 + span);
end

% Above the gap: start exactly at fHigh and fan upward.
spread = 0.8;
for j = 1:nUp
    F(:, nLow + j) = fHigh * ...
        (1 + spread * ((j-1) + span*ramp) / (nUp - 1 + span));
end

sec.F      = F;
sec.k_norm = k_norm;
end
