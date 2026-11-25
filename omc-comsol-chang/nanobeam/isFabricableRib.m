function [ok, margin, why, info] = isFabricableRib(defectAspectRatio, maxdef, P)
%ISFABRICABLERIB  Build-validity and fabricability of a rib-cavity taper, no solve.
%
%   [ok, margin, why, info] = isFabricableRib(defectAspectRatio, maxdef, P)
%
%   Rib counterpart of isFabricable.m, for the hole-style taper of the rib
%   cavity (P.celltype = 'rib', P.ribTaperScheme = 'hole'). Pure arithmetic
%   on the taper parameters: no geometry build, no COMSOL. Vectorized over
%   candidates, because bayesopt's XConstraintFcn passes a table of many rows.
%
%   INPUTS
%     defectAspectRatio - rib aspect ratio w/t of the CENTRE cell (tip-to-tip
%                         width over rib width). Scalar or array.
%     maxdef            - defect ratio. Scalar or array (expanded against the
%                         other input).
%     P                 - rib parameter struct, as passed to
%                         createRibCavityGeom. P.maxdef and P.oblong are
%                         ignored (they are what the two inputs replace).
%
%   MAPPING (the same as optimize_rib_oblong_maxdef)
%     mirrorAspectRatio = w_mirror / t
%     oblong = log(defectAspectRatio/mirrorAspectRatio) / (2*log(1 - maxdef))
%     sig_k  = 1 - maxdef*f((k-1)/ndef)     (1 for mirror cells)
%     a_k = a*sig_k,  t_k = t*sig_k^(1-oblong),  w_k = w_mirror*sig_k^(1+oblong)
%
%   OUTPUTS
%     ok     - logical, size of the inputs: true where the design passes
%              every active rule.
%     margin - slack of the tightest active rule in METRES (negative =
%              violated, and by how much; -Inf for a non-physical mapping).
%     why    - cell array of reason codes naming the tightest rule where
%              ~ok, '' where ok.
%     info   - struct of arrays: .valid (BUILD_* rules only: what
%              createRibCavityGeom accepts), .warn (a WARN_* rule tripped),
%              .fab (FAB_* rules pass), .oblong, .a_ctr, .t_ctr, .w_ctr.
%
%   RULES. Three groups, kept apart on purpose.
%   BUILD_*  hard errors of validateRibParams.m / createRibCavityGeom.m.
%            Mirrored exactly; test_isFabricableRib cross-validates them
%            against the real builder.
%     BUILD_RIB_TOO_WIDE        t < a at both taper ends     (rib:ribTooWide)
%     BUILD_RIB_SHORTER_THAN_SPINE  min L > s/2     (rib:ribShorterThanSpine)
%     BUILD_TIP_FILLET          min t > 2*r2          (rib:tipFilletTooLarge)
%     BUILD_ROOT_FILLET         min (a - t)/2 > r1   (rib:rootFilletTooLarge)
%     BUILD_SEGMENT_TOO_NARROW  every rib inside its spine segment
%                                                    (rib:segmentTooNarrow)
%     BUILD_AIR_TOO_SMALL       air cylinder encloses the beam (rib:airTooSmall)
%     BUILD_SIDEWALL_TOO_STEEP  d > 0 only            (rib:sidewallTooSteep)
%   WARN_*   conditions the builder only WARNS about ("COMSOL may refuse the
%            fillet"). Rejected here unless P.ribFabStrict = 0.
%     WARN_FILLETS_MEET         min L - s/2 > r1 + r2       (rib:filletsMeet)
%     WARN_ROOT_FILLET_TIGHT    per-cell spine beside the rib > r1
%                                                      (rib:rootFilletTight)
%     WARN_STAGING_OVERLAP      4*min L + s/2 - 10 nm > max L
%                                                      (rib:stagingOverlap)
%   FAB_*    lithography floors. Defaults are the rib limits of
%            python-scripts/configs/bounds.yaml:81-83 (50 nm). 0 = rule off.
%     FAB_MIN_FEATURE  t_k, s and the rib arm (w_k - s)/2   P.minRibFeature
%     FAB_MIN_GAP      gap between neighbouring ribs a_k - t_k  P.minRibGap
%     FAB_MIN_THICKNESS  th                               P.minRibThickness
%
%   Example:
%     P = struct('celltype','rib','ribTaperScheme','hole','a',600e-9, ...
%         't',280e-9,'s',100e-9,'th',250e-9,'w_mirror',1300e-9,'airZOffset',0);
%     [ok, margin, why] = isFabricableRib(11.3, 0.2, P)
%
%   Needs ribDefaultParams.m and ribLength.m (omc-comsol-chang root) on the path.
%
%   See also ISFABRICABLE, CREATERIBCAVITYGEOM, VALIDATERIBPARAMS,
%   OPTIMIZE_RIB_OBLONG_MAXDEF, TEST_ISFABRICABLERIB.

narginchk(3, 3);
validateattributes(defectAspectRatio, {'numeric'}, {'real'}, mfilename, 'defectAspectRatio');
validateattributes(maxdef,            {'numeric'}, {'real'}, mfilename, 'maxdef');

if exist('ribDefaultParams','file') ~= 2 || exist('ribLength','file') ~= 2
    error('isFabricableRib:pathMissing', ...
        ['ribDefaultParams.m / ribLength.m (omc-comsol-chang root) are not on ' ...
         'the path. Add the root with addpath(<root>,''-end'').']);
end
Pv = ribDefaultParams(P);
if ~strcmp(Pv.ribTaperScheme, 'hole')
    error('isFabricableRib:needsHoleScheme', ...
        'isFabricableRib maps (aspect ratio, maxdef) to P.oblong, so it needs P.ribTaperScheme = ''hole''.');
end

% Thresholds. Empty or missing -> default.
minFeature   = getfielddef(P, 'minRibFeature',   50e-9);
minGap       = getfielddef(P, 'minRibGap',       50e-9);
minThickness = getfielddef(P, 'minRibThickness', 50e-9);
strictWarn   = getfielddef(P, 'ribFabStrict',    1);

[dAR, md] = expandPair(defectAspectRatio, maxdef);
outSize = size(dAR);
dAR = dAR(:);
md  = md(:);
n   = numel(dAR);

% Mirror-cell half-length through ribLength, so the P.w_mirror / L(a)-rule
% precedence is not duplicated here. (ribLength itself rejects P.w_ctr.)
[Lm, Lsrc] = ribLength(Pv, Pv.a);
if any(strcmp(Lsrc, {'L','w'}))
    error('isFabricableRib:holeSchemeNeedsMirrorWidth', ...
        ['An explicit P.L / P.w cannot be used with P.ribTaperScheme = ''hole''. ' ...
         'Give the mirror-cell width as P.w_mirror, or use the L(a) rule.']);
end
mirrorAspectRatio = 2*Lm / Pv.t;

% --- map the aspect ratio back to the taper exponent ---------------------
% Non-physical inputs become NaN here and are rejected below, so no complex
% or infinite oblong ever reaches the arithmetic.
physical = md > 0 & md < 1 & dAR > 0;
dARs = dAR;
dARs(~physical) = NaN;
sigC = 1 - md;
sigC(~physical) = NaN;
oblong = log(dARs ./ mirrorAspectRatio) ./ (2*log(sigC));
bad = ~isfinite(oblong);

% --- per-cell dimensions, n candidates x N cells (createRibCavityGeom) ----
N = Pv.ncellHalf;
k = 1:N;
h = double(Pv.ribatctr);
f = taperProfile(Pv.taperFunc);
inDef = k <= Pv.ndef;
fk = zeros(1, N);
fk(inDef) = f((k(inDef) - 1)/Pv.ndef);
sig  = 1 - md*fk;
aLoc = Pv.a*sig;
tLoc = Pv.t*sig.^(1 - oblong);
LLoc = Lm*sig.^(1 + oblong);

% validateRibParams checks the two ENDS of the taper, a*(1 - maxdef) and a,
% whether or not a cell actually sits there.
sigE = [sigC, ones(n, 1)];
aE = Pv.a*sigE;
tE = Pv.t*sigE.^(1 - oblong);
LE = Lm*sigE.^(1 + oblong);
LMin = min(LE, [], 2);
LMax = max(LE, [], 2);
tMin = min(tE, [], 2);
gapE = min(aE - tE, [], 2);

% --- cell positions and spine-segment edges (createRibCavityGeom) --------
switch Pv.taperConvention
    case 'centerspacing'
        fs  = zeros(1, N);
        idx = inDef & (k - 1 - h >= 0);
        fs(idx) = f((k(idx) - 1 - h)/(Pv.ndef - h));
        sp = Pv.a*(1 - md*fs);
        if h
            xc = [zeros(n, 1), cumsum(sp(:, 2:end), 2)];
        else
            xc = cumsum([sp(:, 1)/2, sp(:, 2:end)], 2);
        end
        mid = (xc(:, 1:end-1) + xc(:, 2:end))/2;
        if h
            if N > 1
                leftEdge = -mid(:, 1);
            else
                leftEdge = -aLoc(:, 1)/2;
            end
        else
            leftEdge = zeros(n, 1);
        end
        if h && N == 1
            lastSp = aLoc(:, end);
        else
            lastSp = sp(:, end);
        end
        xL = [leftEdge, mid];
        xR = [mid, xc(:, end) + lastSp/2];
    case 'cellwidth'
        if h
            x0 = zeros(n, 1);
        else
            x0 = aLoc(:, 1)/2;
        end
        xc = x0 + [zeros(n, 1), cumsum((aLoc(:, 1:end-1) + aLoc(:, 2:end))/2, 2)];
        xL = xc - aLoc/2;
        xR = xc + aLoc/2;
    otherwise
        error('isFabricableRib:badTaperConvention', ...
            'P.taperConvention = ''%s'' must be ''centerspacing'' or ''cellwidth''.', ...
            Pv.taperConvention);
end
clrMin = min(min(xc - xL, xR - xc) - tLoc/2, [], 2);

% --- the rules: {code, signed margin in metres, group, strict} ------------
% strict = true: the builder errors/warns AT equality, so the margin must be
% > 0. strict = false (fab floors): a feature exactly at the floor passes.
onesN = ones(n, 1);
rules = { ...
    'BUILD_RIB_TOO_WIDE',           gapE,                 'build', true; ...
    'BUILD_RIB_SHORTER_THAN_SPINE', LMin - Pv.s/2,        'build', true; ...
    'BUILD_SEGMENT_TOO_NARROW',     clrMin,               'build', true; ...
    'BUILD_AIR_TOO_SMALL', ...
        (Pv.airrad - abs(Pv.airZOffset)) - max(LMax, Pv.th/2), 'build', true};
if Pv.d > 0
    % Sidewall angle: fillets are zero (a builder error otherwise), and the
    % rib and spine top faces must keep a positive width.
    inset = Pv.th*tan(Pv.d);
    rules = [rules; { ...
        'BUILD_SIDEWALL_TOO_STEEP', min(tMin, Pv.s*onesN) - 2*inset, 'build', true}];
else
    rules = [rules; { ...
        'BUILD_TIP_FILLET',      tMin - 2*Pv.r2,                    'build', true; ...
        'BUILD_ROOT_FILLET',     gapE/2 - Pv.r1,                    'build', true; ...
        'WARN_FILLETS_MEET',     (LMin - Pv.s/2) - (Pv.r1 + Pv.r2), 'warn',  true; ...
        'WARN_STAGING_OVERLAP',  4*LMin + Pv.s/2 - 10e-9 - LMax,    'warn',  true}];
    if Pv.r1 > 0
        rules = [rules; { ...
            'WARN_ROOT_FILLET_TIGHT', clrMin - Pv.r1,               'warn',  true}];
    end
end
tAll = min(min(tLoc, [], 2), tMin);
LAll = min(min(LLoc, [], 2), LMin);
gAll = min(min(aLoc - tLoc, [], 2), gapE);
rules = [rules; { ...
    'FAB_MIN_FEATURE',   min([tAll, Pv.s*onesN, LAll - Pv.s/2], [], 2) - minFeature, 'fab', false; ...
    'FAB_MIN_GAP',       gAll - minGap,                                              'fab', false; ...
    'FAB_MIN_THICKNESS', Pv.th*onesN - minThickness,                                 'fab', false}];

valid  = ~bad;
warnOK = true(n, 1);
fabOK  = true(n, 1);
ok     = ~bad;
margin = inf(n, 1);
why    = repmat({''}, n, 1);
for r = 1:size(rules, 1)
    mk = rules{r, 2};
    if rules{r, 4}
        pass = mk > 0;
    else
        pass = mk >= 0;
    end
    pass(bad) = false;
    active = true;
    switch rules{r, 3}
        case 'build'
            valid = valid & pass;
        case 'warn'
            warnOK = warnOK & pass;
            active = logical(strictWarn);
        case 'fab'
            fabOK = fabOK & pass;
    end
    if ~active
        continue;
    end
    ok = ok & pass;
    tighter = mk < margin;
    margin(tighter) = mk(tighter);
    why(tighter)    = rules(r, 1);
end

% A design whose parameters cannot even be mapped is infeasible by a wide
% margin, so a graded penalty pushes strongly away from it.
margin(bad) = -Inf;
why(bad)    = {'NONPHYSICAL_TAPER'};
why(ok)     = {''};

ok     = reshape(ok, outSize);
margin = reshape(margin, outSize);
why    = reshape(why, outSize);
if nargout > 3
    info = struct( ...
        'valid',  reshape(valid, outSize), ...
        'warn',   reshape(~warnOK & ~bad, outSize), ...
        'fab',    reshape(fabOK, outSize), ...
        'oblong', reshape(oblong, outSize), ...
        'a_ctr',  reshape(Pv.a*sigC, outSize), ...
        't_ctr',  reshape(Pv.t*sigC.^(1 - oblong), outSize), ...
        'w_ctr',  reshape(2*Lm*sigC.^(1 + oblong), outSize));
end
end

% =========================================================================
function v = getfielddef(S, f, d)
if isstruct(S) && isfield(S, f) && ~isempty(S.(f))
    v = S.(f);
else
    v = d;
end
end

function [A, B] = expandPair(a, b)
%EXPANDPAIR  Broadcast two inputs to a common size.
if isscalar(a) && ~isscalar(b)
    A = a*ones(size(b));
    B = b;
elseif isscalar(b) && ~isscalar(a)
    A = a;
    B = b*ones(size(a));
else
    if ~isequal(size(a), size(b))
        error('isFabricableRib:sizeMismatch', ...
            'defectAspectRatio and maxdef must be the same size or scalar.');
    end
    A = a;
    B = b;
end
end

function f = taperProfile(name)
%TAPERPROFILE  Defect profile f(x), identical to createRibCavityGeom.m.
switch name
    case 'quadratic', f = @(x) 1 - x.^2;
    case 'cubic',     f = @(x) 1 - 3*x.^2 + 2*x.^3;
    case 'linear',    f = @(x) 1 - x;
    case 'none',      f = @(x) ones(size(x));
    otherwise
        error('isFabricableRib:badTaperFunc', ...
            'P.taperFunc = ''%s'' is not one of quadratic, cubic, linear, none.', name);
end
end
