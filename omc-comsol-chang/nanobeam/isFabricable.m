function [ok, margin, why] = isFabricable(defectAspectRatio, maxdef, P)
%ISFABRICABLE  Lithography feasibility of a taper design, without any solve.
%
%   [ok, margin, why] = isFabricable(defectAspectRatio, maxdef, P)
%
%   Reproduces the fabrication-tolerance guards that CreateNanobeamGeom throws
%   on (see CreateNanobeamGeom.m:72-84), but as pure arithmetic on the taper
%   parameters -- no geometry build, no COMSOL, microseconds instead of a
%   failed multi-minute solve.
%
%   INPUTS
%     defectAspectRatio - hole aspect ratio at the cavity centre. Scalar or
%                         array; expanded against maxdef.
%     maxdef            - defect ratio. Scalar or array.
%     P                 - parameter struct; uses P.a, P.w, P.hx, P.hy and,
%                         optionally, the tolerance overrides below.
%
%   OUTPUTS
%     ok     - logical array, true where the design is fabricable
%     margin - worst binding margin in METRES: how much slack the tightest
%              rule has. Negative means violated, and its magnitude says by
%              how much -- which lets a caller build a graded penalty rather
%              than a flat one.
%     why    - cell array of strings naming the tightest rule (empty where ok)
%
%   TOLERANCES (defaults match CreateNanobeamGeom exactly; override via P)
%     P.minSidewallGap - beam edge to hole edge, default 148 nm
%     P.minHoleGap     - gap between adjacent holes,  default  50 nm
%     P.minFeature     - smallest hole dimension,     default  50 nm
%
%   The checks are applied at the CAVITY CENTRE, where the taper is most
%   extreme, as well as at the nominal mirror cell:
%     hx_def = hx*(1-maxdef)^(1-oblong)   grows when oblong > 1
%     hy_def = hy*(1-maxdef)^(1+oblong)   shrinks when oblong > -1
%     a_def  = a *(1-maxdef)
%   so a_def - hx_def is normally the binding rule: the lattice contracts
%   while the hole widens, and the solid bridge between holes vanishes.
%
%   Example -- prune a candidate before paying for a solve:
%     if ~isFabricable(1.25, 0.22, P), return; end
%
%   See also CREATENANOBEAMGEOM, OPTIMIZE_OBLONG_MAXDEF.

narginchk(3, 3);
validateattributes(defectAspectRatio, {'numeric'}, {'real'}, mfilename, 'defectAspectRatio');
validateattributes(maxdef,            {'numeric'}, {'real'}, mfilename, 'maxdef');

% Tolerances -- defaulted to CreateNanobeamGeom's hard-coded values so this
% helper can never disagree with the code that actually throws.
minSidewallGap = getfielddef(P, 'minSidewallGap', 148e-9);
minHoleGap     = getfielddef(P, 'minHoleGap',      50e-9);
minFeature     = getfielddef(P, 'minFeature',      50e-9);

% Expand to a common size so scalar/array mixes work.
[dAR, md] = expandPair(defectAspectRatio, maxdef);

a  = P.a;   w  = P.w;
hx = P.hx;  hy = P.hy;

% --- map the aspect ratio back to the taper exponent -------------------
% Identical to the mapping in optimize_oblong_maxdef / sweep_oblong_maxdef.
mirrorAspectRatio = hy / hx;
r      = dAR ./ mirrorAspectRatio;
oneMd  = 1 - md;
oblong = log(r) ./ (2 * log(oneMd));

% A non-physical mapping (maxdef outside (0,1), or a non-positive aspect
% ratio) cannot produce a geometry at all -- treat it as infeasible rather
% than letting a complex or NaN oblong propagate into the geometry builder.
bad = ~isfinite(oblong) | ~isreal(oblong) | md <= 0 | md >= 1 | dAR <= 0;

% --- dimensions at the cavity centre -----------------------------------
a_def  = a  .* oneMd;
hx_def = hx .* oneMd .^ (1 - oblong);
hy_def = hy .* oneMd .^ (1 + oblong);

% --- the rules, each as a signed margin in metres -----------------------
% CreateNanobeamGeom checks the sidewall and minimum-feature rules on the
% NOMINAL hole and the hole-gap rule at the DEFECT. Both locations are
% checked here for every rule: the taper can make either one binding, and a
% design that only fails at the centre would otherwise slip through.
ones_ = ones(size(dAR));
rules = { ...
    'sidewall gap (mirror)', (w - hy) * ones_          - minSidewallGap; ...
    'sidewall gap (defect)', (w - hy_def)              - minSidewallGap; ...
    'hole gap (defect)',     (a_def - hx_def)          - minHoleGap;     ...
    'hx min feature',        min(hx * ones_, hx_def)   - minFeature;     ...
    'hy min feature',        min(hy * ones_, hy_def)   - minFeature};

nRules = size(rules, 1);
margin = inf(size(dAR));
why    = repmat({''}, size(dAR));

for k = 1:nRules
    mk = rules{k, 2};
    tighter = mk < margin;
    margin(tighter) = mk(tighter);
    why(tighter)    = rules(k, 1);
end

% A design whose parameters cannot even be mapped is infeasible by a wide
% margin, so a graded penalty pushes strongly away from it.
margin(bad) = -inf;
why(bad)    = {'non-physical taper parameters'};

ok = margin >= 0;
why(ok) = {''};
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
    A = a * ones(size(b));  B = b;
elseif isscalar(b) && ~isscalar(a)
    A = a;  B = b * ones(size(a));
else
    if ~isequal(size(a), size(b))
        error('isFabricable:sizeMismatch', ...
            'defectAspectRatio and maxdef must be the same size or scalar.');
    end
    A = a;  B = b;
end
end
