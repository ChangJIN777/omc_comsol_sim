function [L, src, t] = ribLength(P, a)
%RIBLENGTH  Rib half-length (beam axis to rib tip) for lattice value(s) a.
%
%   [L, src, t] = ribLength(P, a)
%
%   Precedence, highest first:
%     1. P.L  explicit half-length       -> L = P.L for every a   (src 'L')
%     2. P.w  explicit tip-to-tip width  -> L = P.w/2 for every a (src 'w')
%     3. P.w_mirror (+ optional P.w_ctr): tip-to-tip width 2L of the mirror
%        cell (a = P.a) and, optionally, of the centre cell
%        (a = P.a*(1 - P.maxdef)), linear in a in between:
%            L(a) = P.w_mirror/2 + slope*(a - P.a)                  (src 'w_mirror')
%        slope = (P.w_mirror - P.w_ctr)/(2*P.a*P.maxdef) with P.w_ctr,
%        P.dLda without it. maxdef = 0 needs P.w_ctr == P.w_mirror (slope 0).
%     4. the linear rule  L(a) = P.L_ref + P.dLda*(a - P.a_ref)    (src 'rule')
%   Giving both P.L and P.w with P.w ~= 2*P.L is an error, as is P.w_mirror
%   together with an explicit P.L or P.w, or P.w_ctr without P.w_mirror.
%
%   The builders write P.L / P.w back as OUTPUTS and record them in
%   P.ribDerived. A P.L / P.w equal to its recorded derived value is
%   ignored here, so re-running a builder on its own output, or after
%   changing P.a, still follows the rule. Setting a NEW value makes it
%   explicit again.
%
%   HOLE-STYLE TAPER (P.ribTaperScheme = 'hole'; default 'rib' = all of the
%   above, unchanged). The source above gives the MIRROR-cell length
%   Lm = L(P.a), and CreateNanobeamGeom's maxdef/oblong scaling is applied
%   with hy <-> w = 2L and hx <-> t:
%       L(a) = Lm *(a/P.a)^(1 + P.oblong),   t(a) = P.t*(a/P.a)^(1 - P.oblong)
%   (a/P.a = 1 - maxdef*f for a cavity cell, 1 for a mirror or unit cell).
%   P.oblong defaults to 2.7 (ribDefaultParams); a positive oblong shortens
%   the ribs and widens them toward the centre. P.w_ctr is an error (oblong
%   sets the centre width). src still names the source of Lm.
%
%   t (optional third output) is the rib width along x for each a: P.t in
%   the 'rib' scheme, tapered as above in the 'hole' scheme. Needs P.t.
%
%   L has the size of a. Errors if any L is not positive.

hasL = isExplicit(P, 'L');
hasW = isExplicit(P, 'w');
if hasL && hasW && abs(P.w - 2*P.L) > 1e-15
    error('rib:inconsistentLW', ...
        'Both P.L = %g nm and P.w = %g nm are set but P.w ~= 2*P.L.', ...
        P.L*1e9, P.w*1e9);
end
hasWm = isfield(P,'w_mirror') && ~isempty(P.w_mirror);
hasWc = isfield(P,'w_ctr') && ~isempty(P.w_ctr);
if hasWm && (hasL || hasW)
    error('rib:conflictingLength', ...
        'P.w_mirror cannot be combined with an explicit P.L or P.w; set only one.');
end
if hasWc && ~hasWm
    error('rib:wCtrWithoutMirror', 'P.w_ctr needs P.w_mirror.');
end
holeScheme = isfield(P,'ribTaperScheme') && strcmp(P.ribTaperScheme,'hole');
if holeScheme
    oblong = 2.7;                 % default; keep equal to ribDefaultParams
    if isfield(P,'oblong') && ~isempty(P.oblong)
        oblong = P.oblong;
    end
    if ~(isnumeric(oblong) && isscalar(oblong) && isreal(oblong) && isfinite(oblong))
        error('rib:badOblong', ...
            'P.oblong must be a finite real scalar (P.ribTaperScheme = ''hole'').');
    end
    if hasWc
        error('rib:wCtrWithHoleScheme', ...
            ['P.w_ctr cannot be combined with P.ribTaperScheme = ''hole'': the ' ...
             'centre width is w_mirror*(1 - maxdef)^(1 + oblong). Remove P.w_ctr.']);
    end
    aBase = P.a*ones(size(a));    % the sources below give the mirror-cell length
else
    aBase = a;
end
if hasL
    L = P.L*ones(size(a));
    src = 'L';
elseif hasW
    L = P.w/2*ones(size(a));
    src = 'w';
elseif hasWm
    if ~hasWc
        slope = P.dLda;
    elseif P.maxdef > 0
        slope = (P.w_mirror - P.w_ctr)/(2*P.a*P.maxdef);
    elseif abs(P.w_ctr - P.w_mirror) > 1e-15
        error('rib:wCtrNoTaper', ...
            ['P.w_ctr = %g nm differs from P.w_mirror = %g nm but P.maxdef = 0, ' ...
             'so there is no centre cell to taper to.'], P.w_ctr*1e9, P.w_mirror*1e9);
    else
        slope = 0;
    end
    L = P.w_mirror/2 + slope*(aBase - P.a);
    src = 'w_mirror';
else
    L = P.L_ref + P.dLda*(aBase - P.a_ref);
    src = 'rule';
end
if holeScheme
    sc = a/P.a;
    L = L.*sc.^(1 + oblong);
end
if nargout > 2
    t = P.t*ones(size(a));
    if holeScheme
        t = t.*sc.^(1 - oblong);
    end
end
if any(L(:) <= 0)
    error('rib:nonPositiveL', ...
        ['Rib half-length is not positive (%g nm at a = %g nm; source ''%s''; ' ...
         'L_ref = %g nm, a_ref = %g nm, dLda = %g).'], ...
        L(find(L(:) <= 0, 1))*1e9, a(find(L(:) <= 0, 1))*1e9, src, ...
        P.L_ref*1e9, P.a_ref*1e9, P.dLda);
end
end

% -------------------------------------------------------------------------

function tf = isExplicit(P, f)
% Present, non-empty, and not just the value a builder wrote back.
tf = isfield(P, f) && ~isempty(P.(f));
if tf && isfield(P,'ribDerived') && isfield(P.ribDerived, f)
    tf = ~isequal(P.(f), P.ribDerived.(f));
end
end
