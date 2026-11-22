function [L, src] = ribLength(P, a)
%RIBLENGTH  Rib half-length (beam axis to rib tip) for lattice value(s) a.
%
%   [L, src] = ribLength(P, a)
%
%   Precedence, highest first:
%     1. P.L  explicit half-length       -> L = P.L for every a   (src 'L')
%     2. P.w  explicit tip-to-tip width  -> L = P.w/2 for every a (src 'w')
%     3. the linear rule  L(a) = P.L_ref + P.dLda*(a - P.a_ref)    (src 'rule')
%   Giving both P.L and P.w with P.w ~= 2*P.L is an error.
%
%   The builders write P.L / P.w back as OUTPUTS and record them in
%   P.ribDerived. A P.L / P.w equal to its recorded derived value is
%   ignored here, so re-running a builder on its own output, or after
%   changing P.a, still follows the rule. Setting a NEW value makes it
%   explicit again.
%
%   L has the size of a. Errors if any L is not positive.

hasL = isExplicit(P, 'L');
hasW = isExplicit(P, 'w');
if hasL && hasW && abs(P.w - 2*P.L) > 1e-15
    error('rib:inconsistentLW', ...
        'Both P.L = %g nm and P.w = %g nm are set but P.w ~= 2*P.L.', ...
        P.L*1e9, P.w*1e9);
end
if hasL
    L = P.L*ones(size(a));
    src = 'L';
elseif hasW
    L = P.w/2*ones(size(a));
    src = 'w';
else
    L = P.L_ref + P.dLda*(a - P.a_ref);
    src = 'rule';
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
