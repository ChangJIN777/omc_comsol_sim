function validateRibParams(P, mode)
%VALIDATERIBPARAMS  Reject impossible 'rib' geometry before COMSOL sees it.
%
%   validateRibParams(P, mode)      mode = 'unitcell' | 'cavity'
%
%   P must already carry the defaults (ribDefaultParams). Hard errors are
%   reserved for input that would make COMSOL fail, give a disconnected
%   solid, or silently build the wrong thing; questionable but buildable
%   input only warns. The full table is in rib_README.md.

posScalars = {'a','t','s','th','L_ref','a_ref'};
for k = 1:numel(posScalars)
    validateattributes(P.(posScalars{k}), {'numeric'}, ...
        {'scalar','real','finite','positive'}, mfilename, ['P.',posScalars{k}]);
end
validateattributes(P.dLda, {'numeric'}, {'scalar','real','finite'}, mfilename, 'P.dLda');
validateattributes(P.r1, {'numeric'}, {'scalar','real','finite','nonnegative'}, mfilename, 'P.r1');
validateattributes(P.r2, {'numeric'}, {'scalar','real','finite','nonnegative'}, mfilename, 'P.r2');

% Sidewall angle, radians. A value >= pi/2 is almost certainly degrees.
if ~(isnumeric(P.d) && isscalar(P.d) && isreal(P.d) && isfinite(P.d) && ...
        P.d >= 0 && P.d < pi/2)
    error('rib:badSidewallAngle', ...
        ['P.d must be a sidewall angle in RADIANS in [0, pi/2), e.g. 15*pi/180; ' ...
         'got %s.'], mat2str(P.d));
end

switch mode
    case 'unitcell'
        aMin = P.a;
        aMax = P.a;
    case 'cavity'
        validateattributes(P.ncellHalf, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'P.ncellHalf');
        validateattributes(P.ndef, {'numeric'}, {'scalar','integer','nonnegative'}, mfilename, 'P.ndef');
        validateattributes(P.maxdef, {'numeric'}, {'scalar','real','>=',0,'<',1}, mfilename, 'P.maxdef');
        validateattributes(P.ribatctr, {'numeric','logical'}, {'scalar','binary'}, mfilename, 'P.ribatctr');
        validateattributes(P.airrad, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.airrad');
        validateattributes(P.airZOffset, {'numeric'}, {'scalar','real','finite'}, mfilename, 'P.airZOffset');
        if ~any(strcmp(P.taperConvention, {'centerspacing','cellwidth'}))
            error('rib:badTaperConvention', ...
                'P.taperConvention = ''%s'' must be ''centerspacing'' or ''cellwidth''.', ...
                P.taperConvention);
        end
        if ~any(strcmp(P.mEndBC, {'lowreflecting','fixed'}))
            error('rib:badEndBC', ...
                'P.mEndBC = ''%s'' must be ''lowreflecting'' or ''fixed''.', P.mEndBC);
        end
        aMin = P.a*(1 - P.maxdef);
        aMax = P.a;
    otherwise
        error('rib:badMode', 'mode must be ''unitcell'' or ''cavity''.');
end

% L is either constant (explicit P.L / P.w) or linear in a, so its extremes
% over [aMin, aMax] sit at the ends.
Lends = ribLength(P, [aMin aMax]);
LMin = min(Lends);
LMax = max(Lends);

% Rib must be narrower than the narrowest cell, or neighbouring ribs fuse.
if P.t >= aMin
    error('rib:ribTooWide', ...
        'P.t = %g nm must be smaller than the narrowest cell width %g nm.', ...
        P.t*1e9, aMin*1e9);
end
% Rib must protrude past the spine (also at the top face when d > 0, since
% both the tip and the spine edge move in by th*tan(d)).
if LMin <= P.s/2
    error('rib:ribShorterThanSpine', ...
        'Shortest rib half-length %g nm must exceed the spine half-width %g nm.', ...
        LMin*1e9, P.s/2*1e9);
end

if P.d > 0
    % Sidewall angle: walls lean in by th*tan(d) over the slab height.
    inset = P.th*tan(P.d);
    if P.r1 > 0 || P.r2 > 0
        error('rib:filletWithSidewall', ...
            ['In-plane fillets (P.r1 = %g nm, P.r2 = %g nm) cannot be combined ' ...
             'with a sidewall angle (P.d = %g deg). Set P.r1 = P.r2 = 0, which ' ...
             'is the default when d ~= 0.'], P.r1*1e9, P.r2*1e9, P.d*180/pi);
    end
    if P.t - 2*inset <= 0
        error('rib:sidewallTooSteep', ...
            'Rib top width t - 2*th*tan(d) = %g nm is not positive.', (P.t - 2*inset)*1e9);
    end
    if P.s - 2*inset <= 0
        error('rib:sidewallTooSteep', ...
            'Spine top width s - 2*th*tan(d) = %g nm is not positive.', (P.s - 2*inset)*1e9);
    end
    % Angled walls are not mirror-symmetric in z.
    if strcmp(mode,'unitcell') && isfield(P,'mbevenz') && P.mbevenz ~= 0
        error('rib:zSymWithSidewall', ...
            ['P.mbevenz = %g asks for a z-symmetry plane, but with P.d ~= 0 the ' ...
             'slab is not mirror-symmetric in z. Use P.mbevenz = 0 (solveBands: ' ...
             'TwoSymPlanes = 0, zSymCondition = 0; or solveBands_noSym).'], P.mbevenz);
    end
else
    % Fillet checks (d = 0 only; with d > 0 the fillets are zero, see above).
    if 2*P.r2 >= P.t
        error('rib:tipFilletTooLarge', ...
            '2*P.r2 = %g nm must be smaller than P.t = %g nm.', 2*P.r2*1e9, P.t*1e9);
    end
    if P.r1 >= (aMin - P.t)/2
        error('rib:rootFilletTooLarge', ...
            'P.r1 = %g nm does not fit in the %g nm of spine beside the rib.', ...
            P.r1*1e9, (aMin - P.t)/2*1e9);
    end
    if P.r1 + P.r2 >= LMin - P.s/2
        warning('rib:filletsMeet', ...
            ['P.r1 + P.r2 = %g nm leaves no straight rib flank (shortest flank %g nm). ' ...
             'COMSOL may refuse the second fillet.'], (P.r1+P.r2)*1e9, (LMin-P.s/2)*1e9);
    end
    % The d = 0 path stages each cell at y0 = 4*L before moving it in (see
    % addRibCell). Its lowest fillet box must clear the cells already placed.
    selDelta = 10e-9;
    if 4*LMin + P.s/2 - selDelta <= LMax
        warning('rib:stagingOverlap', ...
            ['Shortest rib (L = %g nm) is staged at y0 = %g nm, inside the band ' ...
             'already occupied by placed cells (|y| <= %g nm). Its fillet ' ...
             'selections may catch a neighbour.'], LMin*1e9, 4*LMin*1e9, LMax*1e9);
    end
end

if strcmp(mode,'cavity')
    if P.airrad - abs(P.airZOffset) <= max(LMax, P.th/2)
        error('rib:airTooSmall', ...
            ['P.airrad = %g nm (axis offset %g nm) does not enclose the beam ' ...
             '(L up to %g nm).'], P.airrad*1e9, P.airZOffset*1e9, LMax*1e9);
    end
elseif isfield(P,'run_optical') && P.run_optical
    validateattributes(P.airradUC, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'P.airradUC');
    if P.airradUC <= max(LMax, P.th/2)
        error('rib:airTooSmall', ...
            'P.airradUC = %g nm does not enclose the rib (L = %g nm).', ...
            P.airradUC*1e9, LMax*1e9);
    end
end
end
