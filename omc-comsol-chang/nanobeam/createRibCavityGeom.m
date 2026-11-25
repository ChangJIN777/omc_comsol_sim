function P = createRibCavityGeom(P)
%CREATERIBCAVITYGEOM  Per-cell dimensions and positions of a rib cavity.
%
%   P = createRibCavityGeom(P)
%
%   The 'rib' counterpart of CreateNanobeamGeom.m: pure MATLAB, called
%   by RunNanobeamFEM before buildRibCavity. Defaults and validation
%   come from ribDefaultParams / validateRibParams in the
%   omc-comsol-chang root folder, which must be on the path (the nanobeam
%   test script adds it with addpath(...,'-end')).
%
%   Cells are numbered k = 1 .. ncellHalf outward from the centre. A rib at
%   the centre (P.ribatctr = 1) is cell 1, the analogue of a hole at the
%   centre (holeatctr = 1) in CreateNanobeamGeom.
%
%   RIB DIMENSIONS (as hx/hy in CreateNanobeamGeom, oblong = 0)
%     k <= ndef:  a_k = a*(1 - maxdef*f((k-1)/ndef))      defect cells
%     k >  ndef:  a_k = a                                  mirror cells
%     L_k = ribLength(P, a_k): an explicit P.L (or P.w/2) for EVERY cell if
%           given; else, with P.w_mirror (+ optional P.w_ctr), the tip-to-tip
%           width 2*L_k tapered linearly in a_k from w_mirror (a_k = a) to
%           w_ctr (a_k = a*(1 - maxdef)); otherwise the rule
%           L_ref + dLda*(a_k - a_ref)
%     The sidewall angle P.d is passed through unchanged to buildRibCavity.
%     f = 1-x^2 (quadratic), 1-3x^2+2x^3 (cubic), 1-x (linear), 1 (none: a
%     step to a*(1 - maxdef) inside the defect), as CreateNanobeamGeom.
%
%   HOLE-STYLE TAPER, P.ribTaperScheme = 'hole' (default 'rib' = the above)
%     CreateNanobeamGeom's maxdef/oblong scheme with hy <-> w = 2L, hx <-> t.
%     With sig_k = a_k/a = 1 - maxdef*f((k-1)/ndef)  (1 for mirror cells):
%         t_k = t  *sig_k^(1 - oblong)
%         w_k = w_m*sig_k^(1 + oblong),   L_k = w_k/2
%     w_m = mirror-cell width from P.w_mirror or the L(a) rule at a = P.a.
%     oblong = 0 scales t and w like a; 1 keeps t; -1 keeps w; > 1 widens t
%     and shortens the ribs toward the centre; < -1 lengthens them. P.oblong
%     defaults to 2.7 (ribDefaultParams); P.w_ctr and an explicit P.L / P.w
%     are errors. s, r1, r2, th and d stay constant. Positions follow
%     P.taperConvention as below.
%
%   POSITIONS, P.taperConvention
%     'centerspacing' (DEFAULT; CreateNanobeamGeom's convention). The spacing
%         sp_k is the centre-to-centre distance from cell k-1 to cell k, and
%             x_k = x_{k-1} + sp_k      (i.e. x_{n+1} = x_n + a_{n+1})
%         with the spacing index shifted by one when a rib sits on the centre:
%             k <= ndef:  xa_k = (k-1-ribatctr)/(ndef-ribatctr)
%                         sp_k = a*(1 - maxdef*f(xa_k)),  NaN if xa_k < 0
%             k >  ndef:  sp_k = a
%             ribatctr = 1:  x = [0, cumsum(sp_2 .. sp_N)]   (sp_1 = NaN, unused)
%             ribatctr = 0:  x = cumsum([sp_1/2, sp_2 .. sp_N])
%         These are CreateNanobeamGeom.m's lines for a_hole and xposHalf, with
%         holeatctr -> ribatctr. The first spacing out of the centre is the
%         full defect a*(1 - maxdef). Each spine segment runs between the
%         midpoints of neighbouring centres (so it is not centred on its rib
%         inside the taper), and the beam ends sp_N/2 past the last centre, as
%         CreateNanobeamGeom's beamLenHalf = xpos(end) + a_hole(end)/2.
%         With the defaults: x = 0, 240, 483.75, 738.75, 1012.5, 1312.5, then
%         +300 nm per cell, up to 3712.5 nm. That is 10.5 nm beyond the
%         template from cell 6 outward, and the beam is 7.725 um long.
%     'cellwidth' (the template's convention). a_k is the WIDTH of cell k's
%         spine segment, segments butt together, and the rib is centred in it:
%             x_1 = 0 (ribatctr) or a_1/2,   x_{k+1} = x_k + (a_k + a_{k+1})/2
%         With the defaults: x = 0, 241.2, 487.2, 742.8, 1012.8, 1302, 1602,
%         ..., 3702 nm and a 7.704 um beam, which are the template's cell positions
%         and air-cylinder length exactly.
%
%   OUTPUT FIELDS
%     P.ribGeomHalf  cells with x >= 0 (centre once), ordered outward
%     P.ribGeom      full beam, mirrored about x = 0
%     P.ribGeomCols  column indices: x (rib centre), a (a_k, sets L_k), t, L,
%                    s, xL / xR (spine segment edges)
%     P.beamLenHalf, P.beamLen
%     P.a_cell, P.L_cell, P.t_cell   half-beam a_k, L_k and t_k
%     P.sp_cell      half-beam centre spacings sp_k ('centerspacing'), or the
%                    actual centre-to-centre distances ('cellwidth')
%     P.w            2*max(L_k), the beam's full y extent. Read as the beam
%                    width by SetupNanobeamFEM / SolveNanobeamFEM. Written
%                    back as an output and recorded in P.ribDerived.w, so a
%                    second call does not mistake it for an explicit width.
%     P.w_hole       2*L_k per half-beam cell, read by RunNanobeamFEM's strain
%                    interpolation (max(w_hole)/2 = y range)
%
%   To hand-tune individual cells (the template's author did this for the
%   centre rib), edit P.ribGeomHalf / P.ribGeom after this call and before the
%   build, the same way RunNanobeamFEM applies P.shave.

if exist('ribDefaultParams','file') ~= 2
    error('rib:pathMissing', ...
        ['ribDefaultParams.m (omc-comsol-chang root) is not on the path. ' ...
         'Add it with addpath(<root>,''-end''), as test_nanobeamRibFEM.m does.']);
end
P = ribDefaultParams(P);
validateRibParams(P, 'cavity');

N = P.ncellHalf;
k = (1:N)';
h = double(P.ribatctr);
f = taperFunc(P.taperFunc);
inDef = k <= P.ndef;

% rib dimensions: defect cells tapered with x = (k-1)/ndef, mirror cells = a
aLoc = P.a*ones(N,1);
if any(inDef)
    aLoc(inDef) = P.a*(1 - P.maxdef*f((k(inDef)-1)/P.ndef));
end
% L_k and t_k: t_k = P.t unless P.ribTaperScheme = 'hole' (see ribLength)
[LLoc, Lsrc, tLoc] = ribLength(P, aLoc);
sLoc = P.s*ones(N,1);

switch P.taperConvention
    case 'centerspacing'
        sp = P.a*ones(N,1);
        idx = inDef & (k-1-h >= 0);          % ndef-h >= 1 whenever idx holds
        sp(idx) = P.a*(1 - P.maxdef*f((k(idx)-1-h)/(P.ndef-h)));
        if h
            sp(1) = NaN;                     % no spacing before the centre cell
            xc = [0; cumsum(sp(2:end))];
        else
            xc = cumsum([sp(1)/2; sp(2:end)]);
        end
        mid = (xc(1:end-1) + xc(2:end))/2;
        if h
            if N > 1
                leftEdge = -mid(1);          % centre segment symmetric about x = 0
            else
                leftEdge = -aLoc(1)/2;       % lone centre cell: no spacing defined
            end
        else
            leftEdge = 0;                    % midpoint of the two central ribs
        end
        lastSp = sp(end);
        if isnan(lastSp)
            lastSp = aLoc(end);
        end
        xL = [leftEdge; mid];
        xR = [mid; xc(end) + lastSp/2];

    case 'cellwidth'
        if h
            x0 = 0;
        else
            x0 = aLoc(1)/2;
        end
        xc = x0 + [0; cumsum((aLoc(1:end-1) + aLoc(2:end))/2)];
        xL = xc - aLoc/2;
        xR = xc + aLoc/2;
        sp = [NaN; diff(xc)];
        if ~h
            sp(1) = 2*xc(1);
        end

    otherwise
        error('rib:badTaperConvention', ...
            'P.taperConvention = ''%s'' must be ''centerspacing'' or ''cellwidth''.', ...
            P.taperConvention);
end

% Per cell: the rib must sit inside its own spine segment on BOTH sides
% ('centerspacing' segments are not centred on their rib, and t_k varies in
% the 'hole' scheme), with room left for the root fillet.
clr = min(xc - xL, xR - xc) - tLoc/2;
[clrMin, kBad] = min(clr);
if clrMin <= 0
    error('rib:segmentTooNarrow', ...
        ['Cell %d: rib (t = %g nm at x = %g nm) is not narrower than its spine ' ...
         'segment [%g, %g] nm.'], kBad, tLoc(kBad)*1e9, xc(kBad)*1e9, ...
        xL(kBad)*1e9, xR(kBad)*1e9);
end
if P.d == 0 && P.r1 > 0 && clrMin <= P.r1
    warning('rib:rootFilletTight', ...
        ['Cell %d: only %g nm of spine beside the rib (t = %g nm) for the root ' ...
         'fillet P.r1 = %g nm. COMSOL may refuse the fillet.'], ...
        kBad, clrMin*1e9, tLoc(kBad)*1e9, P.r1*1e9);
end

C = struct('x',1,'a',2,'t',3,'L',4,'s',5,'xL',6,'xR',7);
P.ribGeomCols = C;
P.ribGeomHalf = [xc aLoc tLoc LLoc sLoc xL xR];

if h
    mirrorRows = flipud(P.ribGeomHalf(2:end,:));   % centre cell appears once
else
    mirrorRows = flipud(P.ribGeomHalf);
end
mirrorRows(:,C.x) = -mirrorRows(:,C.x);
mirrorRows(:,[C.xL C.xR]) = -mirrorRows(:,[C.xR C.xL]);
P.ribGeom = [mirrorRows; P.ribGeomHalf];

P.beamLenHalf = xR(end);
P.beamLen     = 2*P.beamLenHalf;
P.a_cell      = aLoc;
P.L_cell      = LLoc;
P.t_cell      = tLoc;
P.sp_cell     = sp;
P.w           = 2*max(LLoc);
if ~strcmp(Lsrc,'w')
    P.ribDerived.w = P.w;     % computed here, not the caller's explicit width
end
P.w_hole      = 2*LLoc;
end

% -------------------------------------------------------------------------

function f = taperFunc(name)
%TAPERFUNC  Defect profile f(x) on 0 <= x < 1, identical to CreateNanobeamGeom.m.
switch name
    case 'quadratic', f = @(x) 1 - x.^2;
    case 'cubic',     f = @(x) 1 - 3*x.^2 + 2*x.^3;
    case 'linear',    f = @(x) 1 - x;
    case 'none',      f = @(x) ones(size(x));
    otherwise
        error('rib:badTaperFunc', ...
            'P.taperFunc = ''%s'' is not one of quadratic, cubic, linear, none.', name);
end
end
