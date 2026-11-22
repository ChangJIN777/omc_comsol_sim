function P = ribDefaultParams(P)
%RIBDEFAULTPARAMS  Fill in missing fields for the 'rib' celltype.
%
%   P = ribDefaultParams(P)
%   P = ribDefaultParams()          % all defaults
%
%   Only fields that are ABSENT are set, so a script can override any of
%   them first. Every builder, solver branch and file-name branch for 'rib'
%   calls this, so the defaults live in exactly one place. Used by BOTH
%   pipelines: nanobeam/test_nanobeamRibFEM.m puts this folder on the path
%   with addpath(...,'-end') for that reason. The full parameter reference is
%   rib_README.md.
%
%   GEOMETRY (SI units; beam axis = x, ribs along +-y, thickness = z)
%     A continuous spine of width s carries one rib per cell, centred on it.
%     A rib is t wide in x and reaches |y| = L from the beam axis. Fillet r2
%     rounds the two convex rib-tip corners and r1 the two concave rib/spine
%     corners, per half. With a sidewall angle d the walls lean inward toward
%     the top face (z = +th/2), as in the former LN rib builder.
%
%   FIELD         DEFAULT        SOURCE
%   celltype      'rib'
%   xsect         'rect'
%   beamMat       'diamond'
%   anisoMat      1              visible-rib template: anisotropic (Voigt)
%   rxtal         0 deg          template stiffness tensor is unrotated
%   a             300 nm         template mirror pitch
%   t             100 nm         template; rib_optimize_sweep_diamond.m
%   s             100 nm         rib_optimize_sweep_diamond.m; template IR set
%   L_ref         400 nm         rib half-length L at a = a_ref (template centre)
%   a_ref         240 nm         anchor of the L(a) rule, NOT the lattice constant
%   dLda          -0.4998        slope dL/da (template IR set)
%   d             0 rad          sidewall angle in RADIANS (as the LN scripts:
%                                P.d = 15*pi/180). 0 = vertical walls.
%   th            110 nm         EFFECTIVE thickness of the visible template
%                                (z-offset defect, see rib_README.md); nominal 120
%   r1, r2        40 nm if d = 0, 0 if d ~= 0
%                                template IR set (TODO: unknown for the visible
%                                set). Fillets cannot be combined with d ~= 0, so
%                                they default to 0 there: the LN geometry.
%   ncellHalf     14             cavity: cells in x >= 0, centre counted once
%   ribatctr      1              cavity: rib centred on x = 0
%   ndef          5              cavity: defect cells
%   maxdef        0.2            cavity: centre a = a*(1 - maxdef)
%   taperFunc     'quadratic'    cavity taper profile
%   taperConvention 'centerspacing'  cavity: 'centerspacing' | 'cellwidth'
%   airrad        3.1 um         cavity air-cylinder radius
%   airZOffset    60 nm          cavity air-cylinder axis above the mid-plane
%   airradUC      5*a            unit-cell air radius when P.run_optical (LN)
%   mEndBC        'lowreflecting'  cavity beam-end BC ('lowreflecting'|'fixed')
%   lbd0          1.55 um        model parameter lbd0, read by the LN refractive
%                                index in runOpticalBand_1D/_3D (LN builder value)
%   holeatedge 0, plotgeom 0     read unconditionally by runBands
%
%   RIB LENGTH (L and w are NOT defaulted here, on purpose)
%     An explicitly supplied P.L (half-length) or P.w (tip-to-tip width,
%     = 2L) wins over the L(a) rule; the rule is used only when neither is
%     given. See ribLength.m.

if nargin < 1 || isempty(P)
    P = struct();
end

defaults = {
    'celltype',        'rib'
    'xsect',           'rect'
    'beamMat',         'diamond'
    'anisoMat',        1
    'rxtal',           0
    'a',               300e-9
    't',               100e-9
    's',               100e-9
    'L_ref',           400e-9
    'a_ref',           240e-9
    'dLda',            -0.4998
    'd',               0
    'th',              110e-9
    'ncellHalf',       14
    'ribatctr',        1
    'ndef',            5
    'maxdef',          0.2
    'taperFunc',       'quadratic'
    'taperConvention', 'centerspacing'
    'airrad',          3.1e-6
    'airZOffset',      60e-9
    'mEndBC',          'lowreflecting'
    'lbd0',            1.55e-6
    'holeatedge',      0
    'plotgeom',        0
    };
for k = 1:size(defaults,1)
    if ~isfield(P, defaults{k,1})
        P.(defaults{k,1}) = defaults{k,2};
    end
end
% These depend on fields filled above.
if ~isfield(P,'airradUC')
    P.airradUC = 5*P.a;
end
filletDefault = 40e-9*(P.d == 0);
if ~isfield(P,'r1')
    P.r1 = filletDefault;
end
if ~isfield(P,'r2')
    P.r2 = filletDefault;
end
end
