% Front panel script for a nanobeam OMC cavity terminated by a 2D cross-cell
% phononic shield, with a wrap-around mechanical PML for the radiative Q, and
% with the optical mode solved in the SAME COMSOL model so that the
% optomechanical coupling g_OM can be computed.
%
% Pipeline:  RunNanobeamFEM -> CreateNanobeamGeom
%                           -> BuildNanobeamCrossShieldFEM   (P.celltype)
%                           -> SetupNanobeamFEM              (reused unchanged)
%                           -> SolveNanobeamFEM              (reused unchanged)
%                           -> CalcGOM                       (needs BOTH studies
%                                                             on one model)
%
% WHY ONE MODEL: CalcGOM.m:25-28 requires ds.ofem AND ds.mfem, and
% CalcGOM.m:124-155 Joins the optical and mechanical datasets, so the two modes
% have to come out of a single geometry. That is the whole reason the builder
% had to learn optics instead of the optical run being a separate script.
%
% THE AIR CYLINDER IS SCOPED TO THE CAVITY. The full-length air cylinder of
% BuildNanobeamFEM would enclose the shield and the PML frame and destroy the
% phononic bandgap. Here it stops at P.airCylLen, one free-space wavelength
% short of the shield front face, with 15 of the 18 half-beam holes inside it.
% Read WHY THE AIR CYLINDER STOPS SHORT in BuildNanobeamCrossShieldFEM.m before
% changing P.airrad, P.airCylLen, P.nholes or P.shieldPadLen - the builder
% hard-errors on an overlap, but the convergence question (stage S4b) is yours.
%
% BRING-UP ORDER - do not skip straight to a solve. See the staged block at the
% bottom of this file; S0-S3 need no solver and S1 needs only the geometry.

clear all; close all; clc                                   %#ok<CLALL>

%% geometry parameters - beam and holes
P.xsect = 'rect';                       % v1 of the cross-shield builder is rect only
P.celltype = 'crossShield';             % dispatch key in RunNanobeamFEM
P.beamMat = 'diamond';
P.anisoMat = 1;

P.a = 650e-9;                           % nominal HOLE lattice constant
P.w = 800e-9;                           % beam width
P.theta = 45;                           % etch angle (no effect for rect)
P.th = 250e-9;                          % slab thickness - beam AND shield
P.hx = 343e-9;                          % nominal hole height (along x)
P.hy = 617e-9;                          % nominal hole width  (along y)

P.nholes = 18;                          % # holes in 1/2 beam length
P.ndef = 6;                             % # holes in 1/2 defect region
P.maxdef = 0.16;                      % defect percentage
P.oblong = 1.15;                       % oblong parameter
P.aspect_ratio = 1;

% cavity taper params
P.holeatctr = 1;                        % 1/0 for hole/dielectric in middle
P.taperFunc = 'cubic';

% end waveguide mirror taper params
P.wvgmir = 0;                           % the cross shield is the termination
P.wgmTaper.func = 'cubic';
P.wgmTaper.endtype = 'custom';
P.wgmTaper.a_end = 650e-9;
P.wgmTaper.hx_end = 343e-9;
P.wgmTaper.hy_end = 300e-9;

P.asymCav = 0;                          % not supported by this builder in v1

P.lambda = 1550e-9;
P.nbeam = 2.386;

% Disorder
P.stdDev = [0,0];
P.stdDevPos = 0;
P.asym = 0;
P.extractLocMechModes = 1;

%% cross-cell phononic shield
% ---------------------------------------------------------------------------
% THESE ARE PLACEHOLDER CELL DIMENSIONS. Before believing any Q, run stage S10:
% feed the SAME aShield/hShield/wShield and th = P.th through the band
% structure pipeline (buildCrossUnitCell + runBands) and confirm that P.freq
% actually lands inside the cross lattice's bandgap. A shield outside its gap
% gives a Q that does not grow with nShieldX, which is the tell.
%
% Why these numbers:
%   node   = aShield - wShield = 800 nm  >= P.w = 750 nm, so the beam end face
%            lands entirely on a node and the clamp audit stays quiet even at
%            P.shieldPadLen = 0.
%   wall   = aShield - hShield = 100 nm, above a 50 nm lithography limit.
%   solid fraction 0.31, so there is real mass contrast.
%
% The tension to be aware of: a bandgap at 7 GHz in diamond wants a pitch near
% v_shear/(2f) = 0.92 um, and 1.6 um is not that. The tight-pitch alternative
% below is the other horn of the dilemma - it has the right pitch but a 500 nm
% node against a 750 nm beam, so 33% of the beam end face lands on a 50 nm fin
% and the builder says so. Pick one deliberately.
% ---------------------------------------------------------------------------
P.aShield = 894e-9;                    % cross lattice pitch, BOTH x and y
P.hShield = 863e-9;                    % cross arm length (must be < aShield)
P.wShield =  416e-9;                    % cross arm width

% Tight-pitch alternative - expect the beamOverhangsNode warning:
% P.aShield = 900e-9;  P.hShield = 800e-9;  P.wShield = 400e-9;

P.nShieldX = 6;                         % cells along x
P.nShieldY = 5;                         % cells along y IN THE SIMULATED HALF
                                        % -> 4 physical rows; y = 0 bisects a
                                        %    row of nodes. Filename records
                                        %    both as nsyH2_nsyP4.

P.shieldPadLen = 0;                     % 0 = lattice faithful (default).
                                        % ~aShield/2 gives a proper clamp and
                                        % bonds away the free-edge fins.
P.shieldFillet = 'none';                % 'none' | 'analytic' | 'legacy'.
P.r1Shield = 0;                         % ignored while shieldFillet = 'none'
P.r2Shield = 0;
P.shieldMinFeature = 50e-9;             % for the isCrossFabricable audit
% P.shieldHmax = 40e-9;                 % default is min((a-h)/3, th/3)

%% simulation / calculation / plot / save options
P.solveMech = 1;
P.solveOpt = 1;                         % air cylinder is cavity-scoped - see
                                        % the header and P.airCylLen below
P.calcG = 1*(P.solveMech && P.solveOpt);    % -> 1
P.calcS = 0*P.solveMech;
P.solveMechPML = 1;                     % 1 to extract the radiative mechanical Q

P.plotgeom = 1;                         % mphgeom snapshot into datLoc
P.storeMPH = 1;                         % keep the .mph - you will want to look
P.plotMech = 1*P.solveMech;
P.plotOpt = 1*P.solveOpt;               % -> 1, follows solveOpt
P.plotStrCpl = 1*P.calcS;

%% mechanical simulation parameters
P.mevenx = 1;                           % +/-1 even/odd about x
P.meveny = 1;                           % +/-1 required (0 is rejected)
P.mevenz = 1;                           % +/-1 even/odd about z
P.freq = 7e9;                           % target mechanical frequency
P.mneigs = 30;                          % more than the no-PML case: a PML run
                                        % also returns PML-localised modes
P.mMesh = 3;
P.mAdjMesh = 1;

P.rxtal = 0;
P.rxtalInFilename = 1;

%% optical simulation parameters
% Symmetries are the rect-beam fundamental, identical to
% test_nanobeamRectFEM.m:99-102, the closest sibling:
%   oevenx = (-1)^holeatctr = -1 for a hole at the centre. SetupNanobeamFEM.m:39
%           then raises it to the power holeatctr again, so the value that
%           reaches the BC is (-1)^1 = -1 -> PEC on the x = 0 plane.
%   oeveny = -1  TE-like (E mostly along y, so E_y is odd about y = 0 -> PEC).
%   oevenz = +1  even about the slab mid-plane -> PMC on z = 0. This MUST have
%           the same |value| as P.mevenz (both 1 here): there is one z = 0 cut
%           in the geometry and the builder hard-errors if the two physics
%           disagree about whether to make it.
P.oevenx = (-1)^(P.holeatctr);
P.oeveny = -1;
P.oevenz = 1;

% oneigs = 1 matches every rect sibling (test_nanobeamRectFEM.m:102,
% test_nanobeamRectFEM_withPML.m:98, test_nanobeamTri*.m, and
% test_fabricatednanobeamRectFEM.m:137); only test_1DSnowflakeFEM.m:109 asks
% for 5. Keep 1 here: the shift-invert target is set to the exact design
% wavelength (SolveNanobeamFEM.m:64 shifts to -1i*2*pi*c/P.lambda), and
% RunNanobeamFEM.m:163 passes oModes = 1 to CalcGOM unconditionally - so
% asking for more eigenvalues would not change which one the coupling is
% computed for. If S4b shows the solver landing on a band-edge mode instead of
% the fundamental, raise this to 5 and re-check which solution index is the
% high-Q one (SolveNanobeamFEM.m:199 prints ofem.hiQSol).
P.oneigs = 1;
P.oMesh = 3;
P.oAdjMesh = 1;

% Air half-cylinder around the CAVITY. 2*lambda + w/2 is the convention in
% test_nanobeamRectFEM.m:105 and five sibling scripts.
P.airrad = 2*P.lambda + P.w/2;          % 2*1550 + 375 = 3475 nm
%
% P.airCylLen is the +x extent of that cylinder and is the field that makes
% optics compatible with the shield at all. Leave it unset to take the builder
% default, min(beamLenHalf, xS0 - P.lambda):
%
%   beamLenHalf = 8951.8 nm (CreateNanobeamGeom, nholes = 18, a ~ 529 nm)
%   shieldPadLen = 0   ->  xS0 = 8951.8 nm
%   airCylLen   = 8951.8 - 1550 = 7401.8 nm
%   x-gap to the shield front face   = 1550.0 nm  (1.00 lambda)
%   x-gap to the +x PML arm          = 11150.0 nm
%   x-gap to the +y PML arm          = 1550.0 nm, y-gap = 4525.0 nm
%   cap lands between hole 15 (ends 7198.8) and hole 16 (starts 7530.8), so
%   15 of 18 half-beam holes - 6 defect + 9 mirror periods - are in the air.
%
% The cylinder reaches |z| = 3475 nm against a 500 nm slab and y = 3475 nm
% against an 8000 nm tall shield, so x is the ONLY separation. Shorten the beam
% (P.nholes) or lengthen the pad (P.shieldPadLen) and this number moves with it.
% P.airCylLen = 7.4e-6;                 % uncomment to override the default

%% SiV strain coupling (unused while calcS = 0)
P.zSiV = {[1 1 1]};
P.xSlc = 0;
P.ySlc = 0;
P.zSlc = 0;

%% mechanical PML
% lambda_shear = sqrt(c44/rho)/f = sqrt(578e9/3500)/7e9 = 1.84 um for diamond
% at 7 GHz, so 2 um is about one wavelength. Confirm with the S8 sweep.
P.PMLLen  = 2e-6;                       % +x frame thickness
P.PMLLenY = 2e-6;                       % +y frame thickness
P.PMLmeshDiv = 8;                       % elements across the absorbing direction
P.PMLScalingType = 'userDefined';       % FREQUENCY DEPENDENT: the PML stretch
                                        % is keyed to P.freq via
                                        % typicalWavelength = v/P.freq.
                                        % Set 'rational' for a wavelength
                                        % independent cross-check in S8.
% Reference speed for typicalWavelength. Unset = c11 (longitudinal, 17.5 km/s,
% lambda = 2.50 um at 7 GHz), which is the legacy choice. Uncomment for the
% shear branch (12.85 km/s, 1.84 um), which is what actually carries radiation
% out of a thin suspended slab.
% P.PMLWaveSpeed = sqrt(P.D(10)/P.rho);

%% simulation settings
% DOF BUDGET FOR THE COMBINED RUN - READ BEFORE THE FIRST SOLVE.
% SolveNanobeamFEM makes ONE mesh node and both studies share it, but the two
% physics are selected on different domains (emw on beam + shield + air, smech
% on beam + shield + PML), so they do NOT share degrees of freedom and adding
% optics does NOT inflate the mechanical count.
%
% Estimates at these parameters (the builder prints its own version of this
% line; both use the crude V_tet = h^3/6, which understates a real tet fill by
% roughly 1.4x):
%   optical     hmax = lambda/5 = 310 nm over 84.9 um^3 of quarter cylinder
%               -> 1.7e4 tets, ~3.4e5 DOF.  Comfortable.
%   mechanical  shield hmax = min((a-h)/3, th/3) = 33.3 nm over 6.0 um^3 of
%               solid shield -> 9.7e5 tets, plus 5.3e3 in the PML at
%               PMLLen/PMLmeshDiv = 250 nm -> ~4.1e6 DOF, and nearer 5.8e6
%               once the h^3/6 optimism is removed.  THAT is the constraint,
%               and it is already the constraint with P.solveOpt = 0.
%
% If the mechanical side runs out of budget, mAdjMesh cannot rescue it: the
% loop raises hauto, but applyMeshOverrides re-applies hmax = P.shieldHmax to
% the shield on every pass, so the dominant term never moves. Levers, in order:
%   1. P.nShieldX = 3, P.nShieldY = 3 for the combined run (the shield tet
%      count scales with nShieldX*nShieldY), then take Q from a mech-only run
%      at the full 6x5 and g_OM from the small one - the optical mode does not
%      touch the shield.
%   2. P.shieldHmax = 45e-9 (just over two elements across the 100 nm
%      ligament) cuts the shield tets by 2.4x. Do not go to exactly
%      (a-h)/2 = 50e-9: the builder's shieldMeshCoarse warning fires there,
%      and two elements across the feature that carries the entire shield
%      response is already the floor. Check the Q against 33.3 nm either way.
%   3. Raise P.max_dof to 8e6 if the machine has the memory.
% A fourth fix - separate optical and mechanical mesh nodes - would need an
% edit to SolveNanobeamFEM, which this script deliberately does not require.
P.max_dof = 5e6;

%% single run
currentDate = datestr(now,'mmddyyyy');                      %#ok<TNOW1,DATST>
datLoc = [fullfile('.','test','1D_OMC_crossShield',currentDate),filesep];
[ds,model] = RunNanobeamFEM(P,datLoc);

%% ========================================================================
%  STAGED BRING-UP
%  ========================================================================
%
% S0  static analysis, no COMSOL:
%       checkcode('BuildNanobeamCrossShieldFEM.m')
%       checkcode('test_nanobeamCrossShieldFEM.m')
%     Both should be silent.
%
% S0b geometry arithmetic, no COMSOL. P = LoadMaterialParams(P);
%     P = CreateNanobeamGeom(P); then check the numbers the air cylinder
%     depends on:
%       P.beamLenHalf                       expect 8951.8 nm
%       2*P.lambda + P.w/2                   expect 3475.0 nm
%       P.beamLenHalf - P.lambda             expect 7401.8 nm  (= airCylLen)
%       max(P.geomHalf(:,2))/2 vs P.w/2      expect 289 vs 375 nm
%     The last one is load bearing: the builder separates the air's z = 0
%     faces from the beam's by a y bound, and needs max(hy)/2 < w/2.
%
% S1  geometry only, MECHANICS ONLY. Comment out the RunNanobeamFEM call above,
%     set P.solveOpt = 0, P.nShieldX = P.nShieldY = 2 and P.solveMechPML = 0,
%     then run the block below. Check in the GUI that the beam unions cleanly
%     onto the shield at x = beamLenHalf and that no selection warned about
%     being empty. This is the unchanged v1 path and must look exactly as it
%     did before optics existed.
%
% S1b geometry only, OPTICS ON. Set P.solveOpt = 1 and repeat. What to check in
%     the GUI, in order of how badly it bites:
%       - the air cylinder is a QUARTER cylinder (y >= 0 and z >= 0). If the
%         z < 0 half is still there, the symZ cut block did not span far
%         enough in z and the z = 0 optical mirror plane is wrong.
%       - there is bare vacuum, not air, between x = P.airCylLen and the
%         shield front face.
%       - the beam is TWO domains, split at x = P.airCylLen by the end-cap
%         imprint. That is expected; see the builder header.
%       - the console lines "Air cylinder clears the cross shield: x by
%         1550.0 nm", "Air end cap ... N interior (removed)" and
%         "z = 0 plane: ... air (optical only) ... solid (mechanical)".
%         A 0 in any of those is a real failure, not cosmetics.
%       - disp(P.domSel) and disp(P.bndSel): cyl, cylXsym, cylXend, cylYsym,
%         cylZsym, cylCurv all present and non-empty.
%
% S2  fillets. P.shieldFillet = 'analytic' with P.r1Shield = P.r2Shield = 10e-9
%     should build. 'legacy' is EXPECTED to warn and probably fail for
%     nx,ny > 1 - run it once so you have seen the failure mode.
%
% S3  mesh only. Re-enable the PML, then in the block below add
%       mesh = model.mesh.create('mesh', P.geomname);
%       mesh.feature('size').set('hauto', P.mMesh); mesh.run;
%       disp(mphxmeshinfo(model).nelements)
%     and tune P.shieldHmax / P.PMLmeshDiv until the count is affordable. See
%     the DOF BUDGET note next to P.max_dof; at 6x5 cells and 33.3 nm this is
%     expected to be tight or over, with or without optics.
%
% S4  OPTICS ONLY bring-up, and do this before the combined solve. Set
%     P.solveMech = 0... except the builder requires P.solveMech = 1, so
%     instead shrink the mechanics out of the way:
%       P.nShieldX = 1; P.nShieldY = 1; P.solveMechPML = 0; P.mneigs = 3;
%     A 1x1 shield still terminates the beam and still exercises the union,
%     but costs ~3e4 shield tets instead of ~1e6, so the optical study is the
%     only expensive thing in the model. Expect from SolveNanobeamFEM.m:210 a
%     single high-Q mode with lambda within a few tens of nm of 1550 and
%     Q > 1e4; PlotEy should show a mode centred on x = 0 and decaying into
%     the mirror. If lambda is far from 1550 the taper is off, not the
%     cylinder.
%
% S4b optical convergence against the truncated air. Sweep
%     P.airCylLen = [5.5, 6.5, 7.4] um and watch lambda, Q and (once S7 runs)
%     g_OM. lambda and g_OM should be flat - they are set by the field at
%     x ~ 0. Q will NOT be flat and will be OVERESTIMATED at short airCylLen,
%     because the beam beyond the cylinder is clad in the EM interface's
%     default PEC, which is a perfect mirror where the device has vacuum. If
%     lambda or g_OM move by more than a per cent, the air is too short and
%     the cavity is being squeezed - the honest fix is fewer mirror holes
%     between the defect and the shield, not a longer cylinder.
%
% S5  eigenfrequency with P.solveMechPML = 0: real frequencies, cavity mode
%     near P.freq survives, shield modes fail the locRatio filter.
% S6  P.solveMechPML = 1: complex eigenvalues, real(lambda) > 0, finite QAll.
% S7  combined solve with P.calcG = 1. CalcGOM Joins the optical and
%     mechanical datasets, so this is the first stage that can fail on the
%     shared mesh node rather than on physics: the optical study solves with
%     hmax = lambda/5 applied globally and no shield refinement, and the
%     mechanical study then re-meshes with hmax = P.shieldHmax on the shield.
%     Sanity checks on the result: ds.cpl.gMax should be of order 100 kHz for
%     a diamond OMC at 7 GHz, and ds.cpl.optWvl should equal the lambda
%     reported in S4. A g_OM that is orders of magnitude off, or a changed
%     optical wavelength, points at the mesh and not at the geometry.
% S8  Q vs P.PMLLen in {1, 2, 3, 4} um and Q vs P.PMLmeshDiv in {4, 8, 16}.
% S9  Q vs P.nShieldX = 1..5 - should climb roughly exponentially inside the
%     bandgap. Immediate saturation means the mode is not in the gap, or a BC
%     is wrong. This is the single most diagnostic test. g_OM must NOT change
%     with nShieldX; if it does, the optical mode is feeling the shield and
%     P.airCylLen is too long.
% S10 independent bandgap check via buildCrossUnitCell + runBands at th = P.th.
%
%% test the model building function (geometry only - no solve)
% import com.comsol.model.*
% import com.comsol.model.util.*
% 
% ModelUtil.showProgress(true);
% ModelUtil.clear();
% clear ds model
% ds = []; model = [];
% 
% model = ModelUtil.create('model');
% P = LoadMaterialParams(P);
% P = CreateNanobeamGeom(P);
% [model,P] = BuildNanobeamCrossShieldFEM(model,P);
% 
% disp('--- domain selections ---'); disp(P.domSel);
% disp('--- boundary selections ---'); disp(P.bndSel);
% disp('--- shield extents (m) ---'); disp(P.shield);
% if isfield(P,'airCyl'), disp('--- air cylinder (m) ---'); disp(P.airCyl); end
%
% mphlaunch(model);
