% Front panel script for a nanobeam OMC cavity terminated by a 2D cross-cell
% phononic shield, with a wrap-around mechanical PML for the radiative Q.
%
% Pipeline:  RunNanobeamFEM -> CreateNanobeamGeom
%                           -> BuildNanobeamCrossShieldFEM   (P.celltype)
%                           -> SetupNanobeamFEM              (reused unchanged)
%                           -> SolveNanobeamFEM              (reused unchanged)
%
% BRING-UP ORDER - do not skip straight to a solve. See the staged block at the
% bottom of this file; S0-S3 need no solver and S1 needs only the geometry.

clear all; close all; clc                                   %#ok<CLALL>

%% geometry parameters - beam and holes
P.xsect = 'rect';                       % v1 of the cross-shield builder is rect only
P.celltype = 'crossShield';             % dispatch key in RunNanobeamFEM
P.beamMat = 'diamond';
P.anisoMat = 1;

P.a = 529e-9;                           % nominal HOLE lattice constant
P.w = 750e-9;                           % beam width
P.theta = 45;                           % etch angle (no effect for rect)
P.th = 500e-9;                          % slab thickness - beam AND shield
P.hx = 197e-9;                          % nominal hole height (along x)
P.hy = 578e-9;                          % nominal hole width  (along y)

P.nholes = 18;                          % # holes in 1/2 beam length
P.ndef = 6;                             % # holes in 1/2 defect region
P.maxdef = 0.1926;                      % defect percentage
P.oblong = 2.276;                       % oblong parameter
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
% THESE ARE PLACEHOLDER CELL DIMENSIONS. Before believing any Q, run stage S8:
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
P.aShield = 1600e-9;                    % cross lattice pitch, BOTH x and y
P.hShield = 1500e-9;                    % cross arm length (must be < aShield)
P.wShield =  800e-9;                    % cross arm width

% Tight-pitch alternative - expect the beamOverhangsNode warning:
% P.aShield = 900e-9;  P.hShield = 800e-9;  P.wShield = 400e-9;

P.nShieldX = 3;                         % cells along x
P.nShieldY = 2;                         % cells along y IN THE SIMULATED HALF
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
P.solveOpt = 0;                         % optics + shield is not supported in v1
P.calcG = 1*(P.solveMech && P.solveOpt);
P.calcS = 0*P.solveMech;
P.solveMechPML = 1;                     % 1 to extract the radiative mechanical Q

P.plotgeom = 1;                         % mphgeom snapshot into datLoc
P.storeMPH = 1;                         % keep the .mph - you will want to look
P.plotMech = 1*P.solveMech;
P.plotOpt = 1*P.solveOpt;
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

%% optical simulation parameters (unused while solveOpt = 0)
P.oevenx = (-1)^(P.holeatctr);
P.oeveny = -1;
P.oevenz = 1;
P.oneigs = 1;
P.oMesh = 3;
P.oAdjMesh = 1;

%% SiV strain coupling (unused while calcS = 0)
P.zSiV = {[1 1 1]};
P.xSlc = 0;
P.ySlc = 0;
P.zSlc = 0;

%% mechanical PML
% lambda_shear = sqrt(c44/rho)/f = sqrt(578e9/3500)/7e9 = 1.84 um for diamond
% at 7 GHz, so 2 um is about one wavelength. Confirm with the S6 sweep.
P.PMLLen  = 2e-6;                       % +x frame thickness
P.PMLLenY = 2e-6;                       % +y frame thickness
P.PMLmeshDiv = 8;                       % elements across the absorbing direction
P.PMLScalingType = 'rational';          % wavelength independent - right for an
                                        % eigenfrequency study. 'userDefined'
                                        % reproduces the legacy PML node.

%% simulation settings
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
%
% S1  geometry only. Comment out the RunNanobeamFEM call above, set
%     P.nShieldX = P.nShieldY = 2 and P.solveMechPML = 0, then run the block
%     below. Check in the GUI that the beam unions cleanly onto the shield at
%     x = beamLenHalf and that no selection warned about being empty.
%
% S2  fillets. P.shieldFillet = 'analytic' with P.r1Shield = P.r2Shield = 10e-9
%     should build. 'legacy' is EXPECTED to warn and probably fail for
%     nx,ny > 1 - run it once so you have seen the failure mode.
%
% S3  mesh only. Re-enable the PML, then in the block below add
%       mesh = model.mesh.create('mesh', P.geomname);
%       mesh.feature('size').set('hauto', P.mMesh); mesh.run;
%       disp(mphxmeshinfo(model).nelements)
%     and tune P.shieldHmax / P.PMLmeshDiv until the count is affordable.
%
% S4  eigenfrequency with P.solveMechPML = 0: real frequencies, cavity mode
%     near P.freq survives, shield modes fail the locRatio filter.
% S5  P.solveMechPML = 1: complex eigenvalues, real(lambda) > 0, finite QAll.
% S6  Q vs P.PMLLen in {1, 2, 3, 4} um and Q vs P.PMLmeshDiv in {4, 8, 16}.
% S7  Q vs P.nShieldX = 1..5 - should climb roughly exponentially inside the
%     bandgap. Immediate saturation means the mode is not in the gap, or a BC
%     is wrong. This is the single most diagnostic test.
% S8  independent bandgap check via buildCrossUnitCell + runBands at th = P.th.
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
%
% mphlaunch(model);
