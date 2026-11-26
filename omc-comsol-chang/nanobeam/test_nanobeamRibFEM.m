% front panel script for nanobeam FEM simulations of the rib cavity (default:
% the visible design of comsol_templates/visible_rib_omc.m), built by createRibCavityGeom +
% buildRibCavity and run through RunNanobeamFEM / SetupNanobeamFEM /
% SolveNanobeamFEM like every other celltype.

clear all; close all; clc

%% path: the shared rib helpers live in the omc-comsol-chang root
% (ribDefaultParams, validateRibParams, ribLength,
% addRibCell). Added at the END of the path so this folder's own
% CreateFileBase.m / LoadMaterialParams.m / RotateXtalTensor.m still win.
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir)
    thisDir = pwd;
end
addpath(fileparts(thisDir), '-end');

%% geometry parameters (defaults and their sources: ribDefaultParams.m)
P.xsect = 'rect';                       % vertical sidewalls, as in the template
P.celltype = 'rib';             % specify what type of unit cells we are simulating
P.beamMat = 'diamond';                  % beam material name
P.anisoMat = 1;

P.a  = 650e-9;                          % mirror lattice constant (cell width)
P.t  = 280e-9;                          % rib width along x
P.s  = 100e-9;                          % spine width along y
P.th = 250e-9;                          % EFFECTIVE thickness of the template
                                        % (z-offset defect); nominal 120e-9
P.r1 = 40e-9;                           % concave rib/spine fillet (TODO: IR-set value)
P.r2 = 40e-9;                           % convex rib-tip fillet    (TODO: IR-set value)
P.L_ref = 400e-9;                       % rib half-length at a = a_ref (template centre)
P.a_ref = 240e-9;
P.dLda  = -0.4998;                      % dL/da (template IR set); used by the L(a)
                                        % rule, and by P.w_mirror when P.w_ctr is unset

% rib tip-to-tip width w = 2L (at the bottom face; fillets and spine do not change
% it). Precedence: P.L > P.w (one constant width, no taper) > P.w_mirror (+ P.w_ctr,
% linear in a, tapered with the cavity) > the L_ref/a_ref/dLda rule above. Leave all
% unset for the rule (740 nm mirror to 800 nm centre at a = 300 nm, maxdef = 0.2).
P.w_mirror = 1300e-9;                   % w at the mirror cells (a = P.a)
% P.w_ctr  = 800e-9;                    % w at the centre cell (a = P.a*(1 - P.maxdef));
                                        % 'rib' scheme only (an error with 'hole')
% P.w = 740e-9;                         % alternative: one constant w for every cell
                                        % (remove w_mirror/w_ctr)
P.d = 0;                                % sidewall angle in RADIANS; d ~= 0 needs
                                        % r1 = r2 = 0 and mevenz = oevenz = 0

% cavity taper
P.ncellHalf = 30;                       % cells in x >= 0 half, centre counted once
P.ribatctr = 1;                         % 1/0 for rib/spine gap at the centre
P.ndef = 7;                             % taper length in cells
P.maxdef = 0.1765;                         % centre a = a*(1 - maxdef)
P.taperFunc = 'quadratic';              % linear/cubic/quadratic/none
P.taperConvention = 'centerspacing';    % 'centerspacing' (default, CreateNanobeamGeom:
                                        % x_{n+1} = x_n + a_{n+1}) or 'cellwidth' (the
                                        % template's; reproduces its positions); see README

% hole-style taper: CreateNanobeamGeom's maxdef/oblong scheme (as
% test_nanobeamRectFEM.m), with hy <-> w (tip to tip), hx <-> t, a <-> a:
%   t_k = t*sig^(1-oblong),  w_k = w_mirror*sig^(1+oblong),
%   sig = a_k/a = 1 - maxdef*f((k-1)/ndef)   (1 for mirror cells).
% A positive oblong makes the ribs shorter and wider toward the centre.
% Needs P.w_mirror (or the L(a) rule) and NO P.w_ctr / P.w / P.L.
% 'rib' (default) = t constant, w from P.w_mirror/P.w_ctr or the rule above.
% Validity: t_k < a_k, t_k > 2*r2, w_k > s + 2*(r1 + r2), checked by
% validateRibParams (createRibCavityGeom errors before COMSOL starts).
P.ribTaperScheme = 'hole';              % 'rib' | 'hole'
P.oblong = 1;                         % hole scheme only (default 2.7)

P.lambda = 1550e-9;                      % target optical wavelength (template optical
                                        % eigen-shift 400 THz)
P.nbeam = 2.386;                        % diamond refractive index (telecom value, used everywhere)

% symmetric cavity only
P.asymCav = 0;
P.asym = 0;

%% specify simulation/calculation/plot/save options
P.solveMech = 1;                        % 1 to solve for mechanics
P.solveOpt = 1;                         % 1 to solve for optics
P.calcG = 1*(P.solveMech && P.solveOpt);% 1 to calculate optomechanical coupling
P.calcS = 0*P.solveMech;                % 1 to calculate strain coupling
P.solveMechPML = 0;                     % not supported by buildRibCavity
P.extractLocMechModes = 1;              % filter modes by confinement

% plotting & saving
P.plotgeom = 0;                         % 1 to plot the geometry
P.storeMPH = 0;                         % 1 to save COMSOL model file
P.plotMech = 1*P.solveMech;             % 1 to plot displacement and strain profiles
P.plotOpt = 1*P.solveOpt;               % 1 to plot E-field profiles
P.plotStrCpl = 1*P.calcS;               % 1 to plot strain coupling profile

%% mechanical simulation parameters
P.mevenx = 1;                           % +/-1 to find even/odd mode about x; 0 for fixed BC
P.meveny = 1;                           % +/-1 to find even/odd mode about y
P.mevenz = 1;                           % +/-1 to find even/odd mode about z; 0 = full
                                        % slab height (no z cut), as the reference
                                        % optbands1D/test_nanobeamBands.m. The z cut is
                                        % shared: mevenz and oevenz must be 0 together
P.freq = 7e9;                        % target mechanical frequency (template eigen-shift)
P.mneigs = 20;                          % # of eigenvalues to find
P.mMesh = 3;                            % mesh quality for mechanical simulations
P.mAdjMesh = 1;                         % adjust mesh if DOFs exceed max_dof
P.mEndBC = 'lowreflecting';             % beam-end BC: 'lowreflecting' (template) or 'fixed'

% crystal orientation: the template's stiffness tensor is unrotated
P.rxtal = 0;
P.rxtalInFilename = 1;

%% optical simulation parameters
P.oevenx = 1;                           % forced to 1 by SetupNanobeamFEM when
                                        % P.ribatctr = 1 (dielectric at centre)
P.oeveny = -1;                          % -1 == TE-like
P.oevenz = 1;                           % +/-1 to find even/odd optical mode about z;
                                        % 0 = no z symmetry (reference obevenz = 0)
P.oneigs = 1;                           % # of eigenvalues to find
P.oMesh = 3;                            % mesh quality for optical simulations
P.oAdjMesh = 1;                         % adjust mesh if DOFs exceed max_dof
P.airrad = 3.1e-6;                      % air-cylinder radius (template)
P.airZOffset = 0;                       % cylinder axis above the beam mid-plane: 0 =
                                        % centred. 60e-9 only reproduces the th = 110 nm
                                        % template; for thicker slabs it puts the axis
                                        % inside the beam and the mesh fails

%% OM coupling parameters
P.g0min = 80e3;                         % min g0 above which to save plots for

%% SiV strain coupling
P.zSiV = {[1 1 1]};
P.xSlc = 0;                             % relative to center
P.ySlc = 0;
P.zSlc = 0;

%% Simulation settings
P.max_dof = 5e6;                        % max # of degrees of freedom

%% pre-flight (pure MATLAB): the 'cellwidth' taper must reproduce the template
% Cell centres read off the template's feature tags (finalunion_<x_nm>_<k>).
% Run explicitly with 'cellwidth', the template's convention; the default
% 'centerspacing' sits 10.5 nm further out from cell 6 on, by design. Uses a
% fresh struct, so edits above do not trip it.
Pt = createRibCavityGeom(ribDefaultParams(struct('taperConvention','cellwidth')));
xTagRibCtr = [0 241 487 743 1013 1302 1602 1902 2202 2502 2802 3102 3402 3702];
assert(isequal(round(Pt.ribGeomHalf(:,Pt.ribGeomCols.x)'*1e9), xTagRibCtr) && ...
       abs(Pt.beamLen - 7.704e-6) < 1e-12, ...
    'rib cellwidth taper no longer reproduces the template cell positions.');

% %% debugging the cavity geometry (no solve)
% import com.comsol.model.*
% import com.comsol.model.util.*
% ModelUtil.clear(); clear model
% model = ModelUtil.create('model');
% P = createRibCavityGeom(P);
% PlotDefectCells(P);
% P = LoadMaterialParams(P);
% [model,P] = buildRibCavity(model,P);
% disp(P.bndSel); disp(P.domSel);
% mphlaunch(model);

%% single run
currentDate = datestr(now,'mmddyyyy');
datLoc = [fullfile('.','test','1D_OMC_rib',currentDate),filesep];
[ds,model] = RunNanobeamFEM(P,datLoc);
