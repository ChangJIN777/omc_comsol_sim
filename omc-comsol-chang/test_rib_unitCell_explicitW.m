clear all; clc; close all;
clear P;

%% rib unit cell with an explicit beam width (P.w) and fillets
% Drives celltype 'rib' (builder buildRibUnitCell_LN.m: spine + centred rib)
% with the cell dimensions of the original LN rib design, but in diamond,
% with vertical walls (d = 0) and 20 nm fillets. Unlike test_rib_unitCell.m
% (visible design, L(a) rule), this script fixes the rib length through an
% EXPLICIT P.w, which wins over the L(a) rule.
%   P.run_optical = 0  mechanical bands: solveBands (symmetry sectors)
%   P.run_optical = 1  optical bands:    solveOpticalBands -> runOpticalBand_1D
% Formerly test_LN_unitCell.m. Full parameter reference: rib_README.md.
P.xsect = 'rect'; 
P.beamMat = 'diamond';                  % material, set explicitly (also the 'rib' default)
P.celltype = 'rib';                   % merged rib device
P.unitcell = 'rectrangular';                  % read by the optical runners ('hexagonal' or not)
P.a = 550e-9;              % lattice constant 
P.s = 100e-9;              % unit cell spine width 
P.w = 1300e-9;              % EXPLICIT tip-to-tip width (= 2L): wins over the L(a) rule
P.th = 250e-9;             % slab thickness
P.t = 250e-9;       % the rib width 
P.d = 0*pi/180;            % sidewall angle in RADIANS. d ~= 0 needs P.r1 = P.r2 = 0
                           % and mbevenz = 0 (keep zSymCondition = 0 below)

P.r1 = 40e-9;   % concave rib/spine (root) fillet; allowed only with d = 0
P.r2 = 40e-9;   % convex rib-tip fillet;            allowed only with d = 0

P.nperiod = 1;  % no. of periods to simulate for
P.holeatedge = 0;   % read by runBands, unused by this builder
P.mbevenz = 0;      % 1 to find even mechanical mode about z

P.kpts = 10;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 10;                           % no. of bands to solve for

P.solveasym = 0;                        % 1 to solve for antisymmetric bands
P.completeBandGaps = 1;                 % 1 to plot complete bandgaps (across all symmetries)
P.plotgeom = 1;                         % 1 to plot the geometry
P.savedat = 1;                          % 1 to save data structures
P.savebndplot = 1;                      % 1 to save bandstructure plot
P.saveplots = 1;                        % 1 to save displacement and strain profiles
P.saveMPH = 0; 
P.bandStruct_2D = 0;                 % 1 to simulate 2D band structures
P.bandStructureDim = 1;     % Optical band only: 3 to simulate 3d structures 

% for choosing between the optical and the mechanical band structure
% (1 adds the air half-cylinder; runBands refuses it, so mechanics needs 0)
P.run_optical = 0;

% if we are going to save the raw data files 
P.saveRawData = 1;

%% mechanical simulation parameters 
% solid mechanics solver parameters
P.TwoSymPlanes = 1; 
% TODO(confirm): zSymCondition was missing, so solveBands stopped with
% 'Unrecognized field'. The original script (HEAD) solved mechanics with
% solveBands_noSym, i.e. no symmetry plane, and set mbevenz = 0. solveBands
% has no symmetry-free mode, so 0 is chosen: a y-symmetric solve
% (mbeveny = +1, and -1 too if solveasym = 1) that leaves mbevenz = 0,
% which also stays valid if d is made nonzero. Use 1 for z-symmetry
% sectors (d = 0 only), or solveBands_noSym to match the original exactly.
P.zSymCondition = 0;
P.mbeveny = 0;                          % 1 to find even mechanical mode about y
P.mbevenz = 0;                          % 1 to find even mechanical mode about z
P.freq = 3e9;                             % Mechanical target frequency - set to 0 for bandstructure simulations
P.optical_freq = 200;           % optical target frequency - in THz 
P.meshSize = 4;                         % mesh quality for mechanical simulations
P.fixed_bc = 0;                       % 1 to fixed the boundaries for xz planes at y = +/- w/2

P.anisoMat = 1;
P.rxtal = 0;                            % ccw rotation of elasticity matrix in deg 
                                        % from <100> inplane direction about <100> surface normal

%% define the maximum number of degree of freedom to limit the simulation time
P.max_dof = 3e6;                        % max # of degrees of freedom

%% debugging the unit cells 
% % import COMSOL class
% import com.comsol.model.*
% import com.comsol.model.util.*
% 
% ModelUtil.showProgress(true);
% ModelUtil.clear();
% clear model
% 
% % create COMSOL model named 'model' from which COMSOL methods can be called, 
% % e.g. model.save
% model = ModelUtil.create('model');
% 
% [model,P] = buildRibUnitCell_LN(model,P);
% mphlaunch(model);
%% Single solve
if P.run_optical
    currentDate = datestr(now,'mmddyyyy');
    datLoc = [fullfile('.','test','ribUnitCell_optical',currentDate),filesep];
    P.datLoc = datLoc;
    bds = solveOpticalBands(P);
else
    currentDate = datestr(now,'mmddyyyy');
    datLoc = [fullfile('.','test','ribUnitCell',currentDate),filesep];
    P.datLoc = datLoc;
    bds = solveBands(P);
end