clear all; clc; close all;
clear P;

%% unit cell params 
P.xsect = 'rect'; 
P.beamMat = 'diamond';                  % beam material name
P.celltype = 'boomerang_strip';                   % specify the cell type
P.unitcell = 'hexagonal';                  % specify the shape of the unit cell
P.a = 480e-9;              % lattice constant 
P.w = 140e-9;              % unit cell width (along x)
P.r = 177e-9;              % unit cell height (along y)
P.th = 220e-9;             % height (along x) of cross (for celltype = 'hollow')
                            % or of inner block (for celltype = 'solid')
% for the center unit cells v2 
P.wo = 372e-9;
P.wi = 91e-9;
P.ho = 210.7e-9;
P.hi = 196e-9;
P.d = sqrt(3)*P.a/2;
% center unit cell parameters

% for the center unit cells v1 
P.h = 210.7e-9;
P.d1 = 140.7e-9;

P.r1 = 10e-9;             % width (along y) of cross (for celltype = 'hollow')
                            % or of inner block (for celltype = 'solid')
P.r2 = 10e-9;              % height (along x) of each leg in cross (for celltype = 'hollow')
                            % or of outer fins (for celltype = 'solid')
P.nperiod = 1;  % no. of periods to simulate for
P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
P.mbevenz = 1;      % 1 to find even mechanical mode about z

P.kpts = 10;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 20;                           % no. of bands to solve for

P.solveasym = 1;                        % 1 to solve for antisymmetric bands
P.completeBandGaps = 1;                 % 1 to plot complete bandgaps (across all symmetries)
P.plotgeom = 1;                         % 1 to plot the geometry
P.savedat = 1;                          % 1 to save data structures
P.savebndplot = 1;                      % 1 to save bandstructure plot
P.saveplots = 0;                        % 1 to save displacement and strain profiles
P.saveMPH = 0; 
P.bandStruct_2D = 0;                 % 1 to simulate 2D band structures
% the symmetry condition parameters 
P.TwoSymPlanes = 1; % if we are solving for band structures with two symmetry planes
P.zSymCondition = 1;
% for optical simulation
P.bandStructureDim=1;           % specify the dimension of the band structure 
P.optical_freq = 100;       % specify the target frequency (THz)

% for the optical bandgap 
P.add_airDisk = 1;
P.airDiskH = 5000e-9;
P.mbevenz = 0;

%% mechanical simulation parameters 
% solid mechanics solver parameters
P.mbeveny = 0;                          % 1 to find even mechanical mode about y
% P.mbevenz = 1;                          % 1 to find even mechanical mode about z
P.freq = 15e9;                             % target frequency - set to 0 for bandstructure simulations
P.meshSize = 4;                         % mesh quality for mechanical simulations
P.fixed_bc = 0;                       % 1 to fixed the boundaries for xz planes at y = +/- w/2

P.anisoMat = 1;
P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                        % from <100> inplane direction about <100> surface normal

%% define the maximum number of degree of freedom to limit the simulation time
P.max_dof = 3e6;                        % max # of degrees of freedom

% %% debugging the unit cells 
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
% buildBoomerangStrip_3D(model,P);
% mphlaunch(model);
%% Single solve
currentDate = datestr(now,'mmddyyyy');
datLoc = [fullfile('.','test','boomerang_strip',currentDate),filesep];
P.datLoc = datLoc;
% optical band

bds = solveOpticalBands(P);
% % mechanical band 
% P.add_airDisk = 0;
% bds = solveBands(P);