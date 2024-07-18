clear all; clc; close all;
clear P;

%% unit cell params 
P.xsect = 'rect'; 
P.beamMat = 'diamond';                  % beam material name
P.celltype = 'cross';                   % specify the cell type
P.a = 450e-9;              % lattice constant (along x)
P.hc = 378e-9;              % unit cell width (along y)
P.wc = 95e-9;              % unit cell thickness (along z)
P.th = 160e-9;             % height (along x) of cross (for celltype = 'hollow')
                            % or of inner block (for celltype = 'solid')
P.r1 = 30e-9;             % width (along y) of cross (for celltype = 'hollow')
                            % or of inner block (for celltype = 'solid')
P.r2 = 30e-9;              % height (along x) of each leg in cross (for celltype = 'hollow')
                            % or of outer fins (for celltype = 'solid')
P.nperiod = 1;  % no. of periods to simulate for
P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
P.mbevenz = 1;      % 1 to find even mechanical mode about z

P.kpts = 3;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 9;                           % no. of bands to solve for

P.solveasym = 1;                        % 1 to solve for antisymmetric bands
P.completeBandGaps = 1;                 % 1 to plot complete bandgaps (across all symmetries)
P.plotgeom = 1;                         % 1 to plot the geometry
P.savedat = 1;                          % 1 to save data structures
P.savebndplot = 1;                      % 1 to save bandstructure plot
P.saveplots = 1;                        % 1 to save displacement and strain profiles
P.saveMPH = 0; 
P.bandStruct_2D = 1;                 % 1 to simulate 2D band structures

%% mechanical simulation parameters 
% solid mechanics solver parameters
P.mbeveny = 1;                          % 1 to find even mechanical mode about y
P.mbevenz = 1;                          % 1 to find even mechanical mode about z
P.freq = 0;                             % target frequency - set to 0 for bandstructure simulations
P.meshSize = 3;                         % mesh quality for mechanical simulations
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
% DrawCrossUnitCell(model,P);
% mphlaunch(model);
%% Single solve
datLoc = '.\test\cross\071424\';
P.datLoc = datLoc;
bds = solveBands(P);