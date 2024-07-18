clear all; clc; close all;
clear P;

%% unit cell params 
P.xsect = 'rect'; 
P.beamMat = 'diamond';                  % beam material name
P.celltype = '2D_ribs';         % specify the cell type
P.a = 1500e-9;   % the lattice constant of the unit cel 
P.w = 3000e-9;  % the width of the unit cell 
P.wi = 670e-9;     % define the inner width of the rib geometry 
P.ho = 975e-9;     % define the outter height of the rib geometry 
P.hi = P.ho - 20e-9;    % define the inner height of the rib geometry
P.wo = 1000e-9;     % define the outter width of the rib geometry 
P.th = 160e-9;  % define the thickness of the unit cell 
P.ai = 150e-9;     % define the spacing between the upper and lower hole 
% P.w = P.ho+P.ai/2;  % the width of the unit cell 

P.nperiod = 1;  % no. of periods to simulate for
P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
P.mbevenz = 1;      % 1 to find even mechanical mode about z

P.kpts = 9;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 9;                           % no. of bands to solve for

P.solveasym = 1;                        % 1 to solve for antisymmetric bands
P.completeBandGaps = 1;                 % 1 to plot complete bandgaps (across all symmetries)
P.plotgeom = 1;                         % 1 to plot the geometry
P.savedat = 1;                          % 1 to save data structures
P.savebndplot = 1;                      % 1 to save bandstructure plot
P.saveplots = 1;                        % 1 to save displacement and strain profiles
P.saveMPH = 0; 
P.bandStruct_2D = 0;                 % 1 to simulate 2D band structures

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
% buildRibUnitCell(model,P);
% mphlaunch(model);
%% Single solve
datLoc = '.\test\rib\071424\';
P.datLoc = datLoc;
bds = solveBands(P);

