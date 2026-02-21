clear all; clc; close all;
clear P;

%% unit cell params 
P.xsect = 'rect'; 
P.beamMat = 'diamond';                  % beam material name
P.celltype = 'hole';                   % specify the cell type
P.unitcell = 'rectrangular';                  % specify the shape of the unit cell
P.a = 650e-9;              % lattice constant 
P.hx = 343e-9;              % the diameter of the hole in x 
P.hy = 617e-9;              % the diameter of the hole in y
P.beam_width = 800e-9; % the width of the unit cell
P.d_in = 0; % the sidewall angle for the inside
P.d_out = 0; % the sidewall angle for the outside

P.th = 250e-9;             % height (along x) of cross (for celltype = 'hollow')
                            % or of inner block (for celltype = 'solid')

P.nperiod = 1;  % no. of periods to simulate for
P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
P.mbevenz = 0;      % 1 to find even mechanical mode about z

P.kpts = 5;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 10;                           % no. of bands to solve for

% simulation parameters 
P.TwoSymPlanes = 1;                     % 1 to solve the band with both y and z symmetry; 0 to solve with only z/y symmetry
P.zSymCondition = 1;                    % for the single symmetry case: 1 for z symmetry and 0 for y symmetry 
P.solveasym = 1;                        % 1 to solve for antisymmetric bands
P.completeBandGaps = 1;                 % 1 to plot complete bandgaps (across all symmetries)
P.plotgeom = 1;                         % 1 to plot the geometry
P.savedat = 1;                          % 1 to save data structures
P.savebndplot = 1;                      % 1 to save bandstructure plot
P.saveplots = 0;                        % 1 to save displacement and strain profiles
P.saveMPH = 0; 
P.bandStruct_2D = 0;                 % 1 to simulate 2D band structures
P.bandStructureDim = 1;         % 3 to simulate the band structure in 3D 

%% optical simulation parameters 
P.airDiskH = 1000e-9; % the radius of the air disk 
P.run_optical = 0; % 1 for optical simulation and 0 for mechanics 
P.optical_freq = 300; % THz the center frequency of the targeted optical bandgap 

%% mechanical simulation parameters 
% solid mechanics solver parameters
P.mbeveny = 0;                          % 1 to find even mechanical mode about y
P.mbevenz = 0;                          % 1 to find even mechanical mode about z
P.freq = 3e9;                             % target frequency - set to 0 for bandstructure simulations
P.meshSize = 5;                         % mesh quality for mechanical simulations
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
% buildHoleUnitCell(model,P);
% mphlaunch(model);
%% Single solve
currentDate = datestr(now,'mmddyyyy');
datLoc = ['.\test\hole_DiamondMechanical\',currentDate,'\'];
P.datLoc = datLoc;
solveBands(P);