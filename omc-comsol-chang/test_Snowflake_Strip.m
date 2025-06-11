clear all; clc; close all;
clear P;

%% unit cell params 
P.xsect = 'rect'; 
P.beamMat = 'diamond';                  % beam material name
P.celltype = 'Snowflake_strip';                   % specify the cell type
P.unitcell = 'strip';                  % specify the shape of the unit cell
P.a = 660e-9;              % lattice constant 
P.w = 82.5e-9;        % the width of the hole 
P.b = sqrt(3)*P.a/2;           % unit cell shift in the y direction (wvg region)
P.r = 272e-9;              % radius of the uni t cell
P.th = 350e-9;             % height (along x) of cross (for celltype = 'hollow')
                            % or of inner block (for celltype = 'solid')
P.r1 = 10e-9;             % width (along y) of cross (for celltype = 'hollow')
                            % or of inner block (for celltype = 'solid')
P.r2 = 10e-9;              % height (along x) of each leg in cross (for celltype = 'hollow')
                            % or of outer fins (for celltype = 'solid')

P.wo = 378e-9;           % the height of the hole in the lower portion
P.wi = 192e-9;           % the width of the hole in the lower portion                            
P.ho = 260e-9;
P.hi = 190e-9;
P.b = sqrt(3)*P.a/2;        
P.d = 100e-9;

P.nperiod = 1;  % no. of periods to simulate for
P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell

P.bandStructureDim = 1;                 % 1D vs 2D band structure
P.kpts = 7;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 40;                           % no. of bands to solve for

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
P.TwoSymPlanes = 1;
P.mbeveny = 1;                          % 1 to find even mechanical mode about y
P.mbevenz = 1;                          % 1 to find even mechanical mode about z
P.zSymCondition = 0;                     % 1 to setup symmetry in the z direction and 0 in the y direction
P.freq = 10e9;                             % target frequency - set to 0 for bandstructure simulations
P.optical_freq = 180;                % target optical mid band frequency (THz)                        % target frequency - set to 0 for bandstructure simulations
P.meshSize = 8;                         % mesh quality for mechanical simulations
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
% buildSnowflakeStrip_3D(model,P);
% mphlaunch(model);

%% run the simulation
% datLoc = '.\test\snowflake_optical\122124_optical\';
currentDate = datestr(now,'mmddyyyy');
datLoc = ['.\test\snowflake_strip\',currentDate,'\'];
P.datLoc = datLoc;
solveBands(P);