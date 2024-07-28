% front panel script for nanobeam FEM simulations (example for rectangular cross-section, block-tether geometry)
% current implementation (Jun 26 2019): for one defect cell in center
% between mirrors in symmetric geometry

clear all; close all; clc

% script folder containing scripts common to FDTD and FEM simulations
% P.scriptLoc = 'L:\Individuals\cchia\OMC\Utility Scripts\';  % comment out if unused

%% geometry parameters
% unit cell params
P.xsect = 'rect';                           % beam cross sectional shape - 'tri' or 'rect'
P.beamMat = 'diamond';                      % beam material name
P.celltype = 'blockTet';                    % hole shape
P.anisoMat = 1;

P.a = 1000e-9;                              % nominal lattice constant
P.w = 1000e-9;                              % beam width
% P.theta = 45;                             % etch angle in degrees (no effect for rect cross section)
P.th = 230e-9;                              % beam thickness
P.hx = 670e-9;                              % nominal block length (along x-axis)
P.hy = 150e-9;                              % nominal tether width (along y-axis)

% hole params for symmetric cavity / right half of asymmetric cavity
P.nholes = 9;                               % # holes in 1/2 beam length
P.ndef = 1;                                 % # of holes in 1/2 defect region
P.maxdef = 0;                               % defect percentage
P.oblong = 0;                               % oblong parameter (zero if holes are not changed)

% cavity taper params
% scl = 3.4/5;
P.holeatctr = 0;                            % 1/0 for hole/dielectric in middle of cavity
% P.taperFunc = 'cubic';                    % linear/cubic/quadratic taper function to center cell in cavity
P.taperTo = 'custom';                       % taper to custom cell in center of cavity; disable this to taper to maxdef
P.a_ctr = P.a;                          % for taperTo = 'custom': lattice constant of center cell
P.hx_ctr = P.hx;                        % for taperTo = 'custom': length of block of center cell
P.hy_ctr = P.hy;                            % for taperTo = 'custom': width of tethers from center cell
P.w_ctr = 1700e-9;                    % for taperTo = 'custom': width of center cell
% P.cavlen = 0e-9;                          % custom cavity length between two center holes; disable if not used

% functionality for left-right asymmetric cavity not implemented for
% block-tether geometry (as of Jun 26, 2019)
P.asym = 0;
P.asymCav = 0;

%% specify simulation/calculation/plot/save options
P.solveMech = 1;                            % 1 to solve for mechanics
P.solveOpt = 0;                             % 1 to solve for optics
P.calcG = 1*(P.solveMech && P.solveOpt);    % 1 to calculate optomechanical coupling
P.calcS = 1*P.solveMech;                    % 1 to calculate strain coupling
P.solveMechPML = 0;                         % 1 to solve for mechanical Q (future implementation)

% plotting & saving
P.plotgeom = 1;                             % 1 to plot the geometry
P.storeMPH = 1;                             % 1 to save COMSOL model file
P.plotMech = 1*P.solveMech;                 % 1 to plot displacement and strain profiles
P.plotOpt = 1*P.solveOpt;                   % 1 to plot E-field profiles
P.plotStrCpl = 1*P.calcS;                   % 1 to plot strain coupling profile

%% mechanical simulation parameters 
% solid mechanics solver parameters
P.mevenx = 1;                               % +/-1 to find even/odd mode about x; 0 for fixed BC
P.meveny = 1;                               % +/-1 to find even/odd mode about y
P.mevenz = 1;                              % +/-1 to find even/odd mode about z
P.freq = 5e9;                               % target mechanical frequency
P.mneigs = 10;                              % # of eignevalues to find
P.mMesh = 1;                                % mesh quality for mechanical simulations
P.mAdjMesh = 1;                             % adjust mesh if DOFs exceed max_dof

% rotate crystal orientation of elasticity matrix
% ccw rotation in deg from <100> inplane direction about <100> surface normal
P.rxtal = 45;                                
P.rxtalInFilename = 1;                      % option to include crystal orientation in filename

%% SiV strain coupling: susceptibilities and positions
% specify SiV axis - [1 1 1],[-1 1 1],[-1 -1 1],[1 -1 1]
% for multiple axes, specify as arrays in cell
P.zSiV = {[1 1 1],[-1 1 1],[-1 -1 1],[1 -1 1]};

% specify custom coordinates to plot at
P.xSlc = 0;
P.ySlc = 0;
P.zSlc = 0;

% define aperture to run stats (min, max, mean, stddev) on strain coupling
LStats.xmin = 0;
LStats.xmax = P.hx_ctr/2;
LStats.ymin = 0;
LStats.ymax = P.hy_ctr/2+50e-9;
LStats.zmin = P.th/2-100e-9;
LStats.zmax = P.th/2-50e-9;
P.LStats = LStats;
%% Simulation settings
P.max_dof = 5e6;                        % max # of degrees of freedom

%% single run
% data location to save files in
datLoc = 'D:\Files\SpinPhon\BTComp\'; 
[ds,model] = RunNanobeamFEM(P,datLoc);
