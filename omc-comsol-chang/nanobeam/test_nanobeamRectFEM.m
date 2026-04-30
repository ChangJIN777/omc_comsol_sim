% front panel script for nanobeam FEM simulations (example for rect cross-section, beam-hole geometry)

clear all; close all; clc

%% geometry parameters
% unit cell params
P.xsect = 'isoFit';                       % beam cross sectional shape - 'tri' or 'rect' or 'isoFit'
P.celltype = 'hole';                    % specify what type of unit cells we are simulating
P.beamMat = 'diamond';                  % beam material name
P.anisoMat = 1;

P.a = 550e-9;                           % nominal lattice constant
P.w = 600e-9;                           % beam width
P.theta = 45;                           % etch angle in degrees (no effect for rect cross section)
P.th = 400e-9;                          % beam thickness
P.hx = 220e-9;                          % nominal hole height (along x-axis)
P.hy = 330e-9;                          % nominal hole width (along y-axis)

% hole params for symmetric cavity / right half of asymmetric cavity
P.nholes = 18;                          % # holes in 1/2 beam length
P.ndef = 8;                             % # of holes in 1/2 defect region
P.maxdef = 0.133;                     % defect percentage
P.oblong = 0.752;                      % oblong parameter (zero if holes are not changed)

% cavity taper params
P.holeatctr = 0;                        % 1/0 for hole/dielectric in middle
P.taperFunc = 'cubic';                  % linear/cubic/quadratic taper function to center hole in cavity
% P.taperTo = 'custom';                 % taper to custom hole in center of cavity; disable for taper to maxdef
% P.a_ctr = 392e-9;                     % for taperTo = 'custom': lattice constant of center hole
% P.hx_ctr = 189e-9;                    % for taperTo = 'custom': hole height of center hole
% P.hy_ctr = 177e-9;                    % for taperTo = 'custom': hole width of center hole
% P.cavlen = 0e-9;                      % custom cavity length between two center holes; disable if not used

% end waveguide mirror taper params
P.wvgmir = 0;                           % no. of mirrors in end waveguide mirror taper
P.wgmTaper.func = 'cubic';              % taper function - linear, quadratic, cubic
P.wgmTaper.endtype = 'custom';          % taper to custom, maxdef, zero (or '') end hole
P.wgmTaper.a_end = 650e-9;              % for endtype = 'custom': lattice constant of end hole
P.wgmTaper.hx_end = 343e-9;             % for endtype = 'custom': hole height of end hole
P.wgmTaper.hy_end = 300e-9;             % for endtype = 'custom': hole width of end hole
% P.wgmTaper.maxdef_end = 0.094932;     % defect percentage - for endtype = 'maxdef'
% P.wgmTaper.oblong_end = 2.9172;       % oblong parameter (zero if holes are not changed) - for endtype = 'maxdef'

% asymmetric cavity - specify param data struct P.PL with similar fields to
% P, for left half of asymmetric cavity
P.asymCav = 0;                          % 1 to e-nable asymmetric cavity
if P.asymCav                
    P.PL = P;
    P.PL.nholes = 3+P.ndef;
    P.PL.wvgmir = 5;
end

P.lambda = 1500e-9;                     % target optical wavelength
P.nbeam = 2.386;                        % the refractive index of diamond at telecom

% Disorder
P.stdDev = [0,0];                       % standard deviation of hole dimensions (hh,hw)
P.stdDevPos = 0;                        % standard deviation of hole positions
P.asym = 0;                             % cross-section asymmetry (target y-offset in bottom apex position)

%% specify simulation/calculation/plot/save options
P.solveMech = 1;                        % 1 to solve for mechanics
P.solveOpt = 1;                         % 1 to solve for optics
P.calcG = 1*(P.solveMech && P.solveOpt);% 1 to calculate optomechanical coupling
P.calcS = 1*P.solveMech;                % 1 to calculate strain coupling
P.solveMechPML = 0;                     % 1 to solve for mechanical Q (future implementation)

% plotting & saving
P.plotgeom = 0;                         % 1 to plot the geometry
P.storeMPH = 0;                         % 1 to save COMSOL model file
P.plotMech = 1*P.solveMech;             % 1 to plot displacement and strain profiles
P.plotOpt = 1*P.solveOpt;               % 1 to plot E-field profiles
P.plotStrCpl = 1*P.calcS;               % 1 to plot strain coupling profile

%% mechanical simulation parameters 
% solid mechanics solver parameters
P.mevenx = 1;                           % +/-1 to find even/odd mode about x; 0 for fixed BC
P.meveny = 1;                           % +/-1 to find even/odd mode about y
P.mevenz = 1;                           % +/-1 to find even/odd mode about z
P.freq = 10e9;                           % target mechanical frequency
P.mneigs = 20;                          % # of eignevalues to find
P.mMesh = 3;                            % mesh quality for mechanical simulations
P.mAdjMesh = 1;                         % adjust mesh if DOFs exceed max_dof

% rotate crystal orientation of elasticity matrix
% ccw rotation in deg from <100> inplane direction about <100> surface normal
P.rxtal = 0;
P.rxtalInFilename = 1;

%% optical simulation parameters
% rf module solver parameters
P.oevenx = (-1)^(P.holeatctr);            % +/-1 to find even/odd optical mode about x (-1 == fundamental for hole in center)
P.oeveny = -1;                          % +/-1 to find even/odd optical mode about y (-1 == TE-like)
P.oevenz = 1;                           % +/-1 to find even/odd optical mode about z 
P.oneigs = 1;                           % # of eigenvalues to find
P.oMesh = 3;                            % mesh quality for optical simulations
P.oAdjMesh = 1;                         % adjust mesh if DOFs exceed max_dof
P.airrad = 2*P.lambda+P.w/2;            % radius of air cylinder surrounding nanobeam

%% OM coupling parameters
P.g0min = 80e3;                         % min g0 above which to save plots for

%% SiV strain coupling: susceptibilities and positions
% specify SiV axis - [1 1 1],[-1 1 1],[-1 -1 1],[1 -1 1]
% for multiple axes, specify as arrays in cell
P.zSiV = {[1 1 1]};%,[-1 1 1],[-1 -1 1],[1 -1 1]};

% specify custom coordinates to plot at
P.xSlc = 0; % relative to center
P.ySlc = 0;
P.zSlc = 0;

% define aperture to run stats (min, max, mean, stddev) on strain coupling
P.LStats.xmin = []; % leave empty to define based on generated geometry
P.LStats.xmax = [];
P.LStats.ymin = 0;
P.LStats.ymax = 60e-9;
P.LStats.zmin = P.th/2-80e-9;   % z-coords relative to center of beam
P.LStats.zmax = P.th/2;

%% Mechanical PML simulation settings (future implementation)
P.PMLmesh = 5;
P.PMLmeshDiv = 20;
P.PMLLen = 10e-6;
P.PMLstr = 0.008;

%% Simulation settings
P.max_dof = 5e6;                        % max # of degrees of freedom

%% single run
% data location to save files in
currentDate = datestr(now,'mmddyyyy');
datLoc = ['.\test\1D_OMC_hole\',currentDate,'\'];
[ds,model] = RunNanobeamFEM(P,datLoc);
