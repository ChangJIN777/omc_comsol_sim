clear all; clc; close all;
clear P

% unit cell params
P.xsect = 'tri';                       % beam cross sectional shape  - 'tri' or 'rect'
P.beamMat = 'diamond';                  % beam material name
P.a = 565.7e-9;                           % nominal lattice constant
P.w = 1074.8e-9;                           % beam width
P.theta = 45;                           % etch angle in degrees (will not apply for P.xsect = 'rect')
P.th = P.w/(2*tan(P.theta*pi/180));                          % beam thickness (define using P.theta for P.xsect = 'tri') 
P.hx = 226.3e-9;                          % nominal hole height (along x-axis)
P.hy = 698.6e-9;                          % nominal hole width (along y-axis

P.nperiod = 1;                          % no. of periods to simulate for
P.holeatedge = 1;                       % 1/0 for hole at edge/center of unit cell

P.kpts = 9;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 9;                           % no. of bands to solve for
P.hole = 'elps';                        % hole shape

P.solveasym = 1;                        % 1 to solve for antisymmetric bands
P.completeBandGaps = 0;                 % 1 to plot complete bandgaps (across all symmetries)
P.plotgeom = 1;                         % 1 to plot the geometry
P.savedat = 1;                          % 1 to save data structures
P.savebndplot = 1;                      % 1 to save bandstructure plot
P.saveplots = 1;                        % 1 to save displacement and strain profiles

%% mechanical simulation parameters 
% solid mechanics solver parameters
P.mbeveny = 1;                          % 1 to find even mechanical mode about y
P.mbevenz = 0;                          % 1 to find even mechanical mode about z
P.mbfreq = 0;                           % target frequency - set to 0 for bandstructure simulations

P.anisoMat = 1;
P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                        % from <100> inplane direction about <100> surface normal

%%
P.max_dof = 3e6;                        % max # of degrees of freedom

%% Single solve
datLoc = 'D:\Files\OMC-SiV\45oBands\';
bds = solveBands(P,datLoc);
