clear all; clc; close all;
clear P

% unit cell params
P.xsect = 'rect';%'tri';                       % beam cross sectional shape  - 'tri' or 'rect'
P.beamMat = 'diamond';                  % beam material name
P.a = 558e-9;%580e-9;                           % nominal lattice constant
P.w = 763e-9;%929e-9;                           % beam width
P.theta = 35;                           % etch angle in degrees (will not apply for P.xsect = 'rect')
P.th = 160e-9;%P.w/(2*tan(P.theta*pi/180));                          % beam thickness (define using P.theta for P.xsect = 'tri') 
P.hx = 190e-9;%250e-9;                          % nominal hole height (along x-axis)
P.hy = 500e-9;%590e-9;                          % nominal hole width (along y-axis

P.nperiod = 1;                          % no. of periods to simulate for
P.holeatedge = 1;                       % 1/0 for hole at edge/center of unit cell

P.kpts = 50;                            % no. of k-points, EXCLUDING gamma point
P.nbands = 4;                           % no. of bands to solve for
P.hole = 'elps';                        % hole shape

P.solveTM = 1;                          % 1 to solve for TM modes (if running solveBands)
P.completeBandGaps = 0;                 % 1 to plot complete bandgaps (across all symmetries)
P.plotgeom = 1;          
% 1 to plot the geometry
P.savedat = 1;                          % 1 to save data structures
P.savebndplot = 1;                      % 1 to save bandstructure plot
P.saveplots = 1;                        % 1 to save displacement and strain profiles

%% optical simulation parameters 
% solid mechanics solver parameters
P.obeveny = -1;                          % -1 for TE-like
P.obevenz = 0;                          % 1 to find even mechanical mode about z
P.obfreq = 100e12;                      % target frequency - set to 150THz for bands in telecom
P.airrad = 1550e-9;
P.anisoMat = 1;
P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                        % from <100> inplane direction about <100> surface normal


P.max_dof = 3e6;                        % max # of degrees of freedom

%% Single solve
datLoc = fullfile('.','rectOMCfinal','test_062624');
bds = solveBands(P,datLoc);
