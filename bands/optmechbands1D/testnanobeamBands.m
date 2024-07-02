%% testnanobeamBands.m
% June 26, 2024
% The goal of this script is to take in a hole-y unit cell and return the
% optical and mechanical bandstructure. From here, we can extract the
% optical and mechanical bandgap features, which we can later use to
% optimize our unit cell geometry.

% This will basically be mashing together Cleaven's optbands1D and
% mechbands1D script, essentially running them in series. They are very
% similar, so it will just be a matter of figuring out how to format the
% common functions so they can be used for both optical and mechanical
% stuff.

% We will keep this separate by using 3 parameter lists: P, ofem, mfem.
% This will allow us to use the same labelling but with different values
% for the different simulation.

clear all; clc; close all;
clear P

%% Unit cell parameters
P.xsect = 'rect';%'tri';                       % beam cross sectional shape  - 'tri' or 'rect'
P.beamMat = 'diamond';                  % beam material name
P.anisoMat = 1;
P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                        % from <100> inplane direction about <100> surface normal
                                        
P.a = 640e-9;%580e-9;                           % nominal lattice constant
P.w = 870e-9;%929e-9;                           % beam width
P.theta = 35;                           % etch angle in degrees (will not apply for P.xsect = 'rect')
P.th = 160e-9;%P.w/(2*tan(P.theta*pi/180));                          % beam thickness (define using P.theta for P.xsect = 'tri') 
P.hx = 220e-9;%250e-9;                          % nominal hole height (along x-axis)
P.hy = 570e-9;%590e-9;                          % nominal hole width (along y-axis
P.hole = 'elps';                        % hole shape
P.nperiod = 1;                          % no. of periods to simulate for
P.holeatedge = 1;                       % 1/0 for hole at edge/center of unit cell
                                 
P.datLoc = 'E:\michael\omc-comsol-master-files\inital_study';

%% Optical simulation parameters
ofem.doOpt = 1;

ofem.obfreq = 100e12;                      % target frequency - set to 150THz for bands in telecom
ofem.kpts = 50;                            % no. of k-points, EXCLUDING gamma point
ofem.nbands = 4;                           % no. of bands to solve for
ofem.max_dof = 3e6;                        % max # of degrees of freedom
ofem.meshSize = 5;                         % Default meshsize

ofem.airrad = 1550e-9;
ofem.anisoMat = 1;
ofem.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                        % from <100> inplane direction about <100> surface normal
ofem.obeveny = -1;                          % -1 for TE-like
ofem.obevenz = 0;                          % 1 to find even mechanical mode about z
ofem.solveTM = 1;                          % 1 to solve for TM modes (if running solveBands)
ofem.completeBandGaps = 0;                 % 1 to plot complete bandgaps (across all symmetries)

ofem.plotgeom = 0;                         % 1 to plot the geometry
ofem.savedat = 1;                          % 1 to save data structures
ofem.savebndplot = 1;                      % 1 to save bandstructure plot
ofem.saveplots = 0;                        % 1 to save displacement and strain profiles

%% Mechanical simulation parameters
mfem.doMech = 1;

mfem.mbfreq = 0;                           % target frequency - set to 0 for bandstructure simulations
mfem.kpts = 10;                             % no. of k-points, EXCLUDING gamma point
mfem.nbands = 5;                           % no. of bands to solve for
mfem.max_dof = 3e6;                        % max # of degrees of freedom
mfem.meshSize = 5;

mfem.mbeveny = 1;                          % 1 to find even mechanical mode about y
mfem.mbevenz = 0;                          % 1 to find even mechanical mode about z
mfem.solveasym = 1;                        % 1 to solve for antisymmetric bands
mfem.completeBandGaps = 0;                 % 1 to plot complete bandgaps (across all symmetries)

mfem.plotgeom = 0;                         % 1 to plot the geometry
mfem.savedat = 1;                          % 1 to save data structures
mfem.savebndplot = 1;                      % 1 to save bandstructure plot
mfem.saveplots = 0;                        % 1 to save displacement and strain profiles

if ofem.doOpt
    tic
    P.optmech = 1;                      % 1/0 for optical/mechanical simulation
    disp('Starting optical simulation')
    opt_bds = solveOpticalBands(P,ofem);
    toc
end

if mfem.doMech
    tic
    P.optmech = 0;
    disp('Starting mechanical simulation')
    mech_bds = solveMechBands(P,mfem);
    toc
end

disp('Optical Bandgap properties')
disp(['Midgap frequency: ', num2str(opt_bds.midGap/1e12), ' THz']);
disp(['Gap size: ', num2str(100*opt_bds.gapSize/opt_bds.midGap), ' %']);

disp('Mechanical Bandgap properties')
disp(['Midgap frequency: ', num2str(mech_bds.midGap/1e9), ' GHz']);
disp(['Gap size: ', num2str(100*mech_bds.gapSize/mech_bds.midGap), ' %']);