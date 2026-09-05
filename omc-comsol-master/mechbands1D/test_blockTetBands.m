clear all; clc; close all;
clear P

% unit cell params
P.xsect = 'rect';                       % beam cross sectional shape  - 'tri' or 'rect'
P.beamMat = 'diamond';                  % beam material name
P.celltype = 'blockTet';                % hole shape - comment out this line for 
P.a = 1000e-9;                      % nominal lattice constant
P.w = 1500e-9;                          % width of block
P.th = 230e-9;                          % beam thickness
P.hx = 670e-9;                      % length of block (along x-axis)
P.hy = 150e-9;                          % width of tether (along y-axis)

P.filRad = 0e-9;                        % fillet radius to simulate rounding of corners

P.nperiod = 1;                          % no. of periods to simulate for
P.holeatedge = 0;                       % 1/0 for hole at edge/center of unit cell

P.kpts = 10;                            % no. of k-points, including gamma point
P.nbands = 5;                           % no. of bands to solve for


P.solveasym = 1;                        % 1 to solve for antisymmetric bands
P.completeBandGaps = 1;                 % 1 to plot complete bandgaps (across all symmetries)
P.plotgeom = 1;                         % 1 to plot the geometry
P.savedat = 1;                          % 1 to save data structures
P.savebndplot = 1;                      % 1 to save bandstructure plot
P.saveplots = 1;                        % 1 to save displacement and strain profiles

%% mechanical simulation parameters 
% solid mechanics solver parameters
P.mbeveny = 1;                          % 1 to find even mechanical mode about y
P.mbevenz = 1;                          % 1 to find even mechanical mode about z
P.mbfreq = 0;                           % target frequency - set to 0 for bandstructure simulations

P.anisoMat = 1;
P.rxtal = 0;                            % ccw rotation of elasticity matrix in deg 
                                        % from <100> inplane direction about <100> surface normal

%%
P.max_dof = 3e6;                        % max # of degrees of freedom

%% Single solve
datLoc = [fullfile('.','test','BlockTet','ASN'),filesep];
bds = solveBands(P,datLoc);

%% Loop
datLoc = 'D:\Files\PhonXtal\BlockTet\ASN\';
for w = (1.1:0.1:2)
    P.w = w*1e-6;
    bds = solveBands(P,datLoc);
end
