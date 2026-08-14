clear all; clc; close all;
clear P

% unit cell params
P.xsect = 'rect';                       % beam cross sectional shape  - 'tri' or 'rect'
P.beamMat = 'diamond';                  % beam material name
P.a = 528e-9;                           % nominal lattice constant
P.w = 936e-9;                           % beam width
P.theta = 35;                           % etch angle in degrees (will not apply for P.xsect = 'rect')
P.th = 160e-9;                          % beam thickness (define using P.theta for P.xsect = 'tri') 
P.hx = 197e-9;                          % nominal hole height (along x-axis)
P.hy = 578e-9;                          % nominal hole width (along y-axis

P.nperiod = 1;                          % no. of periods to simulate for
P.holeatedge = 1;                       % 1/0 for hole at edge/center of unit cell

P.kpts = 9;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 9;                           % no. of bands to solve for
P.hole = 'elps';                        % hole shape

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
P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                        % from <100> inplane direction about <100> surface normal

%% define the maximum number of degree of freedom to limit the simulation time
P.max_dof = 3e6;                        % max # of degrees of freedom

%% Single solve
datLoc = [fullfile('.','test','Holely','070424'),filesep];
bds = solveBands(P,datLoc);

% %% change thickness
% P.th = 200e-9;
% datLoc = 'D:\Files\OMC-SiV\RectOMC\Bands\th200_';
% 
% %% Sweep hx
% PS = P;
% datLocS = [datLoc,'swphx\'];
% for hxr = 0.1:0.1:0.5
%     PS.hx = PS.a*hxr; 
%     bds = solveBands(PS,datLocS);
% end
% 
% %% Sweep hy
% PS = P;
% datLocS = [datLoc,'swphy\'];
% for hyr = 0.3:0.1:0.8
%     PS.hy = PS.w*hyr; 
%     bds = solveBands(PS,datLocS);
% end
% 
% %% Sweep a, fixed hx/a
% PS = P;
% datLocS = [datLoc,'swpa\'];
% hxr = P.hx/P.a;
% for a = (400:50:650)*1e-9
%     PS.a = a; 
%     PS.hx = PS.a*hxr;
%     bds = solveBands(PS,datLocS);
% end
% 
% %% Sweep w, fixed hy/a
% PS = P;
% datLocS = [datLoc,'swpw\'];
% hyr = P.hy/P.w;
% for w = (550:50:800)*1e-9
%     PS.w = w; 
%     PS.hy = PS.w*hyr;
%     bds = solveBands(PS,datLocS);
% end
% 
% %% Sweep th
% PS = P;
% datLocS = [datLoc,'swpth\'];
% for th = (200:50:400)*1e-9
%     PS.th = th; 
%     bds = solveBands(PS,datLocS);
% end
% 
