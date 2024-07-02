clear all; clc; close all;
clear P

% unit cell params
P.xsect = 'rect';                    % beam cross sectional shape  - 'tri' or 'rect'
P.beamMat = 'diamond';              % beam material name
P.celltype = 'blockTet';                    % hole shape
P.a = 1000e-9;                           % nominal lattice constant
P.w = 1000e-9;                          % beam width
P.th = 230e-9;                       % beam thickness
P.hx = 670e-9;                          % nominal hole height (along x-axis)
P.hy = 150e-9;


P.filRad = 0e-9;

P.nperiod = 1;                      % no. of periods to simulate for
P.holeatedge = 0;                   % 1/0 for hole at edge/center of unit cell

P.kpts = 10;                         % no. of k-points, EXCLUDING gamma point
P.nbands = 5;                       % no. of bands to solve for

P.solveasym = 1;                    % 1 to solve for antisymmetric bands
P.completeBandGaps = 1;             % 1 to plot complete bandgaps (across all symmetries)
P.plotgeom = 1;                     % 1 to plot the geometry
P.savedat = 1;                      % 1 to save data structures
P.savebndplot = 1;                  % 1 to save bandstructure plot
P.saveplots = 1;                    % 1 to save displacement and strain profiles

%% mechanical simulation parameters 
% solid mechanics solver parameters
P.mbeveny = 1;                      % 1 to find even mechanical mode about y
P.mbevenz = 1;                      % 1 to find even mechanical mode about z
P.mbfreq = 0;                       % target frequency

% % material properties
% P.E = 1050e9;                       % Young's modulus
% P.rho = 3500;                       % density
% P.nu = 0.2;                         % Poisson ratio
%  
% % Anisotropic elasticity matrix - COMSOL v4+ ordering: 
% % 11, 12,22, 13,23,33, 14,24,34,44, 15,25,35,45,55, 16,26,36,46,56,66
% % disable P.D is using isotropic material (defined by E, nu, rho)
% c11 = 1076e9;
% c12 = 125e9;
% c44 = 578e9;
% P.D = [c11, c12,c11, c12,c12,c11, 0,0,0,c44, 0,0,0,0,c44, 0,0,0,0,0,c44];
P.anisoMat = 1;
P.rxtal = 45;                       % ccw rotation of elasticity matrix in deg 
                                    % from <100> inplane direction about <100> surface normal

%%
P.max_dof = 3e6;                      % max # of degrees of freedom

%% Reference
% datLoc = 'D:\Files\PhonXtal\BlockTet\RefCell\';
% bds = solveBands(P,datLoc);
datLocMain = 'D:\Files\PhonXtal\BlockTet\';

%% Sweep thickness
% datLoc = [datLocMain,'SweepTh\'];
% thAll = [50,100,200,300,400,500]*1e-9;
% for th = thAll
%     Pd = P;
%     Pd.th = th;
%     bds = solveBands(Pd,datLoc);
%     
%     % compile gamma points
%     
% end

% plot
% for th = thAll
%     load
% end

%% Sweep hx/a
% datLoc = [datLocMain,'SweepHxa\'];
% hxaAll = 0.1:0.1:0.9;
% for hxa = hxaAll
%     Pd = P;
%     Pd.hx = hxa*Pd.a;
%     bds = solveBands(Pd,datLoc);
% end

%% Sweep hy/w
% datLoc = 'D:\Files\PhonXtal\BlockTet\SweepHy_w\';
% hywAll = 0.1:0.1:0.9;
% for hyw = hywAll
%     Pd = P;
%     Pd.hy = hyw*Pd.a;
%     bds = solveBands(Pd,datLoc);
% end

%% Sweep a, keep hx/a
datLoc = [datLocMain,'Sweepa\'];
hxa0 = P.hx/P.a;
aAll = (750:50:1250)*1e-9;
for a = aAll
    Pd = P;
    Pd.a = a;
    Pd.hx = hxa0*Pd.a;
    bds = solveBands(Pd,datLoc);
end

%% Sweep w, keep hy/w
% datLoc = [datLocMain,'Sweepw\'];
% hyw0 = P.hy/P.w;
% wAll = (750:50:1250)*1e-9;
% for w = wAll
%     Pd = P;
%     Pd.w = w;
%     Pd.hy = hyw0*Pd.w;
%     bds = solveBands(Pd,datLoc);
% end

