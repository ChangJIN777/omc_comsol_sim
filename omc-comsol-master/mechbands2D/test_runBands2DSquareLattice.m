% script to run mechanical bandstructure simulations for 2D unit cell in square lattice

clear all; close all;
clear P

P.beamMat = 'diamond';      % beam material name

% type of unit cell: 
% 'solidX' for solid cross (inner block + small fins) 
% 'hollowX' for block with cross-shaped hole (block - cross)
% 'circTets' for circle with tethers
P.celltype = 'sqTets';     

% unit cell geometry params
P.a = 135e-9;              % lattice constant (along x)
P.w = 135e-9;              % unit cell width (along y)
P.th = 75e-9;              % unit cell thickness (along z)

% For solidX, circTets:   
% cx/cy = length along x/y of central square/elliptical (for solidX/circTets) block
% tx/ty = width along x/y at end oftethers
% For hollowX:
% cx/cy = length along x/y of hollow cross
% tx/ty = width along x/y of each leg of hollow cross
P.cx = 115e-9;%*sqrt(2*pi)/2;
P.cy = 115e-9;%*sqrt(2*pi)/2;
P.tx = 10e-9;
P.ty = 10e-9;

P.filRad = 30e-9;           % fillet radius to round inner vertices of unit cell

% symmetry
P.mbevenz = 1;              % +/-1 for symmetry/antisymmetry in z; 0 for no symmetry

% material properties
P.E = 1050e9;                % Young's modulus
P.rho = 3500;               % density in kg/m^3
P.nu = 0.2;                % Poisson's ratio
 
% Anisotropic elasticity matrix - COMSOL v4+ ordering: 
% [11, 12,22, 13,23,33, 14,24,34,44, 15,25,35,45,55, 16,26,36,46,56,66]
c11 = 1076e9;
c12 = 125e9;
c44 = 578e9;
P.D = [c11, c12,c11, c12,c12,c11, 0,0,0,c44, 0,0,0,0,c44, 0,0,0,0,0,c44];   % comment this line if using isotropic material
P.rxtal = 45;           % ccw rotation in deg from <100> inplane direction about <100> surface normal


% bandstructure simulation params
P.nbands = 8;               % no. of bands
P.mbfreq = 0;               % target frequency
P.kpts = 5;                 % no. of k-points to simulate for along each sweep (INCLUDING Gamma/X/M points)
                            % for full GMX sweep, total no. of k-points solved for = 3*(kpts-1)
P.max_dof = 3e6;            % max # of degrees of freedom
P.mbMesh = 5;               % mesh quality
P.solveasym = 1;            % 1 to solve for symmetric and antisymmetric bands

% plotting
P.plotgeom = 1;             % 1 to plot unit cell geometry
P.savegeom = 1;             % 1 to save unit cell geometry
P.saveplots = 1;            % 1 to save displacement plots
P.saveMPH = 1;              % 1 to save COMSOL model file
P.savebndplot = 1;          % 1 to save bandstructure plot
P.savedat = 1;              % 1 to save data structure containing results

%% Run
datLoc = [fullfile('.','test2DsqTets'),filesep];
bds = solveBands2DSquareLattice(P,datLoc);
