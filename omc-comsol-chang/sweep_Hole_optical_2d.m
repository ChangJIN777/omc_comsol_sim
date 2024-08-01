clear all; clc; close all;
clear P;

%% unit cell parameters
r = 115e-9;
r_min = 80e-9; 
r_max = 120e-9;
r_list = linspace(r_min,r_max,5);
a = 342e-9;
a_min = a*0.8;

b = 0;
b_min = 10e-9; 
b_max = 20e-9; 
b_list = linspace(b_min,b_max,5);
datLoc = '.\test\hole\073124_optical_sweep_2\';
P.datLoc = datLoc;

%% sweeping a 
for i=1:length(a_list)
    sweep_Hole_2dUnitCell(P,a_list(i),r,b);
end

%% sweeping r 
for i=1:length(r_list)
    sweep_Hole_2dUnitCell(P,a,r_list(i),b);
end
%% sweeping b
for i=1:length(b_list)
    sweep_Hole_2dUnitCell(P,a,r,b_list(i));
end

function sweep_Hole_2dUnitCell(P,a,r,b)
    %% unit cell params 
    P.xsect = 'rect'; 
    P.beamMat = 'silicon';                  % beam material name
    P.celltype = 'hole';                   % specify the cell type
    P.unitcell = 'hexagonal';                  % specify the shape of the unit cell
    P.a = a;              % lattice constant 
    P.b = b;              % unit cell shift in the y direction
    P.r = r;              % radius of the unit cell
    P.th = 220e-9;             % height (along x) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r1 = 35e-9;             % width (along y) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r2 = 40e-9;              % height (along x) of each leg in cross (for celltype = 'hollow')
                                % or of outer fins (for celltype = 'solid')
    P.nperiod = 1;  % no. of periods to simulate for
    P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
    P.mbevenz = 1;      % 1 to find even mechanical mode about z

    P.bandStructureDim = 2;                 % 1D vs 2D band structure
    P.kpts = 9;                             % no. of k-points, EXCLUDING gamma point
    P.nbands = 9;                           % no. of bands to solve for

    P.solveasym = 1;                        % 1 to solve for antisymmetric bands
    P.completeBandGaps = 1;                 % 1 to plot complete bandgaps (across all symmetries)
    P.plotgeom = 0;                         % 1 to plot the geometry
    P.savedat = 1;                          % 1 to save data structures
    P.savebndplot = 1;                      % 1 to save bandstructure plot
    P.saveplots = 0;                        % 1 to save displacement and strain profiles
    P.saveMPH = 0; 
    P.bandStruct_2D = 1;                 % 1 to simulate 2D band structures

    %% mechanical simulation parameters 
    % solid mechanics solver parameters
    P.mbeveny = 0;                          % 1 to find even mechanical mode about y
    P.mbevenz = 1;                          % 1 to find even mechanical mode about z
    P.optical_freq = 220;                % target optical mid band frequency (THz)                        % target frequency - set to 0 for bandstructure simulations
    P.meshSize = 3;                         % mesh quality for mechanical simulations
    P.fixed_bc = 0;                       % 1 to fixed the boundaries for xz planes at y = +/- w/2

    P.anisoMat = 1;
    P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                            % from <100> inplane direction about <100> surface normal

    %% define the maximum number of degree of freedom to limit the simulation time
    P.max_dof = 3e6;                        % max # of degrees of freedom
    
    solveOpticalBands(P);
end