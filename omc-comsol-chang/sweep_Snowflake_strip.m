clear all; clc; close all;
clear P;

%% strip unit cell parameters
b = 0;
a_lattice = 320e-9; % the lattice constant of the unit cell
b_min = -sqrt(3)*a_lattice/4; 
b_max = 0; 
b_base_list = linspace(b_min,b_max,10);
currentDate = datestr(now,'mmddyyyy');
datLoc = ['.\test\snowflake_strip_sweep\',currentDate,'\'];
P.datLoc = datLoc;

%% sweeping b
for i=1:length(b_base_list)
    sweep_snowflakeStrip_UnitCell(P,b_base_list(i));
end

function sweep_snowflakeStrip_UnitCell(P,b)
   %% unit cell params 
    P.xsect = 'rect'; 
    P.beamMat = 'diamond';                  % beam material name
    P.celltype = 'Snowflake_strip';                   % specify the cell type
    P.unitcell = 'strip';                  % specify the shape of the unit cell
    P.a = 320e-9;              % lattice constant 
    P.w = 50e-9;        % the width of the hole 
    P.b = 0e-9;              % unit cell shift in the y direction
    P.b_base = b;           % unit cell shift in the y direction (wvg region)
    P.r = 138e-9;              % radius of the unit cell
    P.th = 140e-9;             % height (along x) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r1 = 10e-9;             % width (along y) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r2 = 10e-9;              % height (along x) of each leg in cross (for celltype = 'hollow')
                                % or of outer fins (for celltype = 'solid')
    P.nperiod = 1;  % no. of periods to simulate for
    P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
    P.mbevenz = 1;      % 1 to find even mechanical mode about z

    P.bandStructureDim = 1;                 % 1D vs 2D band structure
    P.kpts = 7;                             % no. of k-points, EXCLUDING gamma point
    P.nbands = 70;                           % no. of bands to solve for

    P.solveasym = 1;                        % 1 to solve for antisymmetric bands
    P.completeBandGaps = 1;                 % 1 to plot complete bandgaps (across all symmetries)
    P.plotgeom = 0;                         % 1 to plot the geometry
    P.savedat = 1;                          % 1 to save data structures
    P.savebndplot = 1;                      % 1 to save bandstructure plot
    P.saveplots = 0;                        % 1 to save displacement and strain profiles
    P.saveMPH = 0; 
    P.bandStruct_2D = 0;                 % 1 to simulate 2D band structures

    %% mechanical simulation parameters 
    % solid mechanics solver parameters
    P.mbeveny = 1;                          % 1 to find even mechanical mode about y
    P.mbevenz = 0;                          % 1 to find even mechanical mode about z
    P.freq = 10e9;                             % target frequency - set to 0 for bandstructure simulations
    P.optical_freq = 300;                % target optical mid band frequency (THz)                        % target frequency - set to 0 for bandstructure simulations
    P.meshSize = 4;                         % mesh quality for mechanical simulations
    P.fixed_bc = 0;                       % 1 to fixed the boundaries for xz planes at y = +/- w/2

    P.anisoMat = 1;
    P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                            % from <100> inplane direction about <100> surface normal

    %% define the maximum number of degree of freedom to limit the simulation time
    P.max_dof = 3e6;                        % max # of degrees of freedom
    
    solveBands(P);
end