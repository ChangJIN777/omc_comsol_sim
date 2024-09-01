%% sweeping code 
clear all; clc; close all;
clear P;
P = struct;

w_ratio_min = 0.2;
w_ratio_max = 0.8;
h_ratio_min = 0.2;
h_ratio_max = 0.8;
w_ratio_list = linspace(w_ratio_min,w_ratio_max,10);
h_ratio_list = linspace(h_ratio_min,h_ratio_max,5);
h_ratio_list = [0.8];

%% sweep the lattice
for i=1:length(w_ratio_list)
    for j=1:length(h_ratio_list)
        sweep_boomerang_lower(P,w_ratio_list(i),h_ratio_list(j));
    end
end

function sweep_boomerang_lower(P,w_ratio,h_ratio) 
    P.celltype = 'boomerang_strip_v2';                   % specify the cell type
    P.xsect = 'rect'; 
    P.beamMat = 'diamond';                  % beam material name
    P.unitcell = 'hexagonal';                  % specify the shape of the unit cell
    %% unit cell params 
    P.a = 400e-9;              % lattice constant 
    P.w = 86e-9;              % unit cell width (along x)
    P.r = 160e-9;              % unit cell height (along y)
    P.b = sqrt(3)*P.a/4;        % removing one row of the 2D unit cells
    P.th = 180e-9;             % height (along x) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.wo = P.a*0.8;           % the height of the hole in the lower portion
    P.wi = P.wo*w_ratio; 
    P.ho = P.a*sqrt(3)/2;
    P.hi = h_ratio*P.ho;
    P.r1 = 10e-9;             % width (along y) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r2 = 10e-9;              % height (along x) of each leg in cross (for celltype = 'hollow')
                                % or of outer fins (for celltype = 'solid')
    P.nperiod = 1;  % no. of periods to simulate for
    P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
    P.mbevenz = 1;      % 1 to find even mechanical mode about z
    
    P.kpts = 9;                             % no. of k-points, EXCLUDING gamma point
    P.nbands = 20;                           % no. of bands to solve for
    
    P.solveasym = 1;                        % 1 to solve for antisymmetric bands
    P.completeBandGaps = 1;                 % 1 to plot complete bandgaps (across all symmetries)
    P.plotgeom = 1;                         % 1 to plot the geometry
    P.savedat = 1;                          % 1 to save data structures
    P.savebndplot = 1;                      % 1 to save bandstructure plot
    P.saveplots = 1;                        % 1 to save displacement and strain profiles
    P.saveMPH = 0; 
    P.bandStruct_2D = 0;                 % 1 to simulate 2D band structures
    
    %% mechanical simulation parameters 
    % solid mechanics solver parameters
    P.mbeveny = 0;                          % 1 to find even mechanical mode about y
    P.mbevenz = 1;                          % 1 to find even mechanical mode about z
    P.freq = 0;                             % target frequency - set to 0 for bandstructure simulations
    P.meshSize = 4;                         % mesh quality for mechanical simulations
    P.fixed_bc = 1;                       % 1 to fixed the boundaries for xz planes at y = +/- w/2
    
    P.anisoMat = 1;
    P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                            % from <100> inplane direction about <100> surface normal
    
    %% define the maximum number of degree of freedom to limit the simulation time
    P.max_dof = 3e6;                        % max # of degrees of freedom
    
    %% run the simluation and save the data
    datLoc = '.\test\boomerang_strip_v2\082124_sweep2\';
    P.datLoc = datLoc;
    bds = solveBands(P);
end 

function sweep_boomerang_lower_v2(P,a)    
    %% unit cell params 
    P.a = a;              % lattice constant 
    P.h = 200*a/450;           % the height of the hole in the lower portion
    P.d = 86e-9;           % the width of the hole in the lower portion
    P.w = 86e-9;              % unit cell width (along x)
    P.r = a*180/450;              % unit cell height (along y)
    P.th = 180e-9;             % height (along x) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r1 = 10e-9;             % width (along y) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r2 = 10e-9;              % height (along x) of each leg in cross (for celltype = 'hollow')
                                % or of outer fins (for celltype = 'solid')
    P.nperiod = 1;  % no. of periods to simulate for
    P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
    P.mbevenz = 1;      % 1 to find even mechanical mode about z
    
    P.kpts = 7;                             % no. of k-points, EXCLUDING gamma point
    P.nbands = 15;                           % no. of bands to solve for
    
    P.solveasym = 1;                        % 1 to solve for antisymmetric bands
    P.completeBandGaps = 1;                 % 1 to plot complete bandgaps (across all symmetries)
    P.plotgeom = 0;                         % 1 to plot the geometry
    P.savedat = 1;                          % 1 to save data structures
    P.savebndplot = 1;                      % 1 to save bandstructure plot
    P.saveplots = 1;                        % 1 to save displacement and strain profiles
    P.saveMPH = 0; 
    P.bandStruct_2D = 0;                 % 1 to simulate 2D band structures
    
    %% mechanical simulation parameters 
    % solid mechanics solver parameters
    P.mbeveny = 0;                          % 1 to find even mechanical mode about y
    P.mbevenz = 1;                          % 1 to find even mechanical mode about z
    P.freq = 0;                             % target frequency - set to 0 for bandstructure simulations
    P.meshSize = 4;                         % mesh quality for mechanical simulations
    P.fixed_bc = 1;                       % 1 to fixed the boundaries for xz planes at y = +/- w/2
    
    P.anisoMat = 1;
    P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                            % from <100> inplane direction about <100> surface normal
    
    %% define the maximum number of degree of freedom to limit the simulation time
    P.max_dof = 3e6;                        % max # of degrees of freedom
    
    %% run the simluation and save the data
    datLoc = '.\test\boomerang_lower\082024_sweep1\';
    P.datLoc = datLoc;
    bds = solveBands(P);
end 