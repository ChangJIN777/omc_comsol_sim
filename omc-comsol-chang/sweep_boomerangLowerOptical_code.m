%% sweeping code 
P.celltype = 'boomerang_lower';                   % specify the cell type

a = 480e-9; % lattice constant
h_min = 110e-9;
h_max = 140e-9;
d_min = 86e-9;
d_max = 150e-9;
h_list = linspace(h_min,h_max,5);
d_list = linspace(d_min,d_max,5);

%% run the sweep
for i=1:length(h_list)
    for j=1:length(d_list)
        sweep_boomerang_lower_optical(P,h_list(i),d_list(j));
    end
end

function sweep_boomerang_lower_optical(P,h,d)    
    %% unit cell params 
    P.xsect = 'rect'; 
    P.beamMat = 'diamond';                  % beam material name
    P.celltype = 'boomerang_lower';                   % specify the cell type
    P.unitcell = 'square';                  % specify the shape of the unit cell
    P.a = 400e-9;              % lattice constant 
    P.h = h;           % the height of the hole in the lower portion
    P.d = d;           % the width of the hole in the lower portion
    P.w = 86e-9;              % unit cell width (along x)
    P.r = 160e-9;              % unit cell height (along y)
    P.th = 180e-9;             % height (along x) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r1 = 35e-9;             % width (along y) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r2 = 40e-9;              % height (along x) of each leg in cross (for celltype = 'hollow')
                                % or of outer fins (for celltype = 'solid')
    P.nperiod = 1;  % no. of periods to simulate for
    P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
    P.mbevenz = 1;      % 1 to find even mechanical mode about z
    
    P.kpts = 10;                             % no. of k-points, EXCLUDING gamma point
    P.nbands = 9;                           % no. of bands to solve for
    P.bandStructureDim = 1;                 % 1D vs 2D band structure

    
    P.solveasym = 1;                        % 1 to solve for antisymmetric bands
    P.completeBandGaps = 1;                 % 1 to plot complete bandgaps (across all symmetries)
    P.plotgeom = 1;                         % 1 to plot the geometry
    P.savedat = 1;                          % 1 to save data structures
    P.savebndplot = 1;                      % 1 to save bandstructure plot
    P.saveplots = 0;                        % 1 to save displacement and strain profiles
    P.saveMPH = 0; 
    P.bandStruct_2D = 1;                 % 1 to simulate 2D band structures
    
    %% optical simulation parameters 
    % solid mechanics solver parameters
    P.mbeveny = 0;                          % 1 to find even mechanical mode about y
    P.mbevenz = 1;                          % 1 to find even mechanical mode about z
    P.freq = 0;                             % target frequency - set to 0 for bandstructure simulations
    P.meshSize = 4;                         % mesh quality for mechanical simulations
    P.fixed_bc = 0;                       % 1 to fixed the boundaries for xz planes at y = +/- w/2
    P.optical_freq = 220;                % target optical mid band frequency (THz)                        % target frequency - set to 0 for bandstructure simulations

    
    P.anisoMat = 1;
    P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                            % from <100> inplane direction about <100> surface normal
    
    %% define the maximum number of degree of freedom to limit the simulation time
    P.max_dof = 3e6;                        % max # of degrees of freedom
    
    %% run the simluation and save the data
    datLoc = '.\test\boomerang_lower\072824_optical_sweep\';
    P.datLoc = datLoc;
    bds = solveOpticalBands(P);
end 