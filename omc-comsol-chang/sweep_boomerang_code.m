%% sweeping code 
P.celltype = 'boomerang';                   % specify the cell type

a = 480e-9; % lattice constant
w_min = 0.2*a; 
w_max = 0.3*a;
r_min = 0.2*a;
r_max = 0.45*a;
w_list = linspace(w_min,w_max,3);
r_list = linspace(w_min,w_max,3);

%% run the sweep
for i=1:length(w_list)
    for j=1:length(r_list)

        sweep_boomerang(P,w_list(i),r_list(j))
    end
end

function sweep_boomerang(P,w,r)    
    %% unit cell params 
    P.xsect = 'rect'; 
    P.beamMat = 'diamond';                  % beam material name
    P.celltype = 'boomerang';                   % specify the cell type
    P.unitcell = 'hexagonal';                  % specify the shape of the unit cell
    P.a = 480e-9;              % lattice constant 
    % P.w = 93e-9;              % unit cell width (along x)
    % P.r = 172e-9;              % unit cell height (along y)
    P.w = w;              % unit cell width (along x)
    P.r = r;              % unit cell height (along y)
    P.th = 180e-9;             % height (along x) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r1 = 35e-9;             % width (along y) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r2 = 40e-9;              % height (along x) of each leg in cross (for celltype = 'hollow')
                                % or of outer fins (for celltype = 'solid')
    P.nperiod = 1;  % no. of periods to simulate for
    P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
    P.mbevenz = 1;      % 1 to find even mechanical mode about z
    
    P.kpts = 9;                             % no. of k-points, EXCLUDING gamma point
    P.nbands = 8;                           % no. of bands to solve for
    
    P.solveasym = 1;                        % 1 to solve for antisymmetric bands
    P.completeBandGaps = 1;                 % 1 to plot complete bandgaps (across all symmetries)
    P.plotgeom = 1;                         % 1 to plot the geometry
    P.savedat = 1;                          % 1 to save data structures
    P.savebndplot = 1;                      % 1 to save bandstructure plot
    P.saveplots = 0;                        % 1 to save displacement and strain profiles
    P.saveMPH = 0; 
    P.bandStruct_2D = 1;                 % 1 to simulate 2D band structures
    
    %% mechanical simulation parameters 
    % solid mechanics solver parameters
    P.mbeveny = 0;                          % 1 to find even mechanical mode about y
    P.mbevenz = 1;                          % 1 to find even mechanical mode about z
    P.freq = 0;                             % target frequency - set to 0 for bandstructure simulations
    P.meshSize = 4;                         % mesh quality for mechanical simulations
    P.fixed_bc = 0;                       % 1 to fixed the boundaries for xz planes at y = +/- w/2
    
    P.anisoMat = 1;
    P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                            % from <100> inplane direction about <100> surface normal
    
    %% define the maximum number of degree of freedom to limit the simulation time
    P.max_dof = 3e6;                        % max # of degrees of freedom
    
    %% run the simluation and save the data
    datLoc = '.\test\cross\072324_sweep1\';
    P.datLoc = datLoc;
    bds = solveBands(P);
end 