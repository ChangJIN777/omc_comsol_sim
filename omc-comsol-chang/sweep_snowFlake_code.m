%% sweeping code 
P.celltype = 'snowflake';                   % specify the cell type

w = 100e-9; 
w_min = w*0.8;
w_max = w*1.3;
r_min = 200e-9;
r_max = 300e-9;
a_min = 680e-9;
a_max = 720e-9;
w_list = linspace(w_min,w_max,5);
r_list = linspace(r_min,r_max,5);
a_list = linspace(a_min,a_max,5);

%% run the lattice constatn sweep
for i=1:length(a_list)
    sweep_snowflake_lattice(P,a_list(i));    
end

%% run the hole size sweep
for i=1:length(w_list)
    for j=1:length(r_list)
        sweep_snowflake_holeSize(P,w_list(i),r_list(j));
    end
end
function sweep_snowflake_holeSize(P,w,r)    
    clc; close all;
    
    %% unit cell params 
    P.xsect = 'rect'; 
    P.beamMat = 'diamond';                  % beam material name
    P.celltype = 'snowflake';                   % specify the cell type
    P.unitcell = 'hexagonal';                  % specify the shape of the unit cell
    P.a = 700e-9;              % lattice constant 
    P.w = w;              % unit  cell width (along x)
    P.r = r;              % unit cell height (along y)
    P.th = 300e-9;             % height (along x) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r1 = 10e-9;             % width (along y) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r2 = 10e-9;              % height (along x) of each leg in cross (for celltype = 'hollow')
                                % or of outer fins (for celltype = 'solid')
    P.nperiod = 1;  % no. of periods to simulate for
    P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
    P.mbevenz = 1;      % 1 to find even mechanical mode about z
    
    P.kpts = 7;                             % no. of k-points, EXCLUDING gamma point
    P.nbands = 25;                           % no. of bands to solve for
    
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
    P.TwoSymPlanes = 0; 
    P.zSymCondition = 1;
    P.mbeveny = 0;                          % 1 to find even mechanical mode about y
    P.mbevenz = 1;                          % 1 to find even mechanical mode about z
    P.freq = 0;                             % target frequency - set to 0 for bandstructure simulations
    P.meshSize = 7;                         % mesh quality for mechanical simulations
    P.fixed_bc = 0;                       % 1 to fixed the boundaries for xz planes at y = +/- w/2
    
    P.anisoMat = 1;
    P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                            % from <100> inplane direction about <100> surface normal
    
    %% define the maximum number of degree of freedom to limit the simulation time
    P.max_dof = 3e6;                        % max # of degrees of freedom
    
    % %% debugging the unit cells 
    % % import COMSOL class
    % import com.comsol.model.*
    % import com.comsol.model.util.*
    % 
    % ModelUtil.showProgress(true);
    % ModelUtil.clear();
    % clear model
    % 
    % % create COMSOL model named 'model' from which COMSOL methods can be called, 
    % % e.g. model.save
    % model = ModelUtil.create('model');
    % 
    % buildBoomerangUnitCell(model,P);
    % mphlaunch(model);
    %% Single solve
    currentDate = datestr(now,'mmddyyyy');
    datLoc = ['.\test\snowflake_sweep\',currentDate,'\'];
    P.datLoc = datLoc;
    bds = solveBands(P);
end 

function sweep_snowflake_lattice(P,a)    
    clc; close all;
    
    %% unit cell params 
    P.xsect = 'rect'; 
    P.beamMat = 'diamond';                  % beam material name
    P.celltype = 'snowflake';                   % specify the cell type
    P.unitcell = 'hexagonal';                  % specify the shape of the unit cell
    P.a = a;              % lattice constant 
    P.w = 110e-9;              % unit  cell width (along x)
    P.r = 300e-9;              % unit cell height (along y)
    P.th = 300e-9;             % height (along x) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r1 = 10e-9;             % width (along y) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r2 = 10e-9;              % height (along x) of each leg in cross (for celltype = 'hollow')
                                % or of outer fins (for celltype = 'solid')
    P.nperiod = 1;  % no. of periods to simulate for
    P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
    P.mbevenz = 1;      % 1 to find even mechanical mode about z
    
    P.kpts = 7;                             % no. of k-points, EXCLUDING gamma point
    P.nbands = 25;                           % no. of bands to solve for
    
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
    P.TwoSymPlanes = 0; 
    P.zSymCondition = 1;
    P.mbeveny = 0;                          % 1 to find even mechanical mode about y
    P.mbevenz = 1;                          % 1 to find even mechanical mode about z
    P.freq = 0;                             % target frequency - set to 0 for bandstructure simulations
    P.meshSize = 7;                         % mesh quality for mechanical simulations
    P.fixed_bc = 0;                       % 1 to fixed the boundaries for xz planes at y = +/- w/2
    
    P.anisoMat = 1;
    P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                            % from <100> inplane direction about <100> surface normal
    
    %% define the maximum number of degree of freedom to limit the simulation time
    P.max_dof = 3e6;                        % max # of degrees of freedom
    
    % %% debugging the unit cells 
    % % import COMSOL class
    % import com.comsol.model.*
    % import com.comsol.model.util.*
    % 
    % ModelUtil.showProgress(true);
    % ModelUtil.clear();
    % clear model
    % 
    % % create COMSOL model named 'model' from which COMSOL methods can be called, 
    % % e.g. model.save
    % model = ModelUtil.create('model');
    % 
    % buildBoomerangUnitCell(model,P);
    % mphlaunch(model);
    %% Single solve
    currentDate = datestr(now,'mmddyyyy');
    datLoc = ['.\test\snowflake_sweep\',currentDate,'\'];
    P.datLoc = datLoc;
    bds = solveBands(P);
end 