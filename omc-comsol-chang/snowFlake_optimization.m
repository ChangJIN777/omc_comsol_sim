clear, clc;
% initial parameters
a0 = 220e-9;
w0 = 40e-9;
r0 = 90e-9;
params0 = [a0,w0,r0];
lower_bound = [200e-9,30e-9,40e-9];
upper_bound = [300e-9,80e-9,120e-9];

% run the optimization code 
options = optimset('PlotFcns',@optimplotfval);
func = @(params) runOptimization(params);
params = fminsearch(func,params0,options);

function [fullMidBand, fullGapSize] = run_snowFlake(P,a,w,r)    
    %% unit cell params 
    P.xsect = 'rect'; 
    P.beamMat = 'diamond';                  % beam material name
    P.celltype = 'snowflake';                   % specify the cell type
    P.unitcell = 'hexagonal';                  % specify the shape of the unit cell
    P.a = a;              % lattice constant 
    P.w = w;              % unit cell width (along x)
    P.r = r;              % unit cell height (along y)
    P.th = 100e-9;             % height (along x) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r1 = 10e-9;             % width (along y) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r2 = 10e-9;              % height (along x) of each leg in cross (for celltype = 'hollow')
                                % or of outer fins (for celltype = 'solid')
    P.nperiod = 1;  % no. of periods to simulate for
    P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
    P.mbevenz = 1;      % 1 to find even mechanical mode about z
    
    P.kpts = 7;                             % no. of k-points, EXCLUDING gamma point
    P.nbands = 18;                           % no. of bands to solve for
    
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
    datLoc = '.\test\snowflake\090224_sweep1\';
    P.datLoc = datLoc;
    bds = solveBands(P);
    
    %% output the midgap and the bandgap of the complete bandgaps
    fullBand = bds.full;
    fullMidBand = fullBand.midGap;
    fullGapSize = fullBand.gapSize;
    
end 

function fitness = calFitnessSnow(fullMidBand,fullGapSize)
    targetFreq = 50e9; % the target frequency of the bandgap 
    freqTolerance = 1e9; % the tolerance for frequency offsets 
    % the fitness function associated with the optimization code
    fitness = -fullGapSize.*exp(-((fullMidBand-targetFreq)./freqTolerance).^2);
end

function fitness = runOptimization(params)
    P.celltype = 'snowflake';                   % specify the cell type
    a = params(1);
    w = params(2);
    r = params(3);
    
    %% account for any fabrication intolerance 
    smallestFeature = (a/2-r)-(w/2)/sqrt(3);
    if smallestFeature < 10e-9
        fitness = 0;
    else
        [fullMidBand, fullGapSize] = run_snowFlake(P,a,w,r);
        fitness_list = calFitnessSnow(fullMidBand,fullGapSize);
        fitness = min(fitness_list);
    end
end