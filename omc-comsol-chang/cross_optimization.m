%% optimization code 
clear, clc;
% initial parameters
a0 = 200e-9;
w0 = 50e-9;
h0 = 150e-9;
lower_bound = [290e-9,30e-9,40e-9];
upper_bound = [400e-9,80e-9,120e-9];

% run the optimization code with fixed lattice constant 
fix_lattice = 0;

% run the optimization code 
options = optimset('PlotFcns',@optimplotfval);

if fix_lattice==1
    params0 = [w0,h0];
    func = @(params) runOptimization_fixedLattice(params);
    params = fminsearch(func,params0,options);
else
    params0 = [a0,w0,h0];
    func = @(params) runOptimization(params);
    params = fminsearch(func,params0,options);
end

% %% manual sweep 
% P.celltype = 'snowflake';                   % specify the cell type
% 
% a0 = 260e-9;
% w0 = 50e-9;
% r0 = 100e-9;
% r_list = linspace(90e-9,120e-9,4);
% for i=1:length(r_list)
%     run_snowFlake(P,a0,w0,r_list(i));
% end

%% function definition 
function [fullMidBand, fullGapSize] = run_Cross(P,a,w,h)    
    %% unit cell params 
    P.xsect = 'rect'; 
    P.beamMat = 'diamond';                  % beam material name
    P.celltype = 'cross';                   % specify the cell type
    P.unitcell = 'square';                  % specify the shape of the unit cell
    P.a = a;              % lattice constant 
    P.w = w;              % unit cell width (along x)
    P.h = h;              % unit cell height (along y)
    P.th = 120e-9;             % height (along x) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r1 = 10e-9;             % width (along y) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r2 = 10e-9;              % height (along x) of each leg in cross (for celltype = 'hollow')
                                % or of outer fins (for celltype = 'solid')
    P.nperiod = 1;  % no. of periods to simulate for
    P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
    P.mbevenz = 1;      % 1 to find even mechanical mode about z
    
    P.kpts = 5;                             % no. of k-points, EXCLUDING gamma point
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
    datLoc = '.\test\cross\092624_sweep1\';
    P.datLoc = datLoc;
    bds = solveBands(P);
    
    %% output the midgap and the bandgap of the complete bandgaps
    fullBand = bds.full;
    fullMidBand = fullBand.midGap;
    fullGapSize = fullBand.gapSize;
    
end 


function fitness = calFitnessSnow(fullMidBand,fullGapSize)
    targetFreq = 50e9; % the target frequency of the bandgap 
    freqTolerance = 2e9; % the tolerance for frequency offsets 
    % the fitness function associated with the optimization code
    fitness = -fullGapSize.*exp(-((fullMidBand-targetFreq)./freqTolerance).^2);
end

function fitness = runOptimization(params)
    P.celltype = 'cross';                   % specify the cell type
    a = params(1);
    w = params(2);
    h = params(3);
    %% account for any fabrication intolerance 
    smallestFeature = (a-h)/2;
    if smallestFeature < 20e-9
        fitness = 0;
    else
        [fullMidBand, fullGapSize] = run_Cross(P,a,w,h);
        fitness_list = calFitnessSnow(fullMidBand,fullGapSize);
        fitness = min(fitness_list);
    end
end