% front panel script for nanobeam FEM simulations (example for rect cross-section, beam-hole geometry)
clear all; close all; clc
%% initial parameters (v1)
% wo = 557e-9;
% wi = 200e-9;
% ho = 235e-9;
% hi = 175e-9;
% wo_ctr = 478e-9;
% wi_ctr = 300e-9;
% ho_ctr = 257e-9;
% hi_ctr = 175e-9;
% params0 = [wo,wi,ho,hi,wo_ctr,wi_ctr,ho_ctr,hi_ctr];
% minfitness = runSnowflakeFEM(params0);
% fprintf('the minimum fitness %.2f',minfitness);
% %% initial parameters (v2)
wo = 557e-9;
ho = 235e-9;
wo_ctr = 478e-9;
ho_ctr = 257e-9;
params0 = [wo,ho,wo_ctr,ho_ctr];
% % minfitness = runSnowflakeFEM(params0);
% % % fprintf('the minimum fitness %.2f',minfitness);
%% run the optimization code 
options = optimset('PlotFcns',@optimplotfval);
func = @(params) runSnowflakeFEM_v2(params);
params = fminsearch(func,params0);
% minfitness = runSnowflakeFEM(params0);
% fprintf('the minimum fitness %.2f',minfitness);

function minFitness = runSnowflakeFEM(params)    
    %% optimization file name 
    currentDate = datestr(now,'mmddyyyy');
    optFilename = ['.\testing\optimizeData_',currentDate,'.csv'];
    %% geometry parameters
    % unit cell params
    P.xsect = 'rect';                        % beam cross sectional shape - 'tri' or 'rect'
    P.beamMat = 'diamond';                  % beam material name
    P.celltype = 'snowflake';                   % specify the cell type
    P.anisoMat = 1;
    
    % unit cell geometry
    P.a = 650e-9;              % lattice constant 
    P.w = 80e-9;              % unit cell width (along x)
    P.r = 250e-9;              % unit cell height (along y)
    P.th = 350e-9;             % height (along x) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.wo = params(1);           % the height of the hole in the lower portion
    P.wi = params(2);           % the width of the hole in the lower portion                            
    P.ho = params(3);
    P.hi = params(4);
    P.d = 100e-9;    % 
    P.b = sqrt(3)*P.a/2; 
    P.r1 = 10e-9;             % width (along y) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r2 = 10e-9;              % height (along x) of each leg in cross (for celltype = 'hollow')
                                % or of outer fins (for celltype = 'solid')
    % hole params for symmetric cavity / right half of asymmetric cavity
    P.MN_left = 15;                         % # holes in the left mirror region 
    P.MN_right = 15;                        % # holes in the right mirror region 
    P.TN = 5;                               % # holes in the tapering defect region
    
    % cavity taper params
    P.holeatctr = 1;                        % 1/0 for hole/dielectric in middle
    P.taperFunc = 'cubic';                  % linear/cubic/quadratic taper function to center hole in cavity
    P.taperTo = 'custom';                 % taper to custom hole in center of cavity; disable for taper to maxdef
    P.a_ctr = 650e-9;                     % for taperTo = 'custom': lattice constant of center hole
    P.wo_ctr = params(5);                    % for taperTo = 'custom': hole width of center hole
    P.wi_ctr = params(6);                    % for taperTo = 'custom': hole width of center hole
    P.ho_ctr = params(7);                    % for taperTo = 'custom': hole height of center hole
    P.hi_ctr = params(8);
    % P.cavlen = 0e-9;                      % custom cavity length between two center holes; disable if not used
    
    % end waveguide mirror taper params
    P.wvgmir = 0;                           % no. of mirrors in end waveguide mirror taper
    P.wgmTaper.func = 'cubic';              % taper function - linear, quadratic, cubic
    P.wgmTaper.endtype = 'custom';          % taper to custom, maxdef, zero (or '') end hole
    P.wgmTaper.a_end = 500e-9;              % for endtype = 'custom': lattice constant of end hole
    P.wgmTaper.wo_end = P.wo;             % for endtype = 'custom': hole height of end hole
    P.wgmTaper.wi_end = P.wi;             % for endtype = 'custom': hole width of end hole
    P.wgmTaper.ho_end = P.ho;             % for endtype = 'custom': hole height of end hole
    P.wgmTaper.hi_end = P.hi;             % for endtype = 'custom': hole width of end hole
    % P.wgmTaper.maxdef_end = 0.094932;     % defect percentage - for endtype = 'maxdef'
    % P.wgmTaper.oblong_end = 2.9172;       % oblong parameter (zero if holes are not changed) - for endtype = 'maxdef'
    
    % asymmetric cavity - specify param data struct P.PL with similar fields to
    % P, for left half of asymmetric cavity
    P.asymCav = 0;                          % 1 to enable asymmetric cavity
    if P.asymCav                
        P.PL = P;
        P.PL.nholes = 3+P.ndef;
        P.PL.wvgmir = 5;
    end
    
    P.lambda = 1553e-9;                     % target optical wavelength
    
    % Disorder
    P.stdDev = [0,0];                       % standard deviation of hole dimensions (hh,hw)
    P.stdDevPos = 0;                        % standard deviation of hole positions
    P.asym = 0;                             % cross-section asymmetry (target y-offset in bottom apex position)
    %% impose lower and upper bounds 
    outterFrameParams = [P.wo,P.ho,P.wo_ctr,P.ho_ctr];
    lowerbounds = [100e-9,100e-9,100e-9,100e-9];
    upperbounds = [P.a-100e-9,P.b-P.r*sqrt(3)/2-50e-9,P.a-100e-9,P.b-P.r*sqrt(3)/2-50e-9];
    ifViolateUpperBound = outterFrameParams > upperbounds;
    ifViolateLowerBound = outterFrameParams < lowerbounds;
    if sum(ifViolateLowerBound)>0 || sum(ifViolateUpperBound)>0
        minFitness = 1000;
        return
    end
    % make sure that the structures are fab complient
    if P.wo-P.wi<140e-9 || P.ho-P.hi<50e-9 || P.wo_ctr-P.wi_ctr<140e-9 || P.ho_ctr-P.hi_ctr<50e-9 
        minFitness = 1000;
        return
    end
    %% specify simulation/calculation/plot/save options
    P.solveMech = 1;                        % 1 to solve for mechanics
    P.solveOpt = 1;                         % 1 to solve for optics
    P.nbeam = 2.4028;        % refractive index of the dielectric material
    P.calcG = 1*(P.solveMech && P.solveOpt);% 1 to calculate optomechanical coupling
    P.calcS = 0*P.solveMech;                % 1 to calculate strain coupling
    P.solveMechPML = 0;                     % 1 to solve for mechanical Q (future implementation)
    
    % plotting & saving
    P.plotgeom = 0;                         % 1 to plot the geometry
    P.storeMPH = 0;                         % 1 to save COMSOL model file
    P.plotMech = 0*P.solveMech;             % 1 to plot displacement and strain profiles
    P.plotOpt = 0*P.solveOpt;               % 1 to plot E-field profiles
    P.plotStrCpl = 0*P.calcS;               % 1 to plot strain coupling profile
    
    %% mechanical simulation parameters 
    % solid mechanics solver parameters
    P.mevenx = 1;                           % +/-1 to find even/odd mode about x; 0 for fixed BC
    P.meveny = 1;                           % +/-1 to find even/odd mode about y
    P.mevenz = 1;                           % +/-1 to find even/odd mode about z
    P.freq = 20.5e9;                           % target mechanical frequency
    P.mneigs = 10;                          % # of eignevalues to find
    P.mMesh = 7;                            % mesh quality for mechanical simulations
    P.mAdjMesh = 1;                         % adjust mesh if DOFs exceed max_dof
    % if we are adding the 2D phononic shield
    P.addshield = 0;
    
    % rotate crystal orientation of elasticity matrix
    % ccw rotation in deg from <100> inplane direction about <100> surface normal
    P.rxtal = 45;
    P.rxtalInFilename = 1;
    
    %% optical simulation parameters
    % rf module solver parameters
    P.oevenx = 1;            % +/-1 to find even/odd optical mode about x (-1 == fundamental for hole in center)
    P.oeveny = -1;                          % +/-1 to find even/odd optical mode about y (-1 == TE-like)
    P.oevenz = 1;                           % +/-1 to find even/odd optical mode about z 
    P.oneigs = 10;                           % # of eigenvalues to find
    P.oMesh = 7;                            % mesh quality for optical simulations
    P.oAdjMesh = 1;                         % adjust mesh if DOFs exceed max_dof
    P.airrad = P.lambda+P.a*sqrt(3)*3;            % radius of air cylinder surrounding nanobeam
    
    %% OM coupling parameters
    P.g0min = 80e3;                         % min g0 above which to save plots for
    
    %% SiV strain coupling: susceptibilities and positions
    % specify SiV axis - [1 1 1],[-1 1 1],[-1 -1 1],[1 -1 1]
    % for multiple axes, specify as arrays in cell
    P.zSiV = {[1 1 1],[-1 1 1],[-1 -1 1],[1 -1 1]};
    
    % specify custom coordinates to plot at
    P.xSlc = 0; % relative to center
    P.ySlc = 0;
    P.zSlc = P.th/2 - 50e-9;
    
    % define aperture to run stats (min, max, mean, stddev) on strain coupling
    P.LStats.xmin = []; % leave empty to define based on generated geometry
    P.LStats.xmax = [];
    P.LStats.ymin = 0;
    P.LStats.ymax = 60e-9;
    P.LStats.zmin = P.th/2-80e-9;   % z-coords relative to center of beam
    P.LStats.zmax = P.th/2;
    
    %% Mechanical PML simulation settings (future implementation)
    P.PMLmesh = 7;
    P.PMLmeshDiv = 20;
    P.PMLLen = 10e-6;
    P.PMLstr = 0.008;
    
    %% Simulation settings
    P.max_dof = 5e6;                        % max # of degrees of freedom
    
    % single run
    % data location to save files in
    datLoc = '.\testing'; 
    [ds,model] = RunNanobeamFEM(P,datLoc);
    
    % extract the Q of the optical mode and gOM 
    ofem = ds.ofem;
    cpl = ds.cpl;
    opt_Q = cpl.Q; % the Q factor of the optical mode 
    mech_freq = cpl.mechFreq; % the mechanical frequency 
    opt_lambda = cpl.optWvl; % the wavelength of the corresponding resonances 
    gOM = cpl.gOMmax;
    
    % calculate and find the minimum fitness
    fitness = calFitness(opt_Q,opt_lambda,gOM);
    [minFitness,idx] = min(fitness);
       
    % record the data 
    data = [P.a,P.w,P.r,P.ho,P.hi,P.wo,P.wi,P.d,P.ho_ctr,...
        P.hi_ctr,P.wo_ctr,P.wi_ctr,opt_Q,opt_lambda,mech_freq,gOM,minFitness];
    writematrix(data,optFilename,'WriteMode','append');

end

function minFitness = runSnowflakeFEM_v2(params)    
    %% optimization file name 
    currentDate = datestr(now,'mmddyyyy');
    optFilename = ['.\testing\optimizeData_',currentDate,'.csv'];
    %% impose lower and upper bounds 
    lowerbounds = [450e-9,220e-9,300e-9,250e-9,];
    upperbounds = [600e-9,500e-9,600e-9,500e-9];
    ifViolateUpperBound = params > upperbounds;
    ifViolateLowerBound = params < lowerbounds;
    if sum(ifViolateLowerBound)>0 || sum(ifViolateUpperBound)>0
        minFitness = 1000;
        return
    end
    %% geometry parameters
    % unit cell params
    P.xsect = 'rect';                        % beam cross sectional shape - 'tri' or 'rect'
    P.beamMat = 'diamond';                  % beam material name
    P.celltype = 'snowflake';                   % specify the cell type
    P.anisoMat = 1;
    
    % unit cell geometry
    P.a = 650e-9;              % lattice constant 
    P.w = 80e-9;              % unit cell width (along x)
    P.r = 250e-9;              % unit cell height (along y)
    P.th = 350e-9;             % height (along x) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.wo = params(1);           % the height of the hole in the lower portion
    P.wi = 200e-9;           % the width of the hole in the lower portion                            
    P.ho = params(2);
    P.hi = 175e-9;
    P.d = 100e-9;    % 
    P.b = sqrt(3)*P.a/2; 
    P.r1 = 10e-9;             % width (along y) of cross (for celltype = 'hollow')
                                % or of inner block (for celltype = 'solid')
    P.r2 = 10e-9;              % height (along x) of each leg in cross (for celltype = 'hollow')
                                % or of outer fins (for celltype = 'solid')
    % hole params for symmetric cavity / right half of asymmetric cavity
    P.MN_left = 15;                         % # holes in the left mirror region 
    P.MN_right = 15;                        % # holes in the right mirror region 
    P.TN = 5;                               % # holes in the tapering defect region
    
    % cavity taper params
    P.holeatctr = 1;                        % 1/0 for hole/dielectric in middle
    P.taperFunc = 'cubic';                  % linear/cubic/quadratic taper function to center hole in cavity
    P.taperTo = 'custom';                 % taper to custom hole in center of cavity; disable for taper to maxdef
    P.a_ctr = 650e-9;                     % for taperTo = 'custom': lattice constant of center hole
    P.wo_ctr = params(3);                    % for taperTo = 'custom': hole width of center hole
    P.wi_ctr = 300e-9;                    % for taperTo = 'custom': hole width of center hole
    P.ho_ctr = params(4);                    % for taperTo = 'custom': hole height of center hole
    P.hi_ctr = 175e-9;
    % P.cavlen = 0e-9;                      % custom cavity length between two center holes; disable if not used
    
    % end waveguide mirror taper params
    P.wvgmir = 0;                           % no. of mirrors in end waveguide mirror taper
    P.wgmTaper.func = 'cubic';              % taper function - linear, quadratic, cubic
    P.wgmTaper.endtype = 'custom';          % taper to custom, maxdef, zero (or '') end hole
    P.wgmTaper.a_end = 650e-9;              % for endtype = 'custom': lattice constant of end hole
    P.wgmTaper.wo_end = P.wo;             % for endtype = 'custom': hole height of end hole
    P.wgmTaper.wi_end = P.wi;             % for endtype = 'custom': hole width of end hole
    P.wgmTaper.ho_end = P.ho;             % for endtype = 'custom': hole height of end hole
    P.wgmTaper.hi_end = P.hi;             % for endtype = 'custom': hole width of end hole
    % P.wgmTaper.maxdef_end = 0.094932;     % defect percentage - for endtype = 'maxdef'
    % P.wgmTaper.oblong_end = 2.9172;       % oblong parameter (zero if holes are not changed) - for endtype = 'maxdef'
    
    % asymmetric cavity - specify param data struct P.PL with similar fields to
    % P, for left half of asymmetric cavity
    P.asymCav = 0;                          % 1 to enable asymmetric cavity
    if P.asymCav                
        P.PL = P;
        P.PL.nholes = 3+P.ndef;
        P.PL.wvgmir = 5;
    end
    
    P.lambda = 1553e-9;                     % target optical wavelength
    
    % Disorder
    P.stdDev = [0,0];                       % standard deviation of hole dimensions (hh,hw)
    P.stdDevPos = 0;                        % standard deviation of hole positions
    P.asym = 0;                             % cross-section asymmetry (target y-offset in bottom apex position)
    
    %% specify simulation/calculation/plot/save options
    P.solveMech = 1;                        % 1 to solve for mechanics
    P.solveOpt = 1;                         % 1 to solve for optics
    P.nbeam = 2.4028;        % refractive index of the dielectric material
    P.calcG = 1*(P.solveMech && P.solveOpt);% 1 to calculate optomechanical coupling
    P.calcS = 0*P.solveMech;                % 1 to calculate strain coupling
    P.solveMechPML = 0;                     % 1 to solve for mechanical Q (future implementation)
    
    % plotting & saving
    P.plotgeom = 0;                         % 1 to plot the geometry
    P.storeMPH = 0;                         % 1 to save COMSOL model file
    P.plotMech = 0*P.solveMech;             % 1 to plot displacement and strain profiles
    P.plotOpt = 0*P.solveOpt;               % 1 to plot E-field profiles
    P.plotStrCpl = 0*P.calcS;               % 1 to plot strain coupling profile
    
    %% mechanical simulation parameters 
    % solid mechanics solver parameters
    P.mevenx = 1;                           % +/-1 to find even/odd mode about x; 0 for fixed BC
    P.meveny = 1;                           % +/-1 to find even/odd mode about y
    P.mevenz = 1;                           % +/-1 to find even/odd mode about z
    P.freq = 20.5e9;                           % target mechanical frequency
    P.mneigs = 10;                          % # of eignevalues to find
    P.mMesh = 7;                            % mesh quality for mechanical simulations
    P.mAdjMesh = 1;                         % adjust mesh if DOFs exceed max_dof
    % if we are adding the 2D phononic shield
    P.addshield = 0;
    
    % rotate crystal orientation of elasticity matrix
    % ccw rotation in deg from <100> inplane direction about <100> surface normal
    P.rxtal = 45;
    P.rxtalInFilename = 1;
    
    %% optical simulation parameters
    % rf module solver parameters
    P.oevenx = 1;            % +/-1 to find even/odd optical mode about x (-1 == fundamental for hole in center)
    P.oeveny = -1;                          % +/-1 to find even/odd optical mode about y (-1 == TE-like)
    P.oevenz = 1;                           % +/-1 to find even/odd optical mode about z 
    P.oneigs = 20;                           % # of eigenvalues to find
    P.oMesh = 7;                            % mesh quality for optical simulations
    P.oAdjMesh = 1;                         % adjust mesh if DOFs exceed max_dof
    P.airrad = P.lambda+P.a*sqrt(3)*3;            % radius of air cylinder surrounding nanobeam
    
    %% OM coupling parameters
    P.g0min = 80e3;                         % min g0 above which to save plots for
    
    %% SiV strain coupling: susceptibilities and positions
    % specify SiV axis - [1 1 1],[-1 1 1],[-1 -1 1],[1 -1 1]
    % for multiple axes, specify as arrays in cell
    P.zSiV = {[1 1 1],[-1 1 1],[-1 -1 1],[1 -1 1]};
    
    % specify custom coordinates to plot at
    P.xSlc = 0; % relative to center
    P.ySlc = 0;
    P.zSlc = P.th/2 - 50e-9;
    
    % define aperture to run stats (min, max, mean, stddev) on strain coupling
    P.LStats.xmin = []; % leave empty to define based on generated geometry
    P.LStats.xmax = [];
    P.LStats.ymin = 0;
    P.LStats.ymax = 60e-9;
    P.LStats.zmin = P.th/2-80e-9;   % z-coords relative to center of beam
    P.LStats.zmax = P.th/2;
    
    %% Mechanical PML simulation settings (future implementation)
    P.PMLmesh = 7;
    P.PMLmeshDiv = 20;
    P.PMLLen = 10e-6;
    P.PMLstr = 0.008;
    
    %% Simulation settings
    P.max_dof = 5e6;                        % max # of degrees of freedom
    
    % single run
    % data location to save files in
    datLoc = '.\testing'; 
    [ds,model] = RunNanobeamFEM(P,datLoc);
    
    % extract the Q of the optical mode and gOM 
    ofem = ds.ofem;
    cpl = ds.cpl;
    opt_Q = cpl.Q; % the Q factor of the optical mode 
    mech_freq = cpl.mechFreq; % the mechanical frequency 
    opt_lambda = cpl.optWvl; % the wavelength of the corresponding resonances 
    gOM = cpl.gOMmax;
    
    if isempty(opt_lambda)
        minFitness = 1000;
        return;
    end
    % calculate and find the minimum fitness
    fitness = calFitness(opt_Q,opt_lambda,gOM);
    [minFitness,idx] = min(fitness);
       
    % record the data 
    data = [P.a,P.w,P.r,P.ho,P.hi,P.wo,P.wi,P.d,P.ho_ctr,...
        P.hi_ctr,P.wo_ctr,P.wi_ctr,opt_Q,opt_lambda,mech_freq,gOM,minFitness];
    writematrix(data,optFilename,'WriteMode','append');

end

function fitness = calFitness(opt_QAll,opt_lambdaAll,gOM)
    wavelength_tolerance = 100e-9;
    fitness = -abs(gOM.*opt_QAll.*exp(-((opt_lambdaAll-1550e-9)./(wavelength_tolerance)).^2));
end