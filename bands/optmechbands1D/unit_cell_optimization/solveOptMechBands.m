function results = solveOptMechBands(vars,datLoc)
%SOLVEOPTMECHBANDS Summary of this function goes here
%   Detailed explanation goes here

%% Unit cell parameters
P.a = vars.a*1e-9;%580e-9;                           % nominal lattice constant
P.w = vars.w*1e-9;%929e-9;                           % beam width
P.hx = vars.hx*1e-9;%250e-9;                          % nominal hole height (along x-axis)
P.hy = vars.hy*1e-9;%590e-9;                          % nominal hole width (along y-axis

P.datLoc = datLoc;
% Constants
P.xsect = 'rect';%'tri';                       % beam cross sectional shape  - 'tri' or 'rect'
P.beamMat = 'diamond';                  % beam material name
P.anisoMat = 1;
P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                        % from <100> inplane direction about <100> surface normal
                                        
P.theta = 35;                           % etch angle in degrees (will not apply for P.xsect = 'rect')
P.th = 160e-9;%P.w/(2*tan(P.theta*pi/180));                          % beam thickness (define using P.theta for P.xsect = 'tri') 
P.hole = 'elps';                        % hole shape
P.nperiod = 1;                          % no. of periods to simulate for
P.holeatedge = 1;                       % 1/0 for hole at edge/center of unit cell
                                 


%% Optical simulation parameters
ofem.doOpt = 1;

ofem.obfreq = 100e12;                      % target frequency - set to 150THz for bands in telecom
ofem.kpts = 21;                            % no. of k-points, EXCLUDING gamma point
ofem.nbands = 4;                           % no. of bands to solve for
ofem.max_dof = 3e6;                        % max # of degrees of freedom
ofem.meshSize = 5;                         % Default meshsize

ofem.airrad = 1550e-9;
ofem.anisoMat = 1;
ofem.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                        % from <100> inplane direction about <100> surface normal
ofem.obeveny = -1;                          % -1 for TE-like
ofem.obevenz = 0;                          % 1 to find even mechanical mode about z
ofem.solveTM = 0;                          % 1 to solve for TM modes (if running solveBands)
ofem.completeBandGaps = 0;                 % 1 to plot complete bandgaps (across all symmetries)

ofem.plotgeom = 0;                         % 1 to plot the geometry
ofem.savedat = 0;                          % 1 to save data structures
ofem.savebndplot = 1;                      % 1 to save bandstructure plot
ofem.saveplots = 0;                        % 1 to save displacement and strain profiles

%% Mechanical simulation parameters
mfem.doMech = 0;

mfem.mbfreq = 0;                           % target frequency - set to 0 for bandstructure simulations
mfem.kpts = 5;                             % no. of k-points, EXCLUDING gamma point
mfem.nbands = 5;                           % no. of bands to solve for
mfem.max_dof = 3e6;                        % max # of degrees of freedom
mfem.meshSize = 5;

mfem.mbeveny = 1;                          % 1 to find even mechanical mode about y
mfem.mbevenz = 0;                          % 1 to find even mechanical mode about z
mfem.solveasym = 0;                        % 1 to solve for antisymmetric bands
mfem.completeBandGaps = 0;                 % 1 to plot complete bandgaps (across all symmetries)

mfem.plotgeom = 0;                         % 1 to plot the geometry
mfem.savedat = 0;                          % 1 to save data structures
mfem.savebndplot = 1;                      % 1 to save bandstructure plot
mfem.saveplots = 0;                        % 1 to save displacement and strain profiles

if ofem.doOpt
    tic
    P.optmech = 1;                      % 1/0 for optical/mechanical simulation
    disp('Starting optical simulation')
    opt_bds = solveOpticalBands(P,ofem);
    toc
    
    disp('Optical Bandgap properties')
    disp(['Midgap frequency: ', num2str(opt_bds.midGap/1e12), ' THz']);
    disp(['Gap size: ', num2str(100*opt_bds.gapSize/opt_bds.midGap), '%']);
    
    results.opt_bds = opt_bds;
end

if mfem.doMech
    tic
    P.optmech = 0;
    disp('Starting mechanical simulation')
    mech_bds = solveMechBands(P,mfem);
    toc
    
    disp('Mechanical Bandgap properties')
    disp(['Midgap frequency: ', num2str(mech_bds.midGap/1e9), ' GHz']);
    disp(['Gap size: ', num2str(100*mech_bds.gapSize/mech_bds.midGap), '%']);
    
    results.mech_bds = mech_bds;
end

end

