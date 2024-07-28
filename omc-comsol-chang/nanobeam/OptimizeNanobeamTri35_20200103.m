% Function to optimize nanobeam design for OM coupling or SiV strain
% datLoc must end with \
% Supplementary functions...

% datLoc = 'L:\Individuals\cchia\OMC\COMSOL\nanobeamOMscripts\optWithNomCubicTaper\';
% if ~exist(datLoc,'dir')
%     mkdir(datLoc)
% end

function pp = OptimizeNanobeamTri35_20200103(datLoc)
clear global P
clear global optimds
clear global dtID
%% folders
%append \ to end of datLoc if not present
if ~strcmp(datLoc(end),'\')
    datLoc = [datLoc,'\'];
end
if ~exist(datLoc,'dir')
    mkdir(datLoc)
end

global P
global file
global optimds
global dtID
singlePass = 0;
randomStart = 1; % 1 if we want to pick a random starting point
strPoint0 = 0; % last starting point for which iteration data was saved

maxSPtoRun = 0;     % no. of starting points to run for; 0 to run indefinitely

%% Nominal nanobeam parameters
% ===geometry parameters
% unit cell params
P.xsect = 'tri';                        % beam cross sectional shape - 'tri' or 'rect'
P.beamMat = 'diamond';                  % beam material name
P.anisoMat = 1;

P.a = 565.7e-9;                         % nominal lattice constant
P.w = 1074.8e-9;                        % beam width
P.theta = 35;                           % etch angle in degrees (no effect for rect cross section)
P.th = P.w/(2*tan(P.theta*pi/180));     % beam thickness
P.hx = 226.3e-9;                        % nominal hole height (along x-axis)
P.hy = 698.6e-9;                        % nominal hole width (along y-axis)

% unit cell bands
% P.bandLoc = 'L:\Individuals\cchia\OMC\45oDesign\mechBands\tri_elps_45o_rxt_45o_a_566nm_w_1075nm_hx_226nm_hy_699nm_fullBands.fig';
% P.bandmidgap = 9.26e9;

% hole params for symmetric cavity / right half of asymmetric cavity
P.nholes = 18;                          % # holes in 1/2 beam length
P.ndef = 8;                             % # of holes in 1/2 defect region
P.maxdef = 0.15906;                     % defect percentage
P.oblong = 3.0563;                      % oblong parameter (zero if holes are not changed)

% cavity taper params
P.holeatctr = 1;                        % 1/0 for hole/dielectric in middle
P.taperFunc = 'cubic';                  % linear/cubic/quadratic taper function to center hole in cavity
% P.taperTo = 'custom';                 % taper to custom hole in center of cavity; disable for taper to maxdef
% P.a_ctr = 392e-9;                     % for taperTo = 'custom': lattice constant of center hole
% P.hx_ctr = 189e-9;                    % for taperTo = 'custom': hole height of center hole
% P.hy_ctr = 177e-9;                    % for taperTo = 'custom': hole width of center hole
% P.cavlen = 0e-9;                      % custom cavity length between two center holes; disable if not used

% end waveguide mirror taper params
P.wvgmir = 0;                           % no. of mirrors in end waveguide mirror taper
P.wgmTaper.func = 'cubic';              % taper function - linear, quadratic, cubic
P.wgmTaper.endtype = 'custom';          % taper to custom, maxdef, zero (or '') end hole
P.wgmTaper.a_end = 435e-9;              % for endtype = 'custom': lattice constant of end hole
P.wgmTaper.hx_end = 179e-9;             % for endtype = 'custom': hole height of end hole
P.wgmTaper.hy_end = 272e-9;             % for endtype = 'custom': hole width of end hole
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

P.lambda = 1550e-9;                      % target optical wavelength

%Disorder
P.stdDev = [0,0]; %0.05*[P.hh,P.hw];             % standard deviation of hole dimensions (hh,hw)
P.stdDevPos = 0; %0.05*P.a;                  % standard deviation of hole positions
P.asym = 0; %0.02*P.w;                       % cross-section asymmetry (target y-offset in bottom apex position)

%% specify simulation/calculation/plot/save options
P.solveMech = 1;                        % 1 to solve for mechanics
P.solveOpt = 1;                         % 1 to solve for optics
P.calcG = 1*(P.solveMech && P.solveOpt);% 1 to calculate optomechanical coupling
P.calcS = 1*P.solveMech;                % 1 to calculate strain coupling
P.solveMechPML = 0;                     % 1 to solve for mechanical Q (future implementation)

% plotting & saving
P.plotgeom = 0;                         % 1 to plot the geometry
P.storeMPH = 1;                         % 1 to save COMSOL model file
P.plotMech = 1*P.solveMech;             % 1 to plot displacement and strain profiles
P.plotOpt = 1*P.solveOpt;               % 1 to plot E-field profiles
P.plotStrCpl = 1*P.calcS;               % 1 to plot strain coupling profile

%% mechanical simulation parameters 
% solid mechanics solver parameters
P.mevenx = 1;                           % +/-1 to find even/odd mode about x; 0 for fixed BC
P.meveny = 1;                           % +/-1 to find even/odd mode about y
P.mevenz = 1;                           % +/-1 to find even/odd mode about z
P.freq = 5e9;                           % target mechanical frequency
P.mneigs = 10;                          % # of eignevalues to find
P.mMesh = 3;                            % mesh quality for mechanical simulations
P.mAdjMesh = 1;                         % adjust mesh if DOFs exceed max_dof

% rotate crystal orientation of elasticity matrix
% ccw rotation in deg from <100> inplane direction about <100> surface normal
P.rxtal = 45;
P.rxtalInFilename = 1; 

%% optical simulation parameters
% rf module solver parameters
P.oevenx = -1^(P.holeatctr);            % +/-1 to find even/odd optical mode about x (-1 == fundamental for hole in center)
P.oeveny = -1;                          % +/-1 to find even/odd optical mode about y (-1 == TE-like)
P.oevenz = 1;                           % +/-1 to find even/odd optical mode about z 
P.oneigs = 1;                           % # of eigenvalues to find
P.oMesh = 4;                            % mesh quality for optical simulations
P.oAdjMesh = 1;                         % adjust mesh if DOFs exceed max_dof
P.airrad = 2*P.lambda+P.w/2;            % radius of air cylinder surrounding nanobeam

%% OM coupling parameters
P.g0min = 80e3;                         % min g0 above which to save plots for

%% SiV strain coupling: susceptibilities and positions
% specify SiV axis - [1 1 1],[-1 1 1],[-1 -1 1],[1 -1 1]
% for multiple axes, specify as arrays in cell
P.zSiV = {[1 1 1]};

% specify custom coordinates to plot at
P.xSlc = 0; % relative to center
P.ySlc = 0;
P.zSlc = P.th/2 - 50e-9;

% define aperture to run stats (min, max, mean, stddev) on strain coupling
P.LStats.xmin = []; % leave empty to define based on generated geometry
P.LStats.xmax = [];
P.LStats.ymin = 0;
P.LStats.ymax = 60e-9;
P.LStats.zmin = P.th/2-70e-9;   % z-coords relative to center of beam
P.LStats.zmax = P.th/2-20e-9;

% ===Simulation settings
P.max_dof = 5e6;                      % max # of degrees of freedom

%% bandstructure simulation settings
P.nperiod = 1;                          % no. of periods to simulate for
P.holeatedge = 1;                       % 1/0 for hole at edge/center of unit cell

P.kpts = 7;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 8;                           % no. of bands to solve for
P.hole = 'elps';                        % hole shape

P.solveasym = 1;                        % 1 to solve for antisymmetric bands
P.completeBandGaps = 0;                 % 1 to plot complete bandgaps (across all symmetries)
P.plotgeom = 0;                         % 1 to plot the geometry
P.savedat = 0;                          % 1 to save data structures
P.savebndplot = 1;                      % 1 to save bandstructure plot
P.saveplots = 0;                        % 1 to save displacement and strain profiles

P.mbeveny = 1;                          % 1 to find even mechanical mode about y
P.mbevenz = 0;                          % 1 to find even mechanical mode about z
P.mbfreq = 0;                           % target frequency - set to 0 for bandstructure simulations

%% define fitness function
optimds.fitStr = ['min(abs(omds.cpl.gMax.*omds.cpl.optWvl./299792458),2.0e5/1.96e14)./(2.0e5/1.96e14)',...
            '.*min([omds.cpl.Q,1e6])./1e6',...
            '.*real(omds.cpl.SiV_111.full.LSiVAperMax)./1e6'];
optimds.parVec = '[P.w*1e6,P.a*1e6,P.hx/P.a,P.hy/P.w,P.oblong/10,P.maxdef]';

%% define optimization parameter space
sampbnds = [0.1,0.9]; % limits in which starting point is to be found
                      % e.g. [0.1,0.9] samples 10-90% of bounds
xnom = eval(optimds.parVec);    %nominal params
xmin = [0.90, 0.50, 0.35, 0.60, 0.10, 0.05];
xmax = [1.05, 0.60, 0.60, 0.85, 0.35, 0.30];

optimds.sampbnds = sampbnds;
optimds.xmin = xmin;
optimds.xmax = xmax;
optimds.randomStart = randomStart;

optimds.datLoc = datLoc;
optimds.F0res = {};
dtID = datestr(now,30);
save([datLoc,'optimds',dtID,'.mat'],'optimds');
% save([datLoc,'optimds',datestr(now,30),'.mat'],'optimds');
%% Execute search for optimized design
global itrPath;     %path of txt file to store iteration data
global y;

% Initialize pseudorandom number generator state
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

finished = 0;
strPoint = strPoint0;

while ~finished
    y = zeros(size(xmin));
    
    %% Pick starting point
    if randomStart 
        % Pick a random point within the parameter space, check it has 
        % modes around our target before proceeding with the optimization
        goodStart = 0;
        display('======================');
        display('Finding starting point');
        tStartSP = tic;
        while ~goodStart
            close all;
            
            %sample between min and max sampbnds
            xinit = ((sampbnds(2)-sampbnds(1))*rand(size(xmax))+sampbnds(1)).*(xmax-xmin)+xmin; 
            display(optimds.parVec)
            display(xinit);
            
            % convert param to geom values
            P.w = xinit(1)*1e-6;
            P.a = xinit(2)*1e-6;
            P.hx = P.a*xinit(3);
            P.hy = P.w*xinit(4);
            P.oblong = xinit(5)*10;
            P.maxdef = xinit(6);
            P.airrad = 2*P.lambda+P.w/2;
            P.th = P.w/(2*tan(P.theta*pi/180));   % beam thickness
            
            %% Solve for modes using random starting point
            try 
                addpath('C:\Users\CCHIA\Documents\GitHub\OMC\UnitCell1DFEM\');
                PS = P;
                PS.strPoint = strPoint+1;
                PS.nItr = 0;
                PS.solveMech = 1;                        % 1 to solve for mechanics
                PS.solveOpt = 1;                         % 1 to solve for optics
                PS.calcG = 1*(PS.solveMech && PS.solveOpt);% 1 to calculate optomechanical coupling
                PS.calcS = 0*PS.solveMech;                % 1 to calculate strain coupling
                
                % plotting & saving
                PS.plotgeom = 0;                         % 1 to plot the geometry
                PS.storeMPH = 0;                         % 1 to save COMSOL model file
                PS.plotMech = 0*PS.solveMech;             % 1 to plot displacement and strain profiles
                PS.plotOpt = 0*PS.solveOpt;               % 1 to plot E-field profiles
                PS.plotStrCpl = 0*PS.calcS;               % 1 to plot strain coupling profile
                
                PS.solveasym = 0;                        % 1 to solve for antisymmetric bands
                PS.completeBandGaps = 0;                 % 1 to plot complete bandgaps (across all symmetries)
                PS.plotgeom = 0;                         % 1 to plot the geometry
                PS.savedat = 0;                          % 1 to save data structures
                PS.savebndplot = 0;                      % 1 to save bandstructure plot
                PS.saveplots = 0;                        % 1 to save displacement and strain profiles
                PS.mbeveny = 1;                          % 1 to find even mechanical mode about y
                PS.mbevenz = 0;                          % 1 to find even mechanical mode about z
                
                % filter out designs with hx > hy --> gives weird swelling
                % modes around 6-7 GHz - not of interest for now
                % (added 1/5/2020)
                if PS.hx > PS.hy
                    error('hx > hy - skip this design')
                end
                
                % test geometry creation
                [~] = CreateNanobeamGeom(PS);
                
                % solve for symmetric bands to get midgap
                % use first bandgap above 4GHz since that's usually where
                % flapping mode localizes
                bds0 = solveBands(PS,datLoc);
                [~,gapIdx] = min(abs(bds0.sym.midGap-4e9));
                gapFlap = gapIdx(1);
%                 gap1 = bds0.sym.gapSize(1);
%                 gap2 = bds0.sym.gapSize(2);
                if bds0.sym.gapSize(gapFlap) > 1e9
                    
                    PS.freq = bds0.sym.midGap(gapFlap);
                    [spds,~] = RunNanobeamFEM(PS,datLoc);
                    if spds.cpl.gMax > PS.g0min
                        goodStart = 1;
                    end
                else
                    error(['no band gaps > 1 GHz'])
                end
            catch lasterror
                display(lasterror.message)
                display('Error stack:')
                for lsi=1:length(lasterror.stack)
                    display([num2str(lsi),': line ',...
                             num2str(lasterror.stack(lsi).line),...
                             ' in file: ',lasterror.stack(lsi).file])
                end
                F0arr = {PS.th,PS.a,PS.nholes,PS.ndef,PS.w,PS.hx,PS.hy,...
                         PS.maxdef,PS.oblong,PS.nItr,PS.strPoint,...
                         lasterror.message};
                optimds.F0res(end+1,1:length(F0arr)) = F0arr;
                save([datLoc,'optimds',dtID,'.mat'],'optimds');
                F0files = dir([datLoc,'*_itr0_*']);
                for f0i = 1:length(F0files)
                    delete([datLoc,F0files(f0i).name]);
                end
                errfiles = dir([datLoc,'error_*']);
                for ei = 1:length(errfiles)
                    delete([datLoc,errfiles(ei).name]);
                end
            end
        end
        tEndSP = toc(tStartSP);
        display(['Starting point found, time elapsed/s = ',num2str(tEndSP)]);
        strPoint = strPoint + 1;
        P.strPoint = strPoint;
        P.nItr = 1;
        F0files = dir([datLoc,'*_itr0_*']);
        for f0i = 1:length(F0files)
            delete([datLoc,F0files(f0i).name]);
        end
    else 
        % Use nominal parameters as starting point
        xinit = xnom;
        strPoint = strPoint + 1;
        P.strPoint = strPoint;
        P.nItr = 1;
    end 
    
    %% Create .txt file to assemble iteration results
    
    hdrP = {'th','a','nholes','ndef','w','hx','hy', ...
            'maxdef','oblong','optWvl','Q','wM',...
            'gMB','gPE','gMax'};
    hdrF = {'F','nItr'};
    hdrS = '';
    if P.calcS
        hdrS = {'LSiV_z111','SiVxpos','SiVypos','SiVdepth',...
                'LSiV_z111AperMax','LSiV_z111AperMean'};
    end
    optimds.hdrs = [hdrP,hdrS,hdrF];
    
    save([datLoc,'optimds',dtID,'.mat'],'optimds');
    
    %% Do fmin search
%     P.nItr = 1;
    try
        if singlePass
            %% single run
            xinit = eval(optimds.parVec);    %nominal params
            pp = fitness(xinit,optimds.fitStr,datLoc);
        else
            pp = fminsearchbnd(@(x) -1*fitness(x,optimds.fitStr,datLoc),xinit,xmin,xmax);
        end
        finished = 1;
    catch lasterror
        % catches: error('minimum step size hit') in fitness
        display(lasterror.message)
        display('Error stack:')
        for lsi=1:length(lasterror.stack)
            display([num2str(lsi),': line ',...
                     num2str(lasterror.stack(lsi).line),...
                     ' in file: ',lasterror.stack(lsi).file])
        end
        if randomStart == 0
            finished = 1;   % finish run if rand start not chosen
        end
        pp.lasterror = lasterror;
        %i.e. if random starting point chosen and error, finished = 0
        % then pick new starting point and run fminsearchbnd again
    end
    
    % if max no. of starting points specified and reached, end run
    if (maxSPtoRun > 0) && (P.strPoint-strPoint0 >= maxSPtoRun)
        finished = 1;
    end
end 


end %of function pp = Optimize...


%% FITNESS FUNCTION
function F = fitness(x,fitStr,datLoc)
global P;
global file
global y;
global optimds
global itrPath;
global dtID

% addpath('L:\Individuals\cchia\OMC\FDTD\nanobeamFDTD\FDTD scripts\');
addpath('C:\Users\CCHIA\Documents\GitHub\OMC\UnitCell1DFEM\');
close all;

tStartF = tic;

P.w = x(1)*1e-6;
P.a = x(2)*1e-6;
P.hx = P.a*x(3);
P.hy = P.w*x(4);
P.oblong = x(5)*10;
P.maxdef = x(6);

% redefine terms dependent on parameter vector
P.airrad = 2*P.lambda+P.w/2;
P.zSlc = P.th/2 - 50e-9;
P.LStats.zmin = P.th/2-70e-9;   % z-coords relative to center of beam
P.LStats.zmax = P.th/2-20e-9;
P.th = P.w/(2*tan(P.theta*pi/180));   % beam thickness


% % import COMSOL class
% import com.comsol.model.*
% import com.comsol.model.util.*
% ModelUtil.showProgress(true);

% stop optimization if next iteration changes by less than 1% from previous
% iteration
% display(P.parVec)
% display(x);
parVars = strsplit(optimds.parVec(2:end-1),',');
dp = abs(x-y)./y*100;
display('============================================')
display(['Solving for SP ',num2str(P.strPoint),' itr ',num2str(P.nItr)])
display('Parameter vector and % change from last itr:')
for vi = 1:length(x)
    display([parVars{vi},' = ',num2str(x(vi),'%.4f'),' (',num2str(dp(vi),'%.2f'),'%)'])
end
if(max(abs(x-y)./y)<0.01)
    error('minimum step size hit');
end
% display(['geometry change: ',num2str(max(abs(x-y)./y)*100),'%'])
y = x;  %store param vector x from previous iteration

%% Solve for modes using input parameter vector x
try
    % filter out designs with hx > hy --> gives weird swelling
    % modes around 6-7 GHz - not of interest for now
    % (added 1/5/2020)
    if P.hx > P.hy
        error('hx > hy - skip this design')
    end
    
    % test geometry creation
    [~] = CreateNanobeamGeom(P);

    % get bands
    PB = P;
    bds = solveBands(PB,datLoc);
%     gapIdx = find(bds.sym.midGap>4e9);
    [~,gapIdx] = min(abs(bds.sym.midGap-4e9));
    gapFlap = gapIdx(1);
    if bds.sym.gapSize(gapFlap) > 1e9
                    
        P.freq = bds.sym.midGap(gapFlap);

        % solve
        [omds,~] = RunNanobeamFEM(P,datLoc);

        % calc F
        Fnew = eval(fitStr);
        omds.cpl.F = Fnew;
        [omds.cpl.FMax,omds.cpl.FidxMax] = max(Fnew);
        F = omds.cpl.FMax;
        display(['max fitness F = ',num2str(F,'%.5e')]);
        fprintf('\n')

        cpl = omds.cpl;
        mIdx = cpl.FidxMax;
        resStrCpl = [];
        if P.calcS
            resStrCpl = [cpl.SiV_111.full.LSiVMax(mIdx),cpl.SiV_111.full.LSiVMaxPos(1,mIdx)*1e9,...
                cpl.SiV_111.full.LSiVMaxPos(2,mIdx)*1e9,(P.th/2-cpl.SiV_111.full.LSiVMaxPos(3,mIdx))*1e9,...
                cpl.SiV_111.full.LSiVAperMax(mIdx),cpl.SiV_111.full.LSiVAperMean(mIdx)];
        end
        resArr = [P.th,P.a,P.nholes,P.ndef,P.w,P.hx,P.hy,...
                P.maxdef,P.oblong,...
                cpl.optWvl,cpl.Q,cpl.mechFreq,...
                cpl.gMBmax,cpl.gPEmax,cpl.gOMmax,...
                resStrCpl,...
                cpl.FMax,P.nItr];
        optimds.sp(P.strPoint).res(P.nItr,1:length(resArr)) = resArr(:);
        save([datLoc,'optimds',dtID,'.mat'],'optimds');

        % plot freq on bands
        hdl = open([datLoc,bds.sym.P.fileBase,'_fullBands.fig']);
        hold on
        plot([0,1],[1,1]*cpl.mechFreq*1e-9,'--','Color',[0 102 51]/255,'Linewidth',2)
        ylim([0 10])
        saveas(hdl,[datLoc,omds.P.fileBase,'_fullBands.fig']);
        saveas(hdl,[datLoc,omds.P.fileBase,'_fullBands.png']);
        close all
        clear bds

        display(['Data saved for start point ',num2str(P.strPoint),...
                        ' iteration ',num2str(P.nItr)])
        P.nItr = P.nItr + 1;
    else
        error(['no band gaps > 1 GHz'])
    end
    
    
catch lasterror
    display(lasterror.message)
    display('Error stack:')
    for lsi=1:length(lasterror.stack)
        display([num2str(lsi),': line ',...
                 num2str(lasterror.stack(lsi).line),...
                 ' in file: ',lasterror.stack(lsi).file])
    end
    F = 0;
    display(['fitness value of F = ',num2str(F)]);
    if F == 0
        F0arr = {P.th,P.a,P.nholes,P.ndef,P.w,P.hx,P.hy,...
                 P.maxdef,P.oblong,P.nItr,P.strPoint,...
                 lasterror.message};
        optimds.F0res(end+1,1:length(F0arr)) = F0arr;
        
        save([datLoc,'optimds',dtID,'.mat'],'optimds');
        F0files = dir([datLoc,'*sp',num2str(P.strPoint),...
                        '_itr',num2str(P.nItr),'*']);
        for f0i = 1:length(F0files)
            delete([datLoc,F0files(f0i).name]);
        end
        errfiles = dir([datLoc,'error_*']);
        for ei = 1:length(errfiles)
            delete([datLoc,errfiles(ei).name]);
        end
    end
end


% if F == 0
%     optimds.F0.res(end+1,:) = [P.th,P.a,P.nholes,P.ndef,P.w,P.hx,P.hy,...
%                                P.maxdef,P.oblong,P.nItr,P.strPoint];
% 	save([datLoc,'optimds',dtID,'.mat'],'optimds');
%     F0files = dir([datLoc,'*sp',num2str(P.strPoint),...
%                     '_itr',num2str(P.nItr),'*']);
%     for f0i = 1:length(F0files)
%         delete([datLoc,F0files(f0i).name]);
%     end
% end

tEndF = toc(tStartF);
display(['Fitness evaluated, time elapsed/mins = ',num2str(tEndF/60,'%.2f')]);
end %of function F = fitness(x,datLoc)

%% Helper functions - DON'T NEED TO TOUCH THESE FUNCTIONS
function [x,fval,exitflag,output]=fminsearchbnd(fun,x0,LB,UB,options,varargin)
% FMINSEARCHBND: FMINSEARCH, but with bound constraints by transformation
% usage: x=FMINSEARCHBND(fun,x0)
% usage: x=FMINSEARCHBND(fun,x0,LB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options,p1,p2,...)
% usage: [x,fval,exitflag,output]=FMINSEARCHBND(fun,x0,...)
%
% arguments:
%  fun, x0, options - see the help for FMINSEARCH
%
%  LB - lower bound vector or array, must be the same size as x0
%
%       If no lower bounds exist for one of the variables, then
%       supply -inf for that variable.
%
%       If no lower bounds at all, then LB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
%  UB - upper bound vector or array, must be the same size as x0
%
%       If no upper bounds exist for one of the variables, then
%       supply +inf for that variable.
%
%       If no upper bounds at all, then UB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
% Notes:
%
%  If options is supplied, then TolX will apply to the transformed
%  variables. All other FMINSEARCH parameters should be unaffected.
%
%  Variables which are constrained by both a lower and an upper
%  bound will use a sin transformation. Those constrained by
%  only a lower or an upper bound will use a quadratic
%  transformation, and unconstrained variables will be left alone.
%
%  Variables may be fixed by setting their respective bounds equal.
%  In this case, the problem will be reduced in size for FMINSEARCH.
%
%  The bounds are inclusive inequalities, which admit the
%  boundary values themselves, but will not permit ANY function
%  evaluations outside the bounds. These constraints are strictly
%  followed.
%
%  If your problem has an EXCLUSIVE (strict) constraint which will
%  not admit evaluation at the bound itself, then you must provide
%  a slightly offset bound. An example of this is a function which
%  contains the log of one of its parameters. If you constrain the
%  variable to have a lower bound of zero, then FMINSEARCHBND may
%  try to evaluate the function exactly at zero.
%
%
% Example usage:
% rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
%
% fminsearch(rosen,[3 3])     % unconstrained
% ans =
%    1.0000    1.0000
%
% fminsearchbnd(rosen,[3 3],[2 2],[])     % constrained
% ans =
%    2.0000    4.0000
%
% See test_main.m for other examples of use.
%
%
% See also: fminsearch, fminspleas
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 4
% Release date: 7/23/06

% size checks
xsize = size(x0);
x0 = x0(:);
n=length(x0);

if (nargin<3) || isempty(LB)
    LB = repmat(-inf,n,1);
else
    LB = LB(:);
end
if (nargin<4) || isempty(UB)
    UB = repmat(inf,n,1);
else
    UB = UB(:);
end

if (n~=length(LB)) || (n~=length(UB))
    error 'x0 is incompatible in size with either LB or UB.'
end

% set default options if necessary
if (nargin<5) || isempty(options)
    options = optimset('fminsearch');
end

% stuff into a struct to pass around
params.args = varargin;
params.LB = LB;
params.UB = UB;
params.fun = fun;
params.n = n;
params.OutputFcn = [];

% 0 --> unconstrained variable
% 1 --> lower bound only
% 2 --> upper bound only
% 3 --> dual finite bounds
% 4 --> fixed variable
params.BoundClass = zeros(n,1);
for i=1:n
    k = isfinite(LB(i)) + 2*isfinite(UB(i));
    params.BoundClass(i) = k;
    if (k==3) && (LB(i)==UB(i))
        params.BoundClass(i) = 4;
    end
end

% transform starting values into their unconstrained
% surrogates. Check for infeasible starting guesses.
x0u = x0;
k=1;
for i = 1:n
    switch params.BoundClass(i)
        case 1
            % lower bound only
            if x0(i)<=LB(i)
                % infeasible starting value. Use bound.
                x0u(k) = 0;
            else
                x0u(k) = sqrt(x0(i) - LB(i));
            end
            
            % increment k
            k=k+1;
        case 2
            % upper bound only
            if x0(i)>=UB(i)
                % infeasible starting value. use bound.
                x0u(k) = 0;
            else
                x0u(k) = sqrt(UB(i) - x0(i));
            end
            
            % increment k
            k=k+1;
        case 3
            % lower and upper bounds
            if x0(i)<=LB(i)
                % infeasible starting value
                x0u(k) = -pi/2;
            elseif x0(i)>=UB(i)
                % infeasible starting value
                x0u(k) = pi/2;
            else
                x0u(k) = 2*(x0(i) - LB(i))/(UB(i)-LB(i)) - 1;
                % shift by 2*pi to avoid problems at zero in fminsearch
                % otherwise, the initial simplex is vanishingly small
                x0u(k) = 2*pi+asin(max(-1,min(1,x0u(k))));
            end
            
            % increment k
            k=k+1;
        case 0
            % unconstrained variable. x0u(i) is set.
            x0u(k) = x0(i);
            
            % increment k
            k=k+1;
        case 4
            % fixed variable. drop it before fminsearch sees it.
            % k is not incremented for this variable.
    end
    
end
% if any of the unknowns were fixed, then we need to shorten
% x0u now.
if k<=n
    x0u(k:n) = [];
end

% were all the variables fixed?
if isempty(x0u)
    % All variables were fixed. quit immediately, setting the
    % appropriate parameters, then return.
    
    % undo the variable transformations into the original space
    x = xtransform(x0u,params);
    
    % final reshape
    x = reshape(x,xsize);
    
    % stuff fval with the final value
    fval = feval(params.fun,x,params.args{:});
    
    % fminsearchbnd was not called
    exitflag = 0;
    
    output.iterations = 0;
    output.funcount = 1;
    output.algorithm = 'fminsearch';
    output.message = 'All variables were held fixed by the applied bounds';
    
    % return with no call at all to fminsearch
    return
end

% Check for an outputfcn. If there is any, then substitute my
% own wrapper function.
if ~isempty(options.OutputFcn)
    params.OutputFcn = options.OutputFcn;
    options.OutputFcn = @outfun_wrapper;
end

% now we can call fminsearch, but with our own
% intra-objective function.
[xu,fval,exitflag,output] = fminsearch(@intrafun,x0u,options,params);

% undo the variable transformations into the original space
x = xtransform(xu,params);

% final reshape
x = reshape(x,xsize);

% Use a nested function as the OutputFcn wrapper
    function stop = outfun_wrapper(x,varargin)
        % we need to transform x first
        xtrans = xtransform(x,params);
        
        % then call the user supplied OutputFcn
        stop = params.OutputFcn(xtrans,varargin{1:(end-1)});
        
    end

end % mainline end

% ======================================
% ========= begin subfunctions =========
% ======================================
function fval = intrafun(x,params)
% transform variables, then call original function

% transform
xtrans = xtransform(x,params);

% and call fun
fval = feval(params.fun,xtrans,params.args{:});

end % sub function intrafun end

% ======================================
function xtrans = xtransform(x,params)
% converts unconstrained variables into their original domains

xtrans = zeros(1,params.n);
% k allows some variables to be fixed, thus dropped from the
% optimization.
k=1;
for i = 1:params.n
    switch params.BoundClass(i)
        case 1
            % lower bound only
            xtrans(i) = params.LB(i) + x(k).^2;
            
            k=k+1;
        case 2
            % upper bound only
            xtrans(i) = params.UB(i) - x(k).^2;
            
            k=k+1;
        case 3
            % lower and upper bounds
            xtrans(i) = (sin(x(k))+1)/2;
            xtrans(i) = xtrans(i)*(params.UB(i) - params.LB(i)) + params.LB(i);
            % just in case of any floating point problems
            xtrans(i) = max(params.LB(i),min(params.UB(i),xtrans(i)));
            
            k=k+1;
        case 4
            % fixed variable, bounds are equal, set it at either bound
            xtrans(i) = params.LB(i);
        case 0
            % unconstrained variable.
            xtrans(i) = x(k);
            
            k=k+1;
    end
end

end % sub function xtransform end