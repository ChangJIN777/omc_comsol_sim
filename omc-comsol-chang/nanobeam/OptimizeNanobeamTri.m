% Function to optimize nanobeam design for OM coupling or SiV strain
% datLoc must end with \
% Supplementary functions...

% datLoc = 'L:\Individuals\cchia\OMC\COMSOL\nanobeamOMscripts\optWithNomCubicTaper\';
% if ~exist(datLoc,'dir')
%     mkdir(datLoc)
% end

function pp = OptimizeNanobeamTri(datLoc)

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
singlePass = 0;
randomStart = 1; % 1 if we want to pick a random starting point
strPoint0 = 0; % last starting point for which iteration data was saved

maxSPtoRun = 0;     % no. of starting points to run for; 0 to run indefinitely

%% Nominal nanobeam parameters
% ===geometry parameters
% unit cell params
P.xsect = 'tri';            % beam cross sectional shape
P.beamMat = 'diamond';      % beam material name

P.a = 565.7e-9;                     % nominal lattice constant
P.w = 1074.8e-9;                    % beam width
P.theta = 45;                       % etch angle in degrees
P.th = P.w/(2*tan(P.theta*pi/180)); % beam thickness
P.hx = 226.3e-9;                      % nominal hole height (along x-axis)
P.hy = 698.6e-9;                    % nominal hole width (along y-axis)

% unit cell bands
P.bandLoc = 'D:\Files\OMC-SiV\45oBands\45o_a566nm_w1075nm_hx226nm_hy699nm_rxt45o_fullBands.fig';
P.bandmidgap = 4.69e9;

% hole params for symmetric cavity / right half of asymmetric cavity
P.nholes = 17;                          % # holes in 1/2 beam length
P.ndef = 9;                             % # of holes in 1/2 defect region
P.maxdef = 0.15;                         % defect percentage
P.consthole = 0;                        % 1 if hole size is held constant
P.oblong = 2.5;           % oblong parameter (zero if holes are not changed)
% other good params: [0.1380,3.1553], [0.1514,3.391]

% cavity params
P.holeatctr = 1;            % 1/0 for hole/dielectric in middle
P.taperFunc = 'cubic';      % linear/cubic/quadratic taper function to center hole in cavity
% P.taperTo = 'custom';     % taper to custom hole in center of cavity; disable for taper to maxdef
% P.a_ctr = 392e-9;  %for taperTo = 'custom'
% P.hx_ctr = 189e-9; %for taperTo = 'custom'
% P.hy_ctr = 177e-9; %for taperTo = 'custom'
% P.cavlen = -10e-9;         % custom cavity length between two center holes

% end waveguide mirror taper params
P.wvgmir = 0;
P.wgmTaper.func = 'cubic'; %linear, quadratic, cubic
P.wgmTaper.endtype = 'custom'; %custom, maxdef, zero (or '')
P.wgmTaper.a_end = 426e-9;  %for endtype = 'custom'
P.wgmTaper.hx_end = 213e-9; %for endtype = 'custom'
P.wgmTaper.hy_end = 173e-9; %for endtype = 'custom'
% P.wgmTaper.maxdef_end = 0.094932; % defect percentage - for endtype = 'maxdef'
% P.wgmTaper.oblong_end = 2.9172; % oblong parameter (zero if holes are not changed) - for endtype = 'maxdef'

P.nphonmir = 0;

% central taper support params
% P.tapSup.on = 1;
% P.tapSup.theta = P.theta;
% P.tapSup.wo1 = P.w;   % waveguide width at input
% P.tapSup.wo2 = P.w;   % waveguide width at output
% P.tapSup.ws = 1.3*P.w;    % max width of waveguide support
% P.tapSup.tprlgth = 15e-6; % length (of 1/2-beam) to taper over
% P.tapSup.strghtlgth = 0;  % length of straight sections at both ends of support
% P.tapSup.extrlgth1 = 0;   % length of straight section at input end (in addition to strghtlgth)
% P.tapSup.extrlgth2 = 0;   % length of straight section at output end (in addition to strghtlgth)
% P.tapSup.centLength = 0;  % length of straight section in center of taper
% P.tapSup.nRes = 20;       % resolution for taper

% asymmetric cavity: set with same parameters above
P.asymCav = 0;
if P.asymCav
    P.PL = P;
    P.PL.nholes = 4+P.ndef;
    P.PL.wvgmir=7;
end

P.lambda = 1550e-9;                      % target optical wavelength

%Disorder
P.stdDev = [0,0]; %0.05*[P.hh,P.hw];             % standard deviation of hole dimensions (hh,hw)
P.stdDevPos = 0; %0.05*P.a;                  % standard deviation of hole positions
P.asym = 0; %0.02*P.w;                       % cross-section asymmetry (target y-offset in bottom apex position)

% specify simulations/calculations to run
P.solveMech = 1;        % solve for optics
P.solveOpt = 1;         % solve for mechanics
P.solveTwoCompts = 0;   % to solve both opt and mech in separate components in COMSOL
P.calcG = 1*(P.solveMech && P.solveOpt);    % calculate optomechanical coupling
P.calcS = 1*P.solveMech;    % calculate strain coupling

% plotting & saving
P.plotgeom = 1;                     % 1 to plot the geometry
P.storeMPH = 0;                     % 1 to save COMSOL MPH model
P.plotMech = 1*P.solveMech;         % 1 to plot displacement and strain profiles
P.plotOpt = 1*P.solveOpt;           % 1 to plot E-field profiles
P.plotStrCpl = 1*P.calcS;           % 1 to plot strain coupling profile

% ===mechanical simulation parameters 
% solid mechanics solver parameters
P.mevenx = 1;       % +/-1 to find even/odd mode about x; 0 for fixed BC
P.meveny = 1;       % +/-1 to find even/odd mode about y; 0 for fixed BC
P.mevenz = 1;       % +/-1 to find even/odd mode about z; 0 for fixed BC
P.freq = 8e9;       % target mechanical frequency
% P.freqList = [4.487,7.156,9.055]*1e9;
P.freqList = [4.69,7.711,9.556]*1e9;
P.mneigs = 10;      % # of eignevalues to find
P.mMesh = 3;        % mesh quality
P.mAdjMesh = 1;     % adjust mesh if DOFs exceed max_dof

% material properties
P.E = 1050e9;       % Young's modulus
P.rho = 3500;       % density
P.nu = 0.2;         % Poisson ratio
 
% Anisotropic elasticity matrix - COMSOL v4+ ordering: 
% 11, 12,22, 13,23,33, 14,24,34,44, 15,25,35,45,55, 16,26,36,46,56,66
c11 = 1076e9;
c12 = 125e9;
c44 = 578e9;
P.D = [c11, c12,c11, c12,c12,c11, 0,0,0,c44, 0,0,0,0,c44, 0,0,0,0,0,c44];

P.rxtal = 45;           % ccw rotation in deg from <100> inplane direction about <100> surface normal

P.useStress = 0;        % 1 to include instrinsic stress in calculations
P.stress = 1e9;         % Intrinsic stress (Pa) 

% ===optical simulation parameters
% rf module solver parameters
P.oevenx = -1;      % 1 to find even optical mode about x (-1 == fundamental for hole in center)
P.oeveny = -1;      % 1 to find even optical mode about y (-1 == TE-like)
P.oevenz = 1;       % 1 to find even optical mode about z 
P.oneigs = 1;       % # of eigenvalues to find
P.oMesh = 4;        % mesh quality
P.oAdjMesh = 1;     % adjust mesh if DOFs exceed max_dof

% material properties
P.nbeam = 2.386;                        % index of refraction in material
P.airrad = 2*P.lambda+P.w/2;            % radius of surrounding air cylinder

% ===OM coupling: photoelastic tensor components for g calculation
P.p11 = -0.25;
P.p12 = 0.043;
P.p44 = -0.172;
P.g0min = 80e3; %min g0 above which to save plots for

% ===SiV strain coupling: susceptibilities and positions
P.posSiV = 1;                      % 1 to solve for SiV strain coupling
% parameters for lambda calculation
% P.xSiV depends on maxdef - define within loops
% P.zSiV depends on th - define within loops

P.dSiVmin = 20e-9;           % min depth of SiV
P.xSiV = P.a*(1-P.maxdef)/2; 
P.zSiV = linspace(25e-9,P.th-P.dSiVmin,50);
P.ySiV = 0;

P.dg = 1.3e15;    % strain susceptibilities
P.fg = -250e12;

% ===FDTD settings
P.storeFDTD = 0;                        % 1 to keep solved FDTD file
P.finemeshFDTD = 1;
P.xsymFDTD = -1;
P.ysymFDTD = -1;
P.zsymFDTD = 0;
P.inittime = 40e-15;
P.pulselen = 150e-15;

% ===Simulation settings
P.max_dof = 5e6;                      % max # of degrees of freedom

%% define fitness function
% P.fitStr = ['min([abs(omds.cpl.gMax*omds.cpl.optWvl/299792458),1.8e5/1.96e14])/(1.8e5/1.96e14)',...
%             '*min([omds.cpl.Q,1e6])/1e6',...
%             '*ofdtd.Trans*abs(omds.cpl.lambdaG_111Max)'];
P.fitStr = ['min(abs(omds.cpl.gOM.*omds.cpl.optWvl./299792458),1.0e5/1.96e14)./(1.0e5/1.96e14)',...
            '.*min([omds.cpl.Q,4e5])./4e5',...
            '.*ofdtd.Trans.*abs(omds.cpl.lambdaG_111)'];
P.parVec = '[P.oblong/10,P.maxdef]';

%% define optimization parameter space
sampbnds = [0.1,0.9]; % limits in which starting point is to be found
                      % e.g. [0.1,0.9] samples 10-90% of bounds
% xnom = [P.hx/P.a,P.hy/P.w,P.oblong/10,P.maxdef,P.w*1e6];    %nominal params
% xmax = [0.6,0.8,0.35,0.2,1]; %bounds for search
% xmin = [0.4,0.6,0.25,0.15,0.8];
% xmax = [0.45,0.65,0.3,0.25,0.96]; %bounds for search
% xmin = [0.3,0.5,0.15,0.1,0.91];
xnom = eval(P.parVec);    %nominal params
xmin = [0.15   0.1];
xmax = [0.4    0.2];

paramPath = [datLoc,'optimization_params.txt'];
ptxt = fopen(paramPath,'wt+');
paramTxt = ['Fitness func F = ',P.fitStr,'\r\n',...
            'Param vector x = ',P.parVec,'\r\n',...
            'xmin = ',mat2str(xmin),'\r\n',...
            'xmax = ',mat2str(xmax),'\r\n',...
            'sample bnds = ',mat2str(sampbnds),'\r\n',...
            'randomStart = ',num2str(randomStart),'\r\n'];
fprintf(ptxt,paramTxt);
fclose(ptxt);

%% folders
P.scriptLoc = 'L:\Individuals\cchia\OMC\Utility Scripts\';

% required FDTD file locations
% directory path for location of FDTD files
file.path = 'L:\Individuals\cchia\OMC\FDTD\Nanocavity optimization scripts\FDTD scripts\';

% template file containing all preset FDTD elements/parameters
if P.finemeshFDTD
    file.template = [file.path,'FDTDtemplate_finemesh_xy-sym.fsp']; 
else
    file.template = [file.path,'FDTDtemplate.fsp'];
end

file.runsim = [file.path,'script_runNanobeamCavity.lsf']; % Lumerical script to run simulation and collect end time
file.QVanalysis = [file.path,'script_QVanalysis__sym.lsf']; 
file.plotfields = [file.path,'script_plotfields.lsf']; % Lumerical script to plot field profiles
file.savetmon = [file.path,'script_savetmondata.lsf']; % Lumerical script to plot field profiles

if isfield(P,'tapSup') && P.tapSup.on
    file.build = [file.path,'script_buildNanobeamTapSup.lsf'];
else
    file.build = [file.path,'script_buildNanobeamCavity.lsf']; % Lumerical script to build geomtry
end

addpath('L:\Individuals\cchia\OMC\FDTD\Nanocavity optimization scripts');
addpath('L:\Individuals\cchia\OMC\Utility Scripts');
addpath(P.scriptLoc);



%% import COMSOL class
import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.showProgress(true);

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
        display('Finding starting point');
        tStartSP = tic;
        while ~goodStart
            close all;
            
            %sample between min and max sampbnds
            xinit = ((sampbnds(2)-sampbnds(1))*rand(size(xmax))+sampbnds(1)).*(xmax-xmin)+xmin; 
            display(xinit);
            
            % convert param to geom values
            P.oblong = xinit(1)*10;
            P.maxdef = xinit(2);
            P.th = P.w/(2*tan(P.theta*pi/180));   % re-clac beam thickness
            P.airrad = 2*P.lambda+P.w/2;
            
            % create new COMSOL model for starting point
            ModelUtil.clear();
            clear spmodel
            clear spds
            spmodel = ModelUtil.create('spmodel');
            
            %% Solve for modes using random starting point
            try 
                geomOK = 0;
                % Create geometry
                try 
                    if isfield(P,'asymCav') && P.asymCav && isfield(P,'PL')
                        P = CreateNanobeamGeom_asym(P,P.PL);
                    else
                        P = CreateNanobeamGeom(P);
                    end
                    PlotDefectCells(P);
                    geomOK = 1;
                catch lasterror %if geometry cannot be fabbed
                    display(lasterror.message)
                    display('Error stack:')
                    for lsi=1:length(lasterror.stack)
                        display([num2str(lsi),': line ',...
                                 num2str(lasterror.stack(lsi).line),...
                                 ' in file: ',lasterror.stack(lsi).file])
                    end
                end 
                
                if ~geomOK
                    display('Geometry cannot be fabbed')
                else
                    %% setup simulation
                    P.solveOpt = 1; P.solveMech = 1;
                    [spmodel,P] = DrawNanobeam(spmodel,P);
                    [spmodel,spds] = SetupNanobeamFEM(spmodel,P); 

                    %% Solve for optical modes
                    display('Solving optical modes...')
                    P.solveOpt = 1; P.solveMech = 0;
                    spds.P = P;
                    [spmodel,spds] = SolveNanobeamFEM(spmodel,spds);
                    
                    % solve for symmetric mechanical bands if high-Q optical modes exist
                    if isempty(spds.ofem.hiQSol)
                        display(['No optical modes near target wavelength ',num2str(P.lambda*1e9,'%.0f'),' nm'])
                    else 
                        display(['Optical mode found at ',num2str(spds.ofem.lambda*1e9,'%.3f'),' nm']);
                        %% Solve for mechanical modes at each midgap
                        gMaxSP = 0;

                        for k = 1:length(P.freqList)
                            display(['Solving mechanical modes at wM = ',num2str(P.freqList(k)/1e9,'%.3f'),' GHz...'])
                            P.freq = P.freqList(k);
                            P.solveOpt = 0; P.solveMech = 1;
                            spds.P = P;
                            [spmodel,spds] = SolveNanobeamFEM(spmodel,spds);
                            if isempty(spds.mfem.locFreqs)
                                display('No localized mechanical modes found')
                            else
                                display(['Frequencies of localized modes: (', num2str(spds.mfem.locFreqs/1e9,'%.3f '),') GHz'])
                                %% Calculate optomechanical couplings
                                mLoc = spds.mfem.locInd;

                                fprintf('\n') 
                                display(['Determining g0 for ', ...
                                    num2str(length(spds.mfem.locInd)),' localized mechanical modes'])
                                [spds,spmodel] = CalcGOM(spds,spmodel,1,spds.mfem.locInd);
                                if abs(spds.cpl.gMax) > gMaxSP
                                    gMaxSP = abs(spds.cpl.gMax);
                                end

                            end


                        end
                        if gMaxSP > 80e3
                            goodStart = 1;
                        end
                        
                    end % of isempty(spds.ofem.hiQSol)
                end % of if ~geomOK
                
                ModelUtil.clear();
                clear spmodel
            catch lasterror
                display(lasterror.message)
                display('Error stack:')
                for lsi=1:length(lasterror.stack)
                    display([num2str(lsi),': line ',...
                             num2str(lasterror.stack(lsi).line),...
                             ' in file: ',lasterror.stack(lsi).file])
                end
            end
        end
        tEndSP = toc(tStartSP);
        display(['Starting point found, time elapsed/s = ',num2str(tEndSP)]);
        strPoint = strPoint + 1;
        P.strPoint = strPoint;
        P.nItr = 1;
    else 
        % Use nominal parameters as starting point
        xinit = xnom;
        strPoint = strPoint + 1;
        P.strPoint = strPoint;
        P.nItr = 1;
    end 
    
    %% Create .txt file to assemble iteration results
    
    % txt file name
    itrPath = [datLoc,...
        'optimization_SP_',num2str(strPoint),'.txt'];
    itr = fopen(itrPath,'wt+');
    SiVtxt = '';
    
    if P.posSiV
        SiVtxt = 'lambdaG_111 lambdaG_m111 SiVdepth ';
    end
    txtToPrint = ['theta a nholes ndef w hx hy ', ...
            'maxdef oblong ', ...
            'optWvl Q wM ',...
            'gMB gPE gMax ', ...
            'Qwvg Qsc Qt ',...
            'Qtime Trans ',...
            SiVtxt,...
            'F nItr\r\n'];
    fprintf(itr,txtToPrint);
    
    fclose(itr);
    
    %% Do fmin search
%     P.nItr = 1;
    try
        if singlePass
            %% run sweep over maxdef and oblong
            dmaxAll = 0.1:0.05:0.2;%0.1:0.025:0.3;
            oblAll = 2:0.5:3;%1.5:0.5:3.5;
            for dmax = dmaxAll
                for obl = oblAll
                    simFiles = dir([datLoc,'*obl_',num2str(obl,'%.4f'),'_dmax_',num2str(dmax,'%.4f'),'*']);
                    if isempty(simFiles) %&& ...
                       %~(dmax==0.2 && obl==2) %|| ~(dmax==0.2 && obl==1.75) || ~(dmax==0.2 && obl==2.5))
                        P.oblong = obl;
                        P.maxdef = dmax;
                        xinit = eval(P.parVec);    %nominal params
                        pp = fitness(xinit,P.fitStr,datLoc);
                    else
                        display('simulation files present')
                        pp = 0;
                    end
                end
            end
            
            
        else
            pp = fminsearchbnd(@(x) -1*fitness(x,P.fitStr,datLoc),xinit,xmin,xmax);
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
% x(1) = hx (normalized to a)
% x(2) = hy (normalized to w)
% x(3) = oblong
% x(4) = maxdef
% x(5) = width (in um)
global P;
global file
global y;
global itrPath;

addpath('L:\Individuals\cchia\OMC\FDTD\Nanocavity optimization scripts\');

close all;

tStartF = tic;

% P.hx = x(1)*P.a;
% P.hy = x(2)*x(5)*1e-6;
% P.oblong = x(3)*10;
% P.maxdef = x(4);
% P.w = x(5)*1e-6;

P.oblong = x(1)*10;
P.maxdef = x(2);
P.airrad = 2*P.lambda+P.w/2;
P.th = P.w/(2*tan(P.theta*pi/180));   % beam thickness

% import COMSOL class
import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.showProgress(true);

% stop optimization if next iteration changes by less than 1% from previous
% iteration
display(x);
if(max(abs(x-y)./y)<0.01)
    error('minimum step size hit');
end;
display(['geometry change: ',num2str(max(abs(x-y)./y)*100),'%'])
y = x;  %store param vector x from previous iteration

%% Solve for modes using input parameter vector x
try
    ModelUtil.clear();
    clear ommodel

    % create COMSOL model 
    ommodel = ModelUtil.create('ommodel');
    
    omdsMax = [];

    % Create geometry
    geomOK = 0;
    try 
        if isfield(P,'asymCav') && P.asymCav && isfield(P,'PL')
            P = CreateNanobeamGeom_asym(P,P.PL);
        else
            P = CreateNanobeamGeom(P);
        end
        PlotDefectCells(P);
        geomOK = 1;
    catch lasterror %if geometry cannot be fabbed
        display(lasterror.message)
        display('Error stack:')
        for lsi=1:length(lasterror.stack)
            display([num2str(lsi),': line ',...
                     num2str(lasterror.stack(lsi).line),...
                     ' in file: ',lasterror.stack(lsi).file])
        end
    end 
    
    % if geom can be fabbed, setup simulation and solve optical modes
    if ~geomOK
        display('Geometry cannot be fabbed')
        F = 0;
    else
        %% setup simulation
        P.solveOpt = 1; P.solveMech = 1;
        [ommodel,P] = DrawNanobeam(ommodel,P);
        [ommodel,omds] = SetupNanobeamFEM(ommodel,P); 

        %% Solve for optical modes
        display('Solving optical modes...')
        P.solveOpt = 1; P.solveMech = 0;
        omds.P = P;
        [ommodel,omds] = SolveNanobeamFEM(ommodel,omds);

        % solve for symmetric mechanical bands if high-Q optical modes exist
        if isempty(omds.ofem.hiQSol)    % after solving optical modes
            display(['No optical modes near target wavelength ',num2str(P.lambda*1e9,'%.0f'),' nm'])
            F = 0;
        else
            display(['Optical mode found at ',num2str(omds.ofem.lambda*1e9,'%.3f'),' nm']);

            %% Calculate Q's in FDTD
%             fGeomFDTD = ['sp',num2str(P.strPoint),...
%                     '_itr',num2str(P.nItr),'_',...
%                     'a_',num2str(P.a*1e9,'%.0f'),'nm_', ...
%                     'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
%                     'hx_',num2str(P.hx*1e9,'%.0f'),'nm_', ...
%                     'hy_',num2str(P.hy*1e9,'%.0f'),'nm_', ...
%                     'oblong_',num2str(P.oblong,'%.4f'),'_', ...
%                     'maxdef_',num2str(P.maxdef,'%.4f'),'_FDTD'];
            PFDTD = P;
            PFDTD.lambda = omds.ofem.lambda;
            PFDTD.pulselen = 150e-15;
%             PFDTD.fname = fGeomFDTD;
%             PFDTD.nholes = PFDTD.nholes + 5;
            % if results wrong, delete files and rerun with longer pulse
            % length
            nruns = 0;
            fdtdOK = 0;
            while ~fdtdOK
                display(['FDTD run ',num2str(nruns+1)])
                ofdtd = runNanobeamCavity(PFDTD,file,datLoc);
                nruns = nruns + 1;
                if ofdtd.Qx1 < 0 || ofdtd.Qx2 < 0 || ofdtd.Qy1 < 0 || ...
                   ofdtd.Qy2 < 0 || ofdtd.Qz1 < 0 || ofdtd.Qz2 < 0 || ...
                   ofdtd.Qtime < 0 || ofdtd.Trans > 1
                    fdtdOK = 0;
                    PFDTD.pulselen = PFDTD.pulselen + 50e-15;
                    fdtdFiles = dir([datLoc,'sp',num2str(P.strPoint),'_itr',num2str(P.nItr),'_*wl_*Vmode_*Qt_*']);
                    for fi = 1:length(fdtdFiles)
                        delete([datLoc,fdtdFiles(fi).name])
                    end
                else
                    fdtdOK = 1;
                end
                if nruns == 4
                    fdtdOK = 1;
                end
            end
                
            
            %% Solve for mechanical y-symmetric bands
            copyfile(P.bandLoc,...
                    [datLoc,'fullBandsTmp.fig']);
            
            %% Solve for mechanical modes at each midgap
            % calculate all localized mechanical modes, and keep whichever 
            % has the largest optomechanical coupling
            gMaxSP = 0;
            FMaxAll = 0;
            gMaxAll = 0;

            for k = 1:length(P.freqList)
                
                
                display(['Solving mechanical modes at wM = ',num2str(P.freqList(k)/1e9,'%.3f'),' GHz...'])
                P.freq = P.freqList(k);
                P.solveOpt = 0; P.solveMech = 1;
                omds.P = P;
                [ommodel,omds] = SolveNanobeamFEM(ommodel,omds);
                P.solveOpt = 1; P.solveMech = 1;
                omds.P = P;
                if isempty(omds.mfem.locFreqs)
                    display('No localized mechanical modes found')
                    F = 0;
                else
                    display(['Frequencies of localized modes: (', num2str(omds.mfem.locFreqs/1e9,'%.3f '),') GHz'])
                    %% Calculate optomechanical couplings
                    mLoc = omds.mfem.locInd;

                    fprintf('\n') 
                    display(['Determining g0 for ', ...
                        num2str(length(omds.mfem.locInd)),' localized mechanical modes'])
                    [omds,ommodel] = CalcGOM(omds,ommodel,1,omds.mfem.locInd);
                    [omds,ommodel] = CalcLambdaSiV(omds,ommodel,omds.mfem.locInd);
                    if abs(omds.cpl.gMax) > gMaxSP
                        gMaxSP = abs(omds.cpl.gMax);
                    end
                    
                    %% calc fitness
                    Fnew = eval(fitStr);
                    omds.cpl.F = Fnew;
                    [omds.cpl.FMax,omds.cpl.FidxMax] = max(Fnew);
                    
                    %% Save data for modes with high g and lSiV
                    if abs(omds.cpl.gMax) > 80e3 ...
                       && (real(omds.cpl.lambdaG_111Max) > 0.8e6 || ...
                           real(omds.cpl.lambdaG_m111Max) > 0.8e6)
                        close all;
                        
                        fileBase = [num2str(P.theta,'%.0f'),'o_',...
                                      'a_',num2str(P.a*1e9,'%.0f'),'nm_',...
                                      'nholes_',num2str(P.nholes,'%.0f'),'_',...
                                      'ndef_',num2str(P.ndef,'%.0f'),'_',...
                                      'nwgm_',num2str(P.wvgmir,'%.0f'),'_',...
                                      'w_',num2str(P.w*1e9,'%.0f'),'nm_',...
                                      'hx_',num2str(P.hx*1e9,'%.0f'),'nm_',...
                                      'hy_',num2str(P.hy*1e9,'%.0f'),'nm_',...
                                      'obl_',num2str(P.oblong,'%.4f'),'_',...
                                      'dmax_',num2str(P.maxdef,'%.4f')];
                        if isfield(P,'prefname')
                            fileBase = [P.prefname,'_',fileBase];
                        end
                        if isfield(P,'strPoint') && isfield(P,'nItr')
                            fileBase = ['sp',num2str(P.strPoint),'_',...
                                          'itr',num2str(P.nItr),'_',fileBase];
                        end
                        
                        
                        omds.P.fileBase = fileBase;
                        
                        %plots
                        PlotEy(ommodel,omds,omds.cpl.oSol,datLoc,1);
                        
                        % find modes with high g and lSiV
                        highGIdx = find(abs(real(omds.cpl.gOM))>80e3);
                        highLIdx = find(abs(real(omds.cpl.lambdaG_111))>0.8e6);
                        highLmIdx = find(abs(real(omds.cpl.lambdaG_m111))>0.8e6);
                        highLAllIdx = union(highLIdx,highLmIdx);
                        highGLIdx = intersect(highGIdx,highLAllIdx);
                        highGLIdx = intersect(highGLIdx,omds.mfem.locInd);
                        
                        lGMdir = {omds.cpl.lambdaGMaxDir};
                        
                        PlotDispStr(ommodel,omds,highGLIdx,datLoc,1);
                        if P.posSiV
                            PlotStrCplSiV(ommodel,omds,highGLIdx,lGMdir,datLoc)
                        end
                        
                        %bands
                        hdl = open([datLoc,'fullBandsTmp.fig']);
                        for idx = highGLIdx
                            hold on
                            plot([0,1],[1,1]*omds.mfem.freqs(idx)*1e-9,'--','Color',[0 102 51]/255,'Linewidth',2)
                        end
                        saveas(hdl,[datLoc,'fullBandsTmp.fig']);
                        close
                        
                        %data

                        
                        
                        
                        fileName = fileBase;
                        if P.calcG
                            CfBase = ['wl_',num2str(omds.cpl.optWvl*1e9,'%.0f'),'nm_',...
                                      'wM_',num2str(omds.cpl.mechFreq*1e-9,'%.2f'),'GHz_',...
                                      'g0_',num2str(omds.cpl.gMax*1e-3,'%.0f'),'kHz'];
                            fileName = [fileName,'_',CfBase];
                        end

                        if P.calcS
                            CfBase = ['lSiV_',num2str(omds.cpl.lambdaGMax*1e-6,'%.2f'),'MHz'];
                            fileName = [fileName,'_',CfBase];
                        end
                        save([datLoc,fileName,'.mat'],'omds');
                        if P.storeMPH
                            mphsave(ommodel,[datLoc,fileName,'.mph']);
                        end
                        
                        %% find mode with highest fitness for particular geometry
                        if abs(omds.cpl.FMax) > FMaxAll
    %                     if abs(omds.cpl.gMax) > gMaxAll
    %                         gMaxAll = abs(omds.cpl.gMax);
                            FMaxAll = abs(omds.cpl.FMax);
                            omdsMax = omds;

                            %save as mph, then load again
                            pathMph = [datLoc,'ommodelMax.mph'];
                            mphsave(ommodel,pathMph);
                            display('omds and ommodel stored for saving')
                            %end
                        end
                    end
                    
                    

                end


            end
                
                
                
            %% Save data and plots
            % Once we have found our optimal mechanical mode for this
            % structure, output the results, plot and save the data,
            % and calculate the fitness function
            if ~isempty(omdsMax)
                fileBase = omdsMax.P.fileBase;
                
                % save bandstructure plot
                open([datLoc,'fullBandsTmp.fig']);
                saveas(gcf,[datLoc,fileBase,'_fullBands.png']);
                delete([datLoc,'fullBandsTmp.fig']);
                
                fdtdMat = dir([datLoc,'*_FDTD_*.mat']);

                cpl = omdsMax.cpl;
                mIdx = cpl.FidxMax;
                %ofdtd = omdsMax.ofdtd;
                pathMph = [datLoc,'ommodelMax.mph'];
                ommodelMax = mphload(pathMph);
                fprintf('\n') 
                display(['Optical wvl: ',num2str(cpl.optWvl*1e9),' nm,',...
                    ' Q : ',num2str(cpl.Q), ...
                    ', Mechanical freq: ',num2str(omdsMax.mfem.freqs(mIdx)*1e-9),' GHz']);
                display(['FDTD: transmission: ',num2str(ofdtd.Trans*100,'%.1f'),'%'])
                display(['g0 = ',num2str(abs(real(cpl.gMB(mIdx)*1e-3))),...
                         ' + ',num2str(abs(real(cpl.gPE(mIdx)*1e-3))),...
                         ' = ',num2str(abs(real(cpl.gOM(mIdx)*1e-3))),' kHz']);
                if P.posSiV
                    display(['lambda, [111]: ',num2str(real(cpl.lambdaG_111(mIdx))/1e6),' MHz, ',...
                        '[-111]: ',num2str(real(cpl.lambdaG_m111(mIdx))/1e6),' MHz'])
                    display(['optimal SiV xy-position = (',...
                             num2str(cpl.xSiV(mIdx)*1e9),', ',...
                             num2str(cpl.ySiV(mIdx)*1e9),') nm; depth = ',...
                             num2str((P.th-cpl.zSiV(mIdx))*1e9),' nm']);
                end
                
                fprintf('\n') 

                % Fitness function
                % F = abs(cpl.lambdaMax)/cpl.mechFreq;
%                 F = cpl.F;
                F = cpl.FMax;
                display(['max fitness F = ',num2str(cpl.FMax,'%.5e')]);
                display(['FMaxAll = ',num2str(FMaxAll,'%.5e')]);
                fprintf('\n')

                % Append iteration fitness parameter in .txt file
                itr = fopen(itrPath,'at+');
                
                if P.posSiV
                    fprintf(itr,['%.0f %.3e %.0f %.0f %.3e %.3e %.3e',...
                        ' %.6f %.6f', ...
                        ' %.6e %.6e %.6e',...
                        ' %.6e %.6e %.6e', ...
                        ' %.6e %.6e %.6e',...
                        ' %.6e %.6e',...
                        ' %.6e %.6e',...
                        ' %.2e',...
                        ' %.6e %.0f\r\n'],...
                        P.theta,P.a,P.nholes,P.ndef,P.w,P.hx,P.hy,...
                        P.maxdef,P.oblong,...
                        cpl.optWvl,cpl.Q,cpl.mechFreq,...
                        cpl.gMBmax,cpl.gPEmax,cpl.gOMmax,...
                        ofdtd.Qwvg,ofdtd.Qsc,ofdtd.Qt,...
                        ofdtd.Qtime,ofdtd.Trans,...
                        cpl.lambdaG_111Max,cpl.lambdaG_m111Max,...
                        (P.th-cpl.zSiVMax),...
                        cpl.FMax,P.nItr);
                else
                    fprintf(itr,['%.0f %.3e %.0f %.0f %.6e %.6e %.6e',...
                        ' %.4f %.4f', ...
                        ' %.6e %.6e %.6e',...
                        ' %.6e %.6e %.6e', ...
                        ' %.6e %.6e %.6e',...
                        ' %.6e %.6e',...
                        ' %.6e %.0f\r\n'],...
                        P.theta,P.a,P.nholes,P.ndef,P.w,P.hx,P.hy,...
                        P.maxdef,P.oblong,...
                        cpl.optWvl,cpl.Q,cpl.mechFreq,...
                        cpl.gBndMax,cpl.gStrMax,cpl.gMax,...
                        ofdtd.Qwvg,ofdtd.Qsc,ofdtd.Qt,...
                        ofdtd.Qtime,ofdtd.Trans,...
                        cpl.FMax,P.nItr);
                end
                fclose(itr);  
                
                close all
                % Save plots of the optical/mechanical modes and the coupling data
                %PlotAndSave(ommodelMax,omdsMax,datLoc);
                %mphfile = dir([datLoc,'*ommodelMax.mph']);
                %display([datLoc,fdtdMat.name])

                % housekeeping - delete redundant files
                delete(pathMph);
                
%                 fdtdGeomFile = dir([datLoc,'sp',num2str(P.strPoint),...
%                 '_itr',num2str(P.nItr),'*_FDTD_geom.png']);
%                 if ~isempty(fdtdGeomFile)
%                     delete([datLoc,fdtdGeomFile.name]);
%                 end

                display(['Data saved for start point ',num2str(P.strPoint),...
                    ' iteration ',num2str(P.nItr)])
                P.nItr = P.nItr + 1;
            else
                F = 0;
                display(['fitness value of F = ',num2str(F)]);
                fprintf('\n')
            end
        end
    end
    
    ModelUtil.clear();
    clear ommodel
    clear omds
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
end

if F == 0
    F0files = dir([datLoc,'sp',num2str(P.strPoint),...
                    '_itr',num2str(P.nItr),'*']);
    for f0i = 1:length(F0files)
        delete([datLoc,F0files(f0i).name]);
    end
    
    resPath = [datLoc,'sp',num2str(P.strPoint),'_',...
               'obl_',num2str(P.oblong,'%.4f'),'_',...
               'dmax_',num2str(P.maxdef,'%.4f'),'_F_0.txt'];
    resF = fopen(resPath,'wt+');
    fprintf(resF,'F = 0 for current parameter vector');
    fclose(resF);
end

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
    function stop = outfun_wrapper(x,varargin);
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