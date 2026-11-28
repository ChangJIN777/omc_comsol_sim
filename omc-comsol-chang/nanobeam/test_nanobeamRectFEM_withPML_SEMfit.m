% front panel script for nanobeam FEM simulations (example for rect cross-section, beam-hole geometry)

clear all; close all; clc

%% geometry parameters
% unit cell params
P.xsect = 'isoFit';                       % beam cross sectional shape - 'tri' or 'rect' or 'isoFit'
P.celltype = 'hole';                    % specify what type of unit cells we are simulating
P.beamMat = 'diamond';                  % beam material name
P.anisoMat = 1;

P.a = 529e-9;                           % nominal lattice constant
P.w = 711.76e-9;                           % beam width
P.theta = 30;                           % etch angle in degrees (no effect for rect cross section)
P.th = 310.75e-9+P.w/(2*sqrt(3));                          % beam thickness
P.hx = 197e-9;                          % nominal hole height (along x-axis)
P.hy = 578e-9;                          % nominal hole width (along y-axis)
% fit the sem of the fabricated omc 
% hole params for symmetric cavity / right half of asymmetric cavity
P.nholes = 18;                          % # holes in 1/2 beam length
P.ndef = 6;                             % # of holes in 1/2 defect region
P.maxdef = 0.1926;                     % defect percentage
P.oblong = 2.276;                      % oblong parameter (zero if holes are not changed)
P.aspect_ratio = 1;                    % aspect ratio of the hole at the center (1 if holes are spherical)

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
P.wgmTaper.a_end = 650e-9;              % for endtype = 'custom': lattice constant of end hole
P.wgmTaper.hx_end = 343e-9;             % for endtype = 'custom': hole height of end hole
P.wgmTaper.hy_end = 300e-9;             % for endtype = 'custom': hole width of end hole
% P.wgmTaper.maxdef_end = 0.094932;     % defect percentage - for endtype = 'maxdef'
% P.wgmTaper.oblong_end = 2.9172;       % oblong parameter (zero if holes are not changed) - for endtype = 'maxdef'

% asymmetric cavity - specify param data struct P.PL with similar fields to
% P, for left half of asymmetric cavity
% 1: full beam, no x symmetry (SetupNanobeamFEM.m:21-23 sets mevenx = 0).
%    Needed for the as-built SEM beam, P.semHalf = 'both' below.
% 0: symmetric half model, e.g. with P.semHalf = 'average' below.
P.asymCav = 1;                          % 1 to enable asymmetric cavity
if P.asymCav
    P.PL = P;                           % with P.useManualGeom = 1 the SEM block replaces P.PL
    P.PL.nholes = 3+P.ndef;
    P.PL.wvgmir = 5;
end

%% manual hole geometry from SEM data - switch
% 0: holes follow the maxdef/oblong/taperFunc taper set above.
% 1: holes come from the SEM table P.semFile, read by readSEMHoles.m, and
%    CreateNanobeamGeom ignores maxdef/oblong/taperFunc/taperTo/cavlen.
%    The file also sets P.holeatctr, P.nholes, P.ndef, P.wvgmir and the
%    mirror cell P.a/P.hx/P.hy (each change from the values above is
%    printed). P.maxdef/P.oblong stay: PlotDefectCells prints them.
%    Default: the beam as built. P.semHalf = 'both' with P.asymCav = 1 puts
%    the left half in P.PL (replacing P.PL.nholes/P.PL.wvgmir set above),
%    and P.semEndTaper = 'holes' keeps the one-sided end taper as ordinary
%    holes. Every file row is then one hole at its file position and size,
%    the placeholder rows included; nothing is padded or cut. Each beam end
%    lies half the last hole spacing past the last hole
%    (CreateNanobeamGeom.m:424). 'wvgmir' builds the same holes but takes
%    that margin from the last spacing before the taper (3.4 nm longer here).
%    Symmetric half model instead: P.asymCav = 0 above, P.semHalf =
%    'average' and P.semEndTaper = 'drop'.
P.useManualGeom = 1;

if P.useManualGeom
    % SEM table [nm], header hole,a_nm,hx_nm,hy_nm; hole = position index
    % along the beam, a_nm = spacing from hole k to hole k+1 (readSEMHoles)
    scriptDir = fileparts(mfilename('fullpath'));
    if isempty(scriptDir)
        scriptDir = pwd;                % e.g. run as a section
    end
    P.semFile = fullfile(scriptDir,'SEM_data','omc_v2_dose460_a0MN3w1.csv');
    P.semCenterHole = [];               % file hole index of the cavity centre; [] = from hx/hy symmetry
    P.semHalf = 'both';                 % 'both' (asymCav = 1) | 'average' | 'right' | 'left'
    P.semEndTaper = 'holes';            % one-sided end taper: 'holes' | 'wvgmir' | 'drop'
    P.semNHoles = [];                   % holes per half incl. centre; [] = longest half in the file
    P.semMirror = [];                   % [a hx hy] [m] for padding; [] = the file's placeholder rows

    if strcmp(P.semHalf,'both') ~= logical(P.asymCav)
        error('useManualGeom: P.semHalf = ''both'' goes with P.asymCav = 1, and only with it');
    end
    S = readSEMHoles(P.semFile,'CenterHole',P.semCenterHole,'Half',P.semHalf,...
        'EndTaper',P.semEndTaper,'NHoles',P.semNHoles,'Mirror',P.semMirror);
    if strcmp(P.semHalf,'both')         % as built: one hole per file row
        nBuilt = S.nholes + S.wvgmir + S.L.nholes + S.L.wvgmir - S.holeatctr;
        if nBuilt ~= height(S.file)
            warning('useManualGeom:asBuilt',['%d holes built from %d file rows: ',...
                'P.semNHoles or P.semEndTaper = ''drop'' padded or cut the beam'],...
                nBuilt,height(S.file));
        end
    end

    % scalar fields the file sets; print where it differs from the values above
    semSet = {'holeatctr',S.holeatctr; 'nholes',S.nholes; 'ndef',S.ndef; 'wvgmir',S.wvgmir};
    if ~isempty(S.mirror)
        semSet = [semSet; {'a',S.mirror(1); 'hx',S.mirror(2); 'hy',S.mirror(3)}];
    end
    for kSet = 1:size(semSet,1)
        fSet = semSet{kSet,1};  vSet = semSet{kSet,2};
        if abs(P.(fSet) - vSet) > 1e-15
            if vSet < 1e-3, uSc = 1e9; uStr = ' nm'; else, uSc = 1; uStr = ''; end
            fprintf('SEM file sets P.%s = %g%s (was %g%s)\n',fSet,vSet*uSc,uStr,...
                P.(fSet)*uSc,uStr);
        end
        P.(fSet) = vSet;
    end
    P.a_hole  = S.a_hole;   P.hx_hole = S.hx_hole;   P.hy_hole = S.hy_hole;
    P.useManualWvgMir = S.wvgmir > 0;
    P.a_wvgmir = S.a_wvgmir;  P.hx_wvgmir = S.hx_wvgmir;  P.hy_wvgmir = S.hy_wvgmir;
    if P.asymCav                        % left half for CreateNanobeamGeom_asym
        P.PL = rmfield(P,'PL');
        P.PL.nholes  = S.L.nholes;    P.PL.ndef = S.L.ndef;    P.PL.wvgmir = S.L.wvgmir;
        P.PL.a_hole  = S.L.a_hole;    P.PL.hx_hole = S.L.hx_hole;    P.PL.hy_hole = S.L.hy_hole;
        P.PL.useManualWvgMir = S.L.wvgmir > 0;
        P.PL.a_wvgmir = S.L.a_wvgmir;  P.PL.hx_wvgmir = S.L.hx_wvgmir;  P.PL.hy_wvgmir = S.L.hy_wvgmir;
    end

    % --- fabrication check of every hole, end taper included, per half.
    % CreateNanobeamGeom.m:72-84 only tests the nominal P.hx/P.hy, and its
    % manual branches skip the gap rules. Overlap is an error; lithography
    % limits (same numbers as CreateNanobeamGeom) only warn, since SEM data
    % is as-built.
    semHalves = {P};
    if P.asymCav
        semHalves{2} = P.PL;
    end
    for kh = 1:numel(semHalves)
        Q = semHalves{kh};
        aAll  = [Q.a_hole;  Q.a_wvgmir]';
        hxAll = [Q.hx_hole; Q.hx_wvgmir]';
        hyAll = [Q.hy_hole; Q.hy_wvgmir]';
        bridge = aAll(2:end) - (hxAll(2:end)+hxAll(1:end-1))/2;    % hole k-1 to k
        if ~Q.holeatctr
            bridge = [aAll(1)-hxAll(1), bridge]; %#ok<AGROW> central pair
        end
        wall = Q.w - hyAll;                                         % both sides
        if any(bridge <= 0) || any(wall <= 0)
            error(['useManualGeom: half %d: holes overlap or cut the beam (min ',...
                'bridge %.1f nm, min w-hy %.1f nm)'],kh,min(bridge)*1e9,min(wall)*1e9);
        end
        if min(bridge) < 50e-9 || min(wall) < 100e-9 || min([hxAll,hyAll]) < 50e-9
            warning('useManualGeom:fab',['half %d outside lithography limits: min ',...
                'bridge %.1f nm (<50), min w-hy %.1f nm (<100) or min hole size %.1f nm (<50)'],...
                kh,min(bridge)*1e9,min(wall)*1e9,min([hxAll,hyAll])*1e9);
        end
    end
end

P.lambda = 1550e-9;                     % target optical wavelength
P.nbeam = 2.386;                        % the refractive index of diamond at telecom

% Disorder
P.stdDev = [0,0];                       % standard deviation of hole dimensions (hh,hw)
P.stdDevPos = 0;                        % standard deviation of hole positions
P.asym = 0;                             % cross-section asymmetry (target y-offset in bottom apex position)
P.extractLocMechModes = 1;              % if we are filtering mode based on how confined they are

%% specify simulation/calculation/plot/save options
P.solveMech = 1;                        % 1 to solve for mechanics
P.solveOpt = 1;                         % 1 to solve for optics
P.calcG = 1*(P.solveMech && P.solveOpt);% 1 to calculate optomechanical coupling
P.calcS = 0*P.solveMech;                % 1 to calculate strain coupling
P.solveMechPML = 1;                     % 1 to solve for mechanical Q (future implementation)

% plotting & saving
P.plotgeom = 0;                         % 1 to plot the geometry

P.storeMPH = 0;                         % 1 to save COMSOL model file
P.plotMech = 1*P.solveMech;             % 1 to plot displacement and strain profiles
P.plotOpt = 1*P.solveOpt;               % 1 to plot E-field profiles
P.plotStrCpl = 1*P.calcS;               % 1 to plot strain coupling profile

%% mechanical simulation parameters 
% solid mechanics solver parameters
P.mevenx = 1;                           % +/-1 to find even/odd mode about x; 0 for fixed BC (asymCav = 1 forces 0)
P.meveny = 1;                           % +/-1 to find even/odd mode about y
P.mevenz = 1;                           % +/-1 to find even/odd mode about z
P.freq = 8e9;                           % target mechanical frequency
P.mneigs = 20;                          % # of eignevalues to find
P.mMesh = 1;                            % mesh quality for mechanical simulations
P.mAdjMesh = 1;                         % adjust mesh if DOFs exceed max_dof

% rotate crystal orientation of elasticity matrix
% ccw rotation in deg from <100> inplane direction about <100> surface normal
P.rxtal = 0;
P.rxtalInFilename = 1;

%% optical simulation parameters
% rf module solver parameters
P.oevenx = (-1)^(P.holeatctr);            % +/-1 to find even/odd optical mode about x (-1 == fundamental for hole in center)
P.oeveny = -1;                          % +/-1 to find even/odd optical mode about y (-1 == TE-like)
P.oevenz = 1;                           % +/-1 to find even/odd optical mode about z 
P.oneigs = 1;                           % # of eigenvalues to find
P.oMesh = 3;                            % mesh quality for optical simulations
P.oAdjMesh = 1;                         % adjust mesh if DOFs exceed max_dof
P.airrad = 2*P.lambda+P.w/2;            % radius of air cylinder surrounding nanobeam

%% OM coupling parameters
P.g0min = 80e3;                         % min g0 above which to save plots for

%% SiV strain coupling: susceptibilities and positions
% specify SiV axis - [1 1 1],[-1 1 1],[-1 -1 1],[1 -1 1]
% for multiple axes, specify as arrays in cell
P.zSiV = {[1 1 1]};%,[-1 1 1],[-1 -1 1],[1 -1 1]};

% specify custom coordinates to plot at
P.xSlc = 0; % relative to center
P.ySlc = 0;
P.zSlc = 0;

% define aperture to run stats (min, max, mean, stddev) on strain coupling
P.LStats.xmin = []; % leave empty to define based on generated geometry
P.LStats.xmax = [];
P.LStats.ymin = 0;
P.LStats.ymax = 60e-9;
P.LStats.zmin = P.th/2-80e-9;   % z-coords relative to center of beam
P.LStats.zmax = P.th/2;

%% Mechanical PML simulation settings (future implementation)
% P.PMLmesh = 5;
P.PMLmeshDiv = 20;
P.PMLLen = 10e-6;
P.PMLstr = 0.008;

%% Simulation settings/in
P.max_dof = 5e6;                        % max # of degrees of freedom

%% single run
% data location to save files in
currentDate = datestr(now,'mmddyyyy');
datLoc = [fullfile('.','test','1D_OMC_hole_isoFit',currentDate),filesep];
[ds,model] = RunNanobeamFEM(P,datLoc);


% %% debugging the cavity geometry (no solve)
% import com.comsol.model.*
% import com.comsol.model.util.*
% ModelUtil.clear(); clear model
% model = ModelUtil.create('model');
% % hole arrays first, as RunNanobeamFEM.m:78-82 does (BuildNanobeamFEM takes (model,P))
% if P.asymCav
%     P = CreateNanobeamGeom_asym(P,P.PL);
% else
%     P = CreateNanobeamGeom(P);
% end
% PlotDefectCells(P);
% P = LoadMaterialParams(P);
% [model,P] = BuildNanobeamFEM(model,P);
% disp(P.bndSel); disp(P.domSel);
% mphlaunch(model);