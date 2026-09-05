clear all; clc; close all;
clear P;

%% unit cell params 
P.xsect = 'rect'; 
P.beamMat = 'LN';                  % beam material name
P.celltype = 'rib';                   % specify the cell type
P.unitcell = 'rectrangular';                  % specify the shape of the unit cell
P.a = 550e-9;              % lattice constant 
P.s = 300e-9;              % unit cell spine width 
P.w = 1600e-9;              % the beam width 
P.th = 260e-9;             % height (along x) of cross (for celltype = 'hollow')
                            % or of inner block (for celltype = 'solid')
P.t = 400e-9;       % the rib width 
P.d = 15*pi/180; 

P.nperiod = 1;  % no. of periods to simulate for
P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
P.mbevenz = 0;      % 1 to find even mechanical mode about z

P.kpts = 10;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 15;                           % no. of bands to solve for

P.solveasym = 1;                        % 1 to solve for antisymmetric bands
P.completeBandGaps = 1;                 % 1 to plot complete bandgaps (across all symmetries)
P.plotgeom = 1;                         % 1 to plot the geometry
P.savedat = 1;                          % 1 to save data structures
P.savebndplot = 1;                      % 1 to save bandstructure plot
P.saveplots = 1;                        % 1 to save displacement and strain profiles
P.saveMPH = 0; 
P.bandStruct_2D = 0;                 % 1 to simulate 2D band structures
P.bandStructureDim = 1;     % Optical band only: 3 to simulate 3d structures 

% for the optical bandgap 
P.run_optical = 1;

% if we are going to save the raw data files 
P.saveRawData = 1;

%% mechanical simulation parameters 
% solid mechanics solver parameters
P.mbeveny = 0;                          % 1 to find even mechanical mode about y
P.mbevenz = 0;                          % 1 to find even mechanical mode about z
P.freq = 3e9;                             % Mechanical target frequency - set to 0 for bandstructure simulations
P.optical_freq = 200;           % optical target frequency - in THz 
P.meshSize = 4;                         % mesh quality for mechanical simulations
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
% buildRibUnitCell_LN(model,P);
% mphlaunch(model);
%% Create .txt file to assemble iteration result
% txt file name
currentDate = datestr(now,'mmddyyyy');
datLoc = [fullfile('.','test','LN_holeUnitCell_optimize',currentDate),filesep];
itrPath = [datLoc,...
        'optimization_SP_',currentDate,'.txt'];
% create directory to save files
if ~exist(datLoc,'dir')
    mkdir(datLoc)
end
itr = fopen(itrPath,'wt+');
txtToPrint = ['a th s w t d', ...
        'midGapOptical gapSizeOptical gapFracOptical midGapMech gapSizeMech gapFracMech \r\n'];
fprintf(itr,txtToPrint);
fclose(itr);
%% Single solve
% optical band sim 
P.run_optical = 1;
P.kpts = 10;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 15;                           % no. of bands to solve for
datLoc = [fullfile('.','test','LN_holeUnitCell_optical',currentDate),filesep];
P.datLoc = datLoc;
bds_optical = solveOpticalBands(P);
OpticalBands = bds_optical.opticalBand;
midGap_optical = OpticalBands.midGap(1);
gapSize_optical = OpticalBands.gapSize(1);
gapRat_optical = gapSize_optical./midGap_optical;
% mechanical band sim 
datLoc = [fullfile('.','test','LN_holeUnitCell',currentDate),filesep];
P.datLoc = datLoc;
P.run_optical=0;
P.kpts = 10;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 25;                           % no. of bands to solve for
bds_mechanical_struct = solveBands_noSym(P);
bds_mechanical = bds_mechanical_struct.full;
[gapSize_mechanical,I] = max(bds_mechanical.gapSize);
midGap_mechanical = bds_mechanical.midGap(I);
gapRat_mechanical = gapSize_mechanical./midGap_mechanical;
% save the data file 
itr = fopen(itrPath,'at+');
fprintf(itr,'%.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e \r\n',...
    P.a,P.th,P.s,P.w,P.t,P.d,midGap_optical,gapSize_optical,gapRat_optical,midGap_mechanical,gapSize_mechanical,gapRat_mechanical);
fclose(itr);
