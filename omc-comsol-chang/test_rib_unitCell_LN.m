clear all; clc; close all;
clear P;

%% rib unit cell, lithium niobate with sidewalls (the original LN rib design)
% Drives celltype 'rib' (builder buildRibUnitCell_LN.m: spine + centred rib)
% with the LN rib parameters: LiNbO3, 15 deg sidewall angle, explicit beam
% width P.w, no fillets. With these settings the merged builder reproduces
% the old LN cell's shape. (The rib is now centred in the cell instead of
% split over the cell edges: same lattice, same band structure.)
%   P.run_optical = 1  optical bands: solveOpticalBands -> runOpticalBand_1D
%   P.run_optical = 0  mechanical bands without symmetry: solveBands_noSym
% Formerly test_LN_unitCell_optical.m (optical only). The mechanical branch
% restores what the old test_LN_unitCell.m did before it was repurposed.
% Full parameter reference: rib_README.md.
P.xsect = 'rect';
P.beamMat = 'LN';                       % material: lithium niobate (default would be diamond)
P.celltype = 'rib';                     % merged rib device
P.unitcell = 'rectrangular';            % read by the optical runners; anything but
                                        % 'hexagonal' selects the square-lattice k path
P.a = 550e-9;                           % lattice constant = cell width
P.s = 300e-9;                           % spine width along y (bottom face)
P.w = 1400e-9;                          % EXPLICIT tip-to-tip beam width: wins over the
                                        % L(a) rule (L_ref/a_ref/dLda are then unused)
P.th = 260e-9;                          % slab thickness
P.t = 400e-9;                           % rib width along x (bottom face)
P.d = 15*pi/180;                        % sidewall angle in RADIANS (walls lean in by
                                        % th*tan(d) toward the top face)
P.r1 = 0;                               % no fillets: the old LN shape. With d ~= 0 they
P.r2 = 0;                               % default to 0 anyway and may not be set > 0.

P.nperiod = 1;                          % no. of periods to simulate for
P.holeatedge = 0;                       % read by the runners, unused by this builder

P.kpts = 5;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 25;                          % no. of bands to solve for

P.solveasym = 1;                        % 1 to solve for antisymmetric bands
P.completeBandGaps = 1;                 % 1 to plot complete bandgaps (across all symmetries)
P.plotgeom = 1;                         % 1 to plot the geometry
P.savedat = 1;                          % 1 to save data structures
P.savebndplot = 1;                      % 1 to save bandstructure plot
P.saveplots = 1;                        % 1 to save displacement and strain profiles
P.saveMPH = 0;
P.saveRawData = 1;                      % read by runOpticalBand_1D
P.bandStruct_2D = 0;                    % 1 to simulate 2D band structures
P.bandStructureDim = 1;                 % optical: 1 -> runOpticalBand_1D, 3 -> _3D

% choose the optical (air half-cylinder of radius 5a, y >= 0) or the
% mechanical band structure
P.run_optical = 1;
P.optical_freq = 200;                   % optical target frequency in THz (value of the
                                        % other LN rib scripts; set for your design)

%% mechanical simulation parameters
P.mbeveny = 0;                          % no symmetry planes (solveBands_noSym)
P.mbevenz = 0;                          % MUST be 0 with d ~= 0 (not mirror-symmetric in z)
P.freq = 3e9;                           % mechanical target frequency
P.meshSize = 4;                         % mesh quality
P.fixed_bc = 0;                         % no flat outer y face exists to clamp

P.anisoMat = 1;
P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg
                                        % from <100> inplane direction about <100> surface normal

%% define the maximum number of degree of freedom to limit the simulation time
P.max_dof = 3e6;                        % max # of degrees of freedom

% %% debugging the unit cell geometry (no solve)
% import com.comsol.model.*
% import com.comsol.model.util.*
% ModelUtil.showProgress(true);
% ModelUtil.clear();
% clear model
% model = ModelUtil.create('model');
% [model,P] = buildRibUnitCell_LN(model,P);
% mphlaunch(model);

%% Single solve
% Output folders are those of the original LN rib runs. The mechanical file
% names keep the legacy format, so cached results under test/LN_ribUnitCell/
% are reused when the parameters match.
currentDate = datestr(now,'mmddyyyy');
if P.run_optical
    P.datLoc = [fullfile('.','test','LN_ribUnitCell_optical',currentDate),filesep];
    bds = solveOpticalBands(P);
else
    P.datLoc = [fullfile('.','test','LN_ribUnitCell',currentDate),filesep];
    bds = solveBands_noSym(P);
end
