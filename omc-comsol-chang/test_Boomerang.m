clear all; clc; close all;
clear P;

%% unit cell params 
P.xsect = 'rect'; 
P.beamMat = 'diamond';                  % beam material name
P.celltype = 'boomerang';                   % specify the cell type
P.unitcell = 'hexagonal';                  % specify the shape of the unit cell
P.a = 700e-9;              % lattice constant
P.w = 140e-9;              % unit cell width (along x)
P.r = 177e-9;              % unit cell height (along y)
P.th = 220e-9;             % height (along x) of cross (for celltype = 'hollow')
                            % or of inner block (for celltype = 'solid')                   
P.r1 = 10e-9;             % width (along y) of cross (for celltype = 'hollow')
                            % or of inner block (for celltype = 'solid')
P.r2 = 10e-9;              % height (along x) of each leg in cross (for celltype = 'hollow')
                            % or of outer fins (for celltype = 'solid')
% --- air disk (optional geometry) ---------------------------------------
% Toggles the air region buildBoomerangUnitCell.m adds on top of the slab,
% ported from comsol_templates/BoomerangWithAirDisk.m: a second work plane with
% the SAME rhombic footprint as the slab base, extruded upward and left
% uncomposed, so Form Union partitions the overlap and the etched hole plus the
% space above the slab both become air while the diamond keeps its own domain.
% airDiskHeight is the TOTAL height from the base of the remaining solid (z = 0
% when P.mbevenz is nonzero, otherwise -th/2), not the clearance above the top
% face. P.airDiskH is accepted as an alias.
%
% WHERE THIS IS LIVE: on the optical path with bandStructureDim = 3.
% runOpticalBand_3D:94 dispatches 'boomerang' to buildBoomerangUnitCell (NOT
% buildBoomerangUnitCell_2D), and that builder implements the air disk at its
% line 150. The low-reflecting Scattering condition is placed on the top face of
% this region, so the height set here IS the clearance that governs how well
% out-of-plane radiation is absorbed.
%
% It is inert on the other two paths:
%   - bandStructureDim = 2 reaches runOpticalBand_2D, an in-plane cross-section
%     with no thickness and no air region at all, so neither this nor P.th is
%     read.
%   - The mechanical path (solveBands, commented out below) also reaches
%     buildBoomerangUnitCell, but runBands_2D:139 force-disables the air disk
%     deliberately: Solid Mechanics there is applied to domain 1 only, and an
%     air domain would get no material and could renumber which domain is the
%     diamond.
P.add_airDisk = 1;          % 1 = add the air region, 0 = solid only
P.airDiskHeight = 1500e-9;  % total air region height in m, measured from the
                            % base of the remaining solid. With mbevenz nonzero
                            % the half-slab spans z = 0 to th/2 = 110 nm, so this
                            % leaves 1390 nm of clearance above the slab - about
                            % 4.5 evanescent decay lengths at the M point
                            % (1/kappa = 310 nm at 1550 nm for a = 700 nm), which
                            % is past the point where more air changes anything.
                            % Was 4000 nm, i.e. 12.5 decay lengths - accurate but
                            % paying for roughly twice the air elements needed.
                            % Distance only kills EVANESCENT orders; propagating
                            % diffraction orders hit the boundary at a fixed
                            % angle set by the lattice, so if Q accuracy matters
                            % the lever is a PML, not a taller air disk.

P.nperiod = 1;  % no. of periods to simulate for
P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
P.mbevenz = 1;      % 1 to find even mechanical mode about z

P.kpts = 9;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 10;                           % no. of bands to solve for

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
P.TwoSymPlanes = 0;             % if we are solving for band structures with two symmetry planes
P.mbeveny = 0;                          % 1 to find even mechanical mode about y
P.mbevenz = 1;                          % 1 to find even mechanical mode about z
P.freq = 0;                             % target frequency - set to 0 for bandstructure simulations
P.meshSize = 3;                         % mesh quality for mechanical simulations
P.fixed_bc = 0;                       % 1 to fixed the boundaries for xz planes at y = +/- w/2

P.anisoMat = 1;
P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg 
                                        % from <100> inplane direction about <100> surface normal
%% optical simulation parameters
% Chain for celltype 'boomerang' with bandStructureDim = 3:
%   solveOpticalBands -> runOpticalBand_3D -> buildBoomerangUnitCell
%                     -> findGaps_optical
% Of everything that chain reads, P.optical_freq below is the only input the
% blocks above do not already supply. The rest is shared with the mechanical
% run: a, w, r, r1, r2 (geometry), th (3D path only), beamMat, kpts, nbands,
% meshSize, max_dof, plotgeom, savebndplot, saveMPH, and optional prefname.
% Note P.freq above (Hz) is NOT used here - the optical solvers read
% P.optical_freq instead.
%
% 3 -> runOpticalBand_3D, 2 -> runOpticalBand_2D, anything else -> _1D. There is
% no validation, so a typo silently selects the 1D strip solver. The suffix is
% the MODEL dimension, not the k-path: 2 and 3 sweep the identical hexagonal
% Gamma-X-M-Gamma circuit and differ in physics (2 = ewfd in-plane
% cross-section, no thickness and no radiation; 3 = full-vector emw slab with a
% PML and a z-symmetry plane, and the only one that reads P.th).
P.bandStructureDim = 3;
% Target optical mid-gap frequency as a BARE NUMBER IN THz, not Hz - the
% repo-wide convention, and the opposite of P.freq above. Sets the eigenvalue
% search shift. 200 THz is 1499 nm free-space; use 193.4 for 1550 nm.
% runOpticalBand_3D sets eigunit to 'THz' before applying it; runOpticalBand_2D
% never sets eigunit, so on that path the same number falls back to COMSOL's
% default frequency unit.
P.optical_freq = 200;                   % target optical mid-gap frequency [THz]

% NOT a blocker (this note used to claim otherwise). runOpticalBand_3D reads
% P.zEnd for its SymmetryPlane selection, and the 'boomerang' branch dispatches
% to buildBoomerangUnitCell - not buildBoomerangUnitCell_2D - which records
% zEnd at its line 249 via a bndindex lookup at z = 0. buildBoomerangUnitCell_2D
% would indeed be missing zEnd, but nothing routes to it from here.
%
% If you switch to dim = 2 instead: it gets through the solve but then errors in
% the plot, because the 2D branch of solveOpticalBands builds its title from a
% celltype chain that has no 'boomerang' case and no else, leaving bandtitle
% undefined. Set P.savebndplot = 0, or add 'boomerang' to that chain. The 1D/3D
% plot branch has an else fallback, so dim = 3 is fine on that count.
%
% Cache note: the optical fileBase for 'boomerang' encodes only a, r, w, r1, r2
% - NOT th - so two thicknesses collide on one cached filename, which matters
% exactly on this dim = 3 path since it is the only one that reads P.th.
%% define the maximum number of degree of freedom to limit the simulation time
P.max_dof = 3e6;                        % max # of degrees of freedom

% % debugging the unit cells 
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
% mphlaunch(model)
%% Single solve

currentDate = datestr(now,'mmddyyyy');
datLoc = [fullfile('.','test','boomerang',currentDate),filesep];
P.datLoc = datLoc;
% %% solving mechanical modes
% bds = solveBands(P);
%% solving optical bands
bds = solveOpticalBands(P);