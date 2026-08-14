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
% READ THIS BEFORE RELYING ON IT: with the run currently wired up at the bottom
% of this script, setting this to 1 changes nothing. Two separate reasons:
%
%   - The active path is solveOpticalBands with bandStructureDim = 3, which
%     dispatches 'boomerang' to buildBoomerangUnitCell_2D. That builder does not
%     read add_airDisk at all - only buildBoomerangUnitCell,
%     buildBoomerangUnitCellStrip_v2 and buildBoomerangStrip_3D do. So the flag
%     is silently ignored on the optical path.
%   - The mechanical path (solveBands, commented out below) reaches
%     buildBoomerangUnitCell, which DOES implement it - but runBands_2D:139
%     force-disables it, deliberately, because Solid Mechanics there is applied
%     to domain 1 only and an air domain would get no material and could
%     renumber which domain is the diamond.
%
% So this is a declaration of intent that becomes live once the air-disk block
% is ported into buildBoomerangUnitCell_2D (the builder the optical runs
% actually use), which is where an air region belongs anyway.
P.add_airDisk = 1;          % 1 = add the air region, 0 = solid only
P.airDiskHeight = 4000e-9;  % total air region height in m, measured from the
                            % base of the remaining solid (4 um ~ 2.7 free-space
                            % wavelengths at the 200 THz target below)

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
%   solveOpticalBands -> runOpticalBand_3D -> buildBoomerangUnitCell_2D
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

% KNOWN BLOCKER on the dim = 3 path selected above. runOpticalBand_3D reads
% P.zEnd for its SymmetryPlane selection, but buildBoomerangUnitCell_2D - the
% builder the 'boomerang' branch dispatches to - writes only xEnd1/xEnd2/yEnd1/
% yEnd2 and never sets zEnd, so this dies with "Unrecognized field name 'zEnd'"
% before it solves. The MECHANICAL builder buildBoomerangUnitCell does set it,
% which is why the phononic run is unaffected. Fix either by having the 2D
% builder record zEnd (a bndindex lookup at z = 0, as the mechanical one does)
% or by guarding that SymmetryPlane on isfield(P,'zEnd').
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
datLoc = ['.\test\boomerang\',currentDate,'\'];
P.datLoc = datLoc;
% %% solving mechanical modes
% bds = solveBands(P);
%% solving optical bands
bds = solveOpticalBands(P);