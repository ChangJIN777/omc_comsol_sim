clear all; clc; close all;
clear P;

%% unit cell params
% Geometry: a square slab with a cross-shaped VOID etched through it. The two
% arm rectangles are SUBTRACTED from the cell square -- buildCrossUnitCell.m:53
% composes 'r_ucell-r1-r2' -- so the solid is what remains around the cross,
% and the narrowest solid feature is the ligament between an arm end and the
% cell edge, (a-h)/2.
%
% In the xy plane the void is the union of two crossed bars, both centred on
% the cell:
%     horizontal arm :  |x| <= h/2 ,  |y| <= w/2
%     vertical   arm :  |x| <= w/2 ,  |y| <= h/2
% so the solid is the square minus that plus-shape. With h < a the void does
% not reach the boundary, and the solid stays connected all the way round
% through a ligament of (a-h)/2 at each of the four arm ends. At h = a those
% ligaments vanish and the cell falls into four disconnected corner islands.
%
% All dimensions verified against buildCrossUnitCell.m, not inherited from the
% previous comments here, which described a different (solid-cross) geometry
% and had a, h, w, th, r1 and r2 all mislabelled.
P.xsect = 'rect';
P.beamMat = 'diamond';                  % beam material name
P.celltype = 'cross';                   % specify the cell type
P.unitcell = 'square';                  % specify the shape of the unit cell
P.a = 500e-9;       % unit cell side - the lattice constant along BOTH x and y.
                    % ucellPlane is a square of size [a a] (buildCrossUnitCell.m:38).
P.h = 420e-9;       % LENGTH of each cross arm. rec_1 has size [h w] (x by y) and
                    % rec_2 the transpose [w h], giving two crossed bars (:43,:48).
                    % h < a is required, or the arms reach the cell boundary and
                    % split the solid into four disconnected corner islands.
P.w = 170e-9;        % WIDTH of each cross arm (the etched gap width).
P.th = 250e-9;      % slab THICKNESS along z. The work plane sits at z = -th/2
                    % and the profile is extruded a distance th (:33, :65).
P.r1 = 10e-9;       % fillet radius applied by addFillet - fil1.set('radius',r1)
P.r2 = 10e-9;       % fillet radius applied by addFillet - fil2.set('radius',r2)
                    % (buildCrossUnitCell.m:137,:147). These round the inner
                    % corners of the void; they are NOT cross dimensions.
P.nperiod = 1;  % no. of periods to simulate for
P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell
P.mbevenz = 1;      % 1 to find even mechanical mode about z
                    % (nonzero also halves the cell in z: the block below
                    % z = 0 is subtracted, buildCrossUnitCell.m:72-92)

P.kpts = 5;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 10;                           % no. of bands to solve for

%% symmetry decomposition (required by solveBands)
% solveBands gates on these two BEFORE any solve; neither has a default, so
% leaving them out errors out immediately at solveBands.m:138 / :186.
%
% TwoSymPlanes = 0 -> exploit ONE symmetry plane, solving two sectors
%                     (sym, then asym when P.solveasym = 1). This matches the
%                     mbeveny/mbevenz settings below: y is not a symmetry
%                     plane here (mbeveny = 0), z is (mbevenz = 1).
%              = 1 -> exploit BOTH y and z, solving four sectors
%                     (symy_symz, symy_asymz, asymy_symz, asymy_asymz) at 4x
%                     the solve cost. Only correct if the cell really is
%                     symmetric about both planes.
P.TwoSymPlanes = 0;

% Which plane the single symmetry condition applies to, when TwoSymPlanes = 0.
%   1 -> z is the symmetry plane, solveBands drives P.mbevenz = +1 then -1
%   0 -> y is the symmetry plane, solveBands drives P.mbeveny instead
% Set to 1 for this thin (w = 50 nm along z) square cross cell.
P.zSymCondition = 1;

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
% buildCrossUnitCell(model,P);
% mphlaunch(model);
%% Single solvecurrentDate = datestr(now,'mmddyyyy');
currentDate = datestr(now, 'mmddyyyy');
datLoc = [fullfile('.','test','cross',currentDate),filesep];
P.datLoc = datLoc;
bds = solveBands(P);