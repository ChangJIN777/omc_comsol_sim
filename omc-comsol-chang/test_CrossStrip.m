clear all; clc; close all;
clear P;

%% strip params
% Geometry: a STRIP (supercell) of P.ncell square cross cells, built by
% buildCrossStrip.m. It is to test_CrossUnitCell.m what buildHoleStrip_3D is to
% a single hole cell -- but on a SQUARE lattice, not the hexagonal one.
%
% A base rectangle with its corner at [-a/2, 0] and size [a, Ly]
% (buildCrossStrip.m:187) is etched with P.ncell cross-shaped VOIDS, then
% extruded by th (:232). Each void is the union of two crossed bars, exactly as
% in the single cell:
%     horizontal arm :  |x| <= h/2 ,  |y - y_i| <= w/2      (:201)
%     vertical   arm :  |x| <= w/2 ,  |y - y_i| <= h/2      (:206)
% and all of them are SUBTRACTED from the base rectangle by the Compose at
% :214, so the solid is what remains around the crosses.
%
% Square lattice -- every cross sits at x = 0, spaced a apart in y (:157):
%     y_i = b_wvg + (i - 1/2)*a + b*(i == 1),   i = 1 .. P.ncell
%     Ly  = b_wvg + P.ncell*a                                   (:159)
% Contrast buildHoleStrip_3D, whose hexagonal lattice has a row pitch of
% sqrt(3)*a/2 and alternates x = 0 with x = +-a/2, so that half its holes are
% deliberately cut in two by the Floquet planes. Here nothing is cut by
% x = +-a/2: with h < a each void is entirely interior to its own cell, and the
% two x faces are plain rectangles.
%
% NO MIRROR AT y = 0. buildCrossStrip borrows buildHoleStrip_3D's half-strip
% FOOTPRINT (y from 0 to Ly rather than -Ly/2 to +Ly/2) but does not treat
% y = 0 as a symmetry plane. Both y faces are emitted -- y = 0 and y = Ly go
% into the 'yboundaries' cumulative selection (:309) and come back as P.yEnd1
% and P.yEnd2 (:336,:337) -- and the choice of boundary condition is left to
% the caller. With b = b_wvg = 0 (the defaults, and what is set below) y = 0
% and y = Ly land exactly on cell edges, one lattice vector apart, so they are
% a legitimate periodic pair and the strip is exactly P.ncell whole cells.
%
% The narrowest solid feature is unchanged from the single cell: the ligament
% (a - h)/2 = 15.5 nm between an arm end and the cell edge. In the strip it
% appears at the y = 0 and y = Ly faces and, back to back across each interior
% cell boundary, as an (a - h) = 31 nm bridge between neighbouring crosses.
% Watch the mesh there.
%
% Physical parameters below are IDENTICAL to test_CrossUnitCell.m so that the
% strip is directly comparable with the single cell.
P.xsect = 'rect';
P.beamMat = 'diamond';                  % beam material name
P.celltype = 'cross_strip';             % specify the cell type
P.unitcell = 'square';                  % specify the shape of the unit cell
P.a = 914e-9;       % lattice constant along BOTH x and y. Sets the x width of
                    % the footprint (:187) and the y pitch of the cells (:157).
P.h = 847e-9;       % LENGTH of each cross arm. recH has size [h w] (x by y)
                    % and recV the transpose [w h] (:201,:206).
                    % h < a is REQUIRED and enforced -- buildCrossStrip errors
                    % with :armTooLong (:387) rather than letting COMSOL build
                    % a chain of disconnected corner islands.
P.w = 184e-9;       % WIDTH of each cross arm (the etched gap width).
P.th = 250e-9;      % slab THICKNESS along z. The work plane sits at z = -th/2
                    % and the profile is extruded a distance th (:180,:232).
P.r1 = 10e-9;       % fillet radius, applied via fil1.set('radius',r1) (:543)
P.r2 = 10e-9;       % fillet radius, applied via fil2.set('radius',r2) (:557)
                    % See the KNOWN ISSUE block at buildCrossStrip.m:95. These
                    % are the LEGACY selections from buildCrossUnitCell,
                    % reproduced unchanged so results stay comparable with
                    % already-simulated designs. At these parameters that means
                    % two warnings you should expect and not be alarmed by:
                    %   :filletR2NoOp        disksel2 spans [426.5, 436.5] nm
                    %                        but the arm tips sit at 479.0 nm,
                    %                        so r2 does nothing
                    %   :filletSelectionBleed disksel1 spans [208, 624] nm,
                    %                        which reaches into the NEIGHBOURING
                    %                        cells (they are only a = 894 nm
                    %                        away), so a cell's fillet rounds
                    %                        corners belonging to its neighbour
                    % The second one only exists in the strip. Inspect the
                    % corners in the GUI before trusting a long strip -- see
                    % the note under P.ncell.
P.ncell = 7;        % NUMBER OF CROSS CELLS along the strip (y). This is the
                    % field that replaces buildHoleStrip_3D's 13 hardcoded
                    % circles. Validated at :381 as a positive integer.
                    % Kept at 3 deliberately for a first run: it is cheap, and
                    % it is the smallest N that has a cell with neighbours on
                    % BOTH sides, which is what makes the :filletSelectionBleed
                    % case above visible. With N = 3:
                    %     y_i = 447, 1341, 2235 nm      Ly = 2682 nm
                    % Raise it only after the GUI shows the fillets are sane.
P.b = 0;            % extra y shift applied to CELL 1 ONLY, mirroring the role
                    % of P.b in buildHoleStrip_3D. Range-checked at :412 so it
                    % cannot silently push cell 1 outside the footprint.
P.b_wvg = 0;        % y offset applied to the WHOLE array, widening the gap
                    % between y = 0 and the first cell. Leave at 0 unless you
                    % want a line defect: with b_wvg ~= 0 the y = 0 face is no
                    % longer a lattice translation of the y = Ly face, so a
                    % Floquet pair across them is no longer the square lattice.
P.nperiod = 1;      % no. of periods to simulate for -- still means periods
                    % along X. buildCrossStrip does not read it; P.ncell counts
                    % the cells along y.
P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell. Not read by
                    % buildCrossStrip (with h < a the void is always interior);
                    % kept for parity with the other test scripts.
P.mbevenz = 1;      % 1 to find even mechanical mode about z
                    % (nonzero also halves the cell in z: the block below
                    % z = 0 is subtracted, buildCrossStrip.m:240-262, and the
                    % z = 0 plane comes back as geom1_ZsymSel / P.bndSel.Zsym)

P.kpts = 5;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 70;                          % no. of bands to solve for

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
%
% For this strip, 0 is the right choice for a second reason: the y faces are a
% PERIODIC pair here, not a symmetry plane (see the NO MIRROR note above), so
% driving mbeveny would apply a symmetry condition to boundaries that are meant
% to carry a Floquet condition.
P.TwoSymPlanes = 0;

% Which plane the single symmetry condition applies to, when TwoSymPlanes = 0.
%   1 -> z is the symmetry plane, solveBands drives P.mbevenz = +1 then -1
%   0 -> y is the symmetry plane, solveBands drives P.mbeveny instead
% Set to 1: z is the only genuine mirror plane of this strip.
P.zSymCondition = 1;

P.solveasym = 1;                        % 1 to solve for antisymmetric bands
P.completeBandGaps = 1;                 % 1 to plot complete bandgaps (across all symmetries)
P.plotgeom = 1;                         % 1 to plot the geometry
P.savedat = 1;                          % 1 to save data structures
P.savebndplot = 1;                      % 1 to save bandstructure plot
P.saveplots = 0;                        % 1 to save displacement and strain profiles
P.saveMPH = 0;
P.bandStruct_2D = 0;                    % 0 to simulate 1D band structures.
                                        % DIFFERENT from test_CrossUnitCell.m,
                                        % which sets 1. A strip is periodic in
                                        % X ONLY -- the supercell spans the
                                        % whole structure in y -- so the band
                                        % structure is a 1D sweep along kx.
                                        % 0 also routes solveBands to runBands
                                        % rather than runBands_2D, which is the
                                        % solver buildCrossStrip was written
                                        % for (it emits exactly the
                                        % geom1_xboundaries_bnd /
                                        % geom1_yboundaries_bnd / geom1_ZsymSel
                                        % that runBands.m:232,:257,:262
                                        % consumes).

%% mechanical simulation parameters
% solid mechanics solver parameters
P.mbeveny = 0;                          % 1 to find even mechanical mode about y
P.mbevenz = 1;                          % 1 to find even mechanical mode about z
P.freq = 0;                             % target frequency - set to 0 for bandstructure simulations
P.meshSize = 4;                         % mesh quality for mechanical simulations
P.fixed_bc = 0;                         % 1 to fixed the boundaries for xz planes at y = +/- w/2
                                        % (this is the only path that reads
                                        % P.xEnd1/P.xEnd2/P.yEnd1/P.yEnd2 --
                                        % runBands.m:212-219)

P.anisoMat = 1;
P.rxtal = 45;                           % ccw rotation of elasticity matrix in deg
                                        % from <100> inplane direction about <100> surface normal

%% define the maximum number of degree of freedom to limit the simulation time
P.max_dof = 3e6;                        % max # of degrees of freedom

%% WHICH PATH THIS SCRIPT TAKES
% 1 -> build the geometry only, and open it in the COMSOL GUI. THIS IS THE
%      WORKING DEFAULT and it runs today.
% 0 -> run the full band-structure solve through solveBands.
%
% READ THIS BEFORE SETTING IT TO 0. The 'cross_strip' celltype is NOT yet
% wired into the solver chain, so a solve cannot work until two existing files
% are edited (buildCrossStrip.m was delivered as a geometry builder only):
%
%   1. solveBands.m  -- add a 'cross_strip' branch to the fileBase chain,
%      alongside the 'cross' branch at :39. Without it P.fileBase is never
%      assigned and solveBands dies at :134 (fBase = P.fileBase) with an
%      unhelpful "Unrecognized field name 'fileBase'".
%
%   2. runBands.m    -- add a 'cross_strip' branch to the geometry dispatch at
%      :113-125 calling buildCrossStrip(model,P). This one matters MORE than it
%      looks: that if/elseif chain ends in a bare `else` that falls through to
%      buildBoomerangStrip_3D, so an unrecognised celltype does not raise
%      "unknown celltype" -- it quietly tries to build a BOOMERANG strip
%      instead, and whatever it then complains about will point you at the
%      wrong file.
%
%   CreateFileBase.m does NOT need touching: solveBands builds fBase with its
%   own inline chain and never calls it.
%
% The guard below checks for both branches and, if either is missing, stops
% with those instructions rather than letting you hit the errors above. Once
% you have added them the guard passes on its own.
P.geomOnly = 0;

if P.geomOnly
    %% GEOMETRY ONLY -- build and inspect in the COMSOL GUI
    % COMSOL with LiveLink for MATLAB must already be running.
    import com.comsol.model.*
    import com.comsol.model.util.*

    ModelUtil.showProgress(true);
    ModelUtil.clear();
    clear model

    % create COMSOL model named 'model' from which COMSOL methods can be
    % called, e.g. model.save
    model = ModelUtil.create('model');

    [model,P] = buildCrossStrip(model,P);

    % Expect two warnings at these parameters -- :filletR2NoOp and
    % :filletSelectionBleed. Both are the reproduced legacy fillet selections
    % reporting themselves; neither changes the geometry. See the P.r1/P.r2
    % comments above.
    fprintf('\nbuildCrossStrip: %d cells, Ly = %.1f nm\n', P.ncell, P.Ly*1e9);
    fprintf('  cell centres [nm] : %s\n', num2str(P.ycentres*1e9, '%.1f  '));
    fprintf('  xEnd1 / xEnd2     : [%s] / [%s]\n', num2str(P.xEnd1), num2str(P.xEnd2));
    fprintf('  yEnd1 / yEnd2     : [%s] / [%s]\n', num2str(P.yEnd1), num2str(P.yEnd2));
    fprintf('  zEnd              : [%s]\n', num2str(P.zEnd));
    if isfield(P,'bndSel') && isfield(P.bndSel,'Zsym')
        fprintf('  bndSel.Zsym       : [%s]\n', num2str(P.bndSel.Zsym));
    end
    fprintf(['\nInspect the filleted corners around cell 2 (the one with ' ...
             'neighbours\non both sides) before raising P.ncell.\n\n']);

    figure;
    mphgeom(model);
    mphlaunch(model);
else
    %% FULL SOLVE -- requires the wiring described above
    wiringGaps = {};
    solveBandsSrc = fileread(which('solveBands'));
    if ~contains(solveBandsSrc, P.celltype)
        wiringGaps{end+1} = ['solveBands.m: add a ''',P.celltype, ...
            ''' branch to the fileBase chain (next to the ''cross'' branch ' ...
            'at :39), or solveBands dies at :134 on P.fileBase.'];
    end
    runBandsSrc = fileread(which('runBands'));
    if ~contains(runBandsSrc, P.celltype)
        wiringGaps{end+1} = ['runBands.m: add a ''',P.celltype, ...
            ''' branch to the geometry dispatch at :113-125 calling ' ...
            'buildCrossStrip(model,P). The chain''s bare else silently ' ...
            'falls through to buildBoomerangStrip_3D.'];
    end
    if ~isempty(wiringGaps)
        error('test_CrossStrip:solverNotWired', ...
            ['The ''%s'' celltype is not wired into the solver chain yet:\n' ...
             '  - %s\n' ...
             'Set P.geomOnly = 1 to build and inspect the geometry instead.'], ...
            P.celltype, strjoin(wiringGaps, '\n  - '));
    end

    currentDate = datestr(now, 'mmddyyyy');
    datLoc = [fullfile('.','test','cross_strip',currentDate),filesep];
    P.datLoc = datLoc;
    bds = solveBands(P);
end
