clear all; clc; close all;
clear P;

%% strip params
% Geometry: a STRIP (supercell) of P.ncell square cross cells, built by
% buildCrossStrip.m. It is to test_CrossUnitCell.m what buildHoleStrip_3D is to
% a single hole cell -- but on a SQUARE lattice, not the hexagonal one.
%
% A base rectangle with its corner at [-a/2, yLo] and size [a, Ly]
% (buildCrossStrip.m:256) is etched with P.ncell cross-shaped VOIDS, then
% extruded by th (:302). Each void is the union of two crossed bars, exactly as
% in the single cell:
%     horizontal arm :  |x| <= h/2 ,  |y - y_i| <= w/2      (:271)
%     vertical   arm :  |x| <= w/2 ,  |y - y_i| <= h/2      (:276)
% and all of them are SUBTRACTED from the base rectangle by the Compose at
% :285, so the solid is what remains around the crosses.
%
% Square lattice -- every cross sits at x = 0, spaced a apart in y (:213):
%     y_i  = b_wvg + (i - 1/2)*a + b*(i == 1),  i = 1 .. P.ncell
%     yTop = b_wvg + P.ncell*a                                  (:215)
%     yLo  = 0, or y_1 when P.cutBottomHalfCell (:221-225)
%     Ly   = yTop - yLo                                         (:226)
% Note Ly is the EXTENT, not the y of the top face: once the footprint starts
% above 0 the two stop being the same number.
% Contrast buildHoleStrip_3D, whose hexagonal lattice has a row pitch of
% sqrt(3)*a/2 and alternates x = 0 with x = +-a/2, so that half its holes are
% deliberately cut in two by the Floquet planes. Here nothing is cut by
% x = +-a/2: with h < a each void is entirely interior to its own cell, and the
% two x faces are plain rectangles.
%
% NO MIRROR AT y = 0. buildCrossStrip borrows buildHoleStrip_3D's half-strip
% FOOTPRINT (y from yLo to yTop rather than -Ly/2 to +Ly/2) but does not treat
% y = 0 as a symmetry plane. Both y faces are emitted -- y = yLo and y = yTop
% go into the 'yboundaries' cumulative selection (:384) and come back as
% P.yEnd1 and P.yEnd2 (:465,:466) -- and the choice of boundary condition is
% left to the caller. THIS SCRIPT uses that freedom twice, see
% P.cutBottomHalfCell and P.fixed_bc / P.fixed_faces below: the bottom face is
% cut through the middle of cell 1 and the top face is clamped.
%
% WHAT THAT MEANS FOR PERIODICITY. With b = b_wvg = 0 and no truncation the two
% y faces sit exactly one lattice vector apart and are a legitimate periodic
% pair. With P.cutBottomHalfCell = 1 they are half a lattice vector apart and
% no longer a pair -- which is fine here, because this is a 1D kx sweep
% (P.bandStruct_2D = 0): runBands only ever Floquet-pairs the two X faces
% (runBands.m:334-338). The y faces get nothing at all unless P.fixed_bc asks
% for it, and P.mbeveny = 0 with the cross_strip guard at runBands.m:396-399
% keeps the y symmetry/antisymmetry features switched off.
%
% THIN FEATURES / MESH. The narrowest solid feature is the ligament
% (a - h)/2 = 33.5 nm at these parameters (a = 914, h = 847 nm) between an arm
% end and the cell edge. It appears at both y faces and, back to back across
% each interior cell boundary, as an (a - h) = 67 nm bridge between
% neighbouring crosses. Watch the mesh there.
%   The half-cell truncation does not change that width: the cut at y = y_1
% passes through the WIDEST part of cell 1's void (the horizontal arm spans
% |x| <= h/2 there), so the bottom face is two rectangles of exactly the same
% (a - h)/2 = 33.5 nm width, by th/2 = 125 nm tall with mbevenz cutting z < 0 -
% the ligament that is already in the model, just exposed as a boundary rather
% than a bridge.
%   What DOES get thinner is a side effect of the legacy fillets: the cut
% corners land inside disksel2, so r2 = 10 nm rounds them and each bottom face
% is left with (a - h)/2 - r2 = 23.5 nm of flat plus a 10 nm arc. That 23.5 nm
% is the thinnest flat feature in the model. See the P.cutBottomHalfCell and
% P.r2 comments below; P.r2 = 0 removes it at no other cost here.
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
                    % and recV the transpose [w h] (:270,:275).
                    % h < a is REQUIRED and enforced -- buildCrossStrip errors
                    % with :armTooLong (:516) rather than letting COMSOL build
                    % a chain of disconnected corner islands.
P.w = 184e-9;       % WIDTH of each cross arm (the etched gap width).
P.th = 250e-9;      % slab THICKNESS along z. The work plane sits at z = -th/2
                    % and the profile is extruded a distance th (:250,:302).
P.r1 = 10e-9;       % fillet radius, applied via fil1.set('radius',r1) (:719)
P.r2 = 10e-9;       % fillet radius, applied via fil2.set('radius',r2) (:733)
                    % See the KNOWN ISSUE block at buildCrossStrip.m:150. These
                    % are the LEGACY selections from buildCrossUnitCell,
                    % reproduced unchanged so results stay comparable with
                    % already-simulated designs. AT THESE PARAMETERS
                    % (a = 914, h = 847, w = 184 nm) expect exactly these
                    % warnings and do not be alarmed by them:
                    %   :filletR2NoOp        disksel2 spans [418.5, 428.5] nm
                    %                        but the arm tips sit at 433.4 nm,
                    %                        so r2 does nothing TO THE TIPS
                    %   :filletSelectionBleed fires only because
                    %                        P.cutBottomHalfCell = 1 below:
                    %                        the two vertices the cut leaves at
                    %                        (+-h/2, y_1) sit at exactly
                    %                        h/2 = 423.5 nm, inside disksel2, so
                    %                        r2 rounds those instead
                    % disksel1's annulus is only [92, 276] nm here, well short
                    % of the a = 914 nm to the next cell, so the neighbour-cell
                    % bleed that the KNOWN ISSUE block warns about does NOT
                    % happen at this w. It would if P.w grew past about a/3.
                    % Inspect the corners in the GUI before trusting a long
                    % strip -- see the note under P.ncell.
P.ncell = 7;        % NUMBER OF CROSS CELLS along the strip (y). This is the
                    % field that replaces buildHoleStrip_3D's 13 hardcoded
                    % circles. Validated at :510 as a positive integer.
                    % Any N >= 3 has a cell with neighbours on BOTH sides,
                    % which is what makes the :filletSelectionBleed case above
                    % visible; N = 3 is the cheap version to start from.
                    % With N = 7, a = 914 nm, b = b_wvg = 0:
                    %     y_i  = 457, 1371, 2285, 3199, 4113, 5027, 5941 nm
                    %     yTop = 6398 nm
                    %     Ly   = 6398 nm             (cutBottomHalfCell = 0)
                    %     Ly   = 5941 nm = 6.5*a     (cutBottomHalfCell = 1,
                    %                                 yLo = y_1 = 457 nm)
                    % Raise it only after the GUI shows the fillets are sane.
P.b = 0;            % extra y shift applied to CELL 1 ONLY, mirroring the role
                    % of P.b in buildHoleStrip_3D. Range-checked at :547 so it
                    % cannot silently push cell 1 outside the footprint.
P.b_wvg = 0;        % y offset applied to the WHOLE array, widening the gap
                    % between y = 0 and the first cell. Leave at 0 unless you
                    % want a line defect: with b_wvg ~= 0 the lower face is no
                    % longer a lattice translation of the y = yTop face, so a
                    % Floquet pair across them is no longer the square lattice.
P.cutBottomHalfCell = 1;
                    % REMOVE THE LOWER HALF OF THE BOTTOM UNIT CELL.
                    % 0 (default, and what every other caller of
                    %   buildCrossStrip gets) -> footprint starts at y = 0, so
                    %   the strip is exactly P.ncell whole cells.
                    % 1 -> footprint starts at yLo = y_1, the CENTRE of the
                    %   first cross, so the strip is P.ncell - 1/2 cells long
                    %   (6.5 cells, Ly = 5941 nm at the numbers above) and the
                    %   bottom face cuts cell 1 through the middle.
                    % Only the base rectangle's lower extent and Ly change:
                    % the cell centres y_i and every cross void stay at the
                    % SAME absolute y, so cells 2..N are untouched and results
                    % stay comparable cell for cell. The y = yLo face is still
                    % emitted into 'yboundaries' and still returned as
                    % P.yEnd1 -- as TWO boundary indices now, because the cut
                    % plane crosses the void and leaves one ligament face
                    % either side of it.
                    % Validated in buildCrossStrip:readCrossStripParams --
                    % :truncationLeavesNothing if the cut is at or past the end
                    % of the strip, :truncationCutsSecondCell if a large P.b
                    % has moved cell 1 far enough up that the cut would land
                    % inside cell 2's void.
                    % ONE EXTRA WARNING TO EXPECT at these parameters: the cut
                    % leaves two new vertices at (+-h/2, y_1), i.e. at radius
                    % exactly h/2 = 423.5 nm from the cell-1 centre, which is
                    % the middle of the legacy disksel2 annulus
                    % [418.5, 428.5] nm. So P.r2 -- a no-op on the arm tips,
                    % see :filletR2NoOp above -- DOES round those two corners,
                    % and :filletSelectionBleed now names them.
                    % MESH CONSEQUENCE, the one thing here that is genuinely
                    % new: that fillet eats into the bottom face, which goes
                    % from (a-h)/2 = 33.5 nm of flat to
                    % (a-h)/2 - r2 = 23.5 nm of flat plus a 10 nm arc. 23.5 nm
                    % is then the thinnest flat feature in the model. If you
                    % would rather keep the full 33.5 nm, set P.r2 = 0: at THIS
                    % parameter set that changes nothing else at all, because
                    % disksel2 selects nothing without the truncation (that is
                    % exactly what :filletR2NoOp is telling you), and
                    % addCrossFillet skips disksel2/fil2 entirely when r2 = 0.
P.nperiod = 1;      % no. of periods to simulate for -- still means periods
                    % along X. buildCrossStrip does not read it; P.ncell counts
                    % the cells along y.
P.holeatedge = 0;   % 1/0 for hole at edge/center of unit cell. Not read by
                    % buildCrossStrip (with h < a the void is always interior);
                    % kept for parity with the other test scripts.
P.mbevenz = 1;      % 1 to find even mechanical mode about z
                    % (nonzero also halves the cell in z: the block below
                    % z = 0 is subtracted, buildCrossStrip.m:310-345, and the
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
                                        % that runBands.m:337,:362,:373
                                        % consumes).

%% mechanical simulation parameters
% solid mechanics solver parameters
P.mbeveny = 0;                          % 1 to find even mechanical mode about y
P.mbevenz = 1;                          % 1 to find even mechanical mode about z
P.freq = 0;                             % target frequency - set to 0 for bandstructure simulations
P.meshSize = 4;                         % mesh quality for mechanical simulations
P.fixed_bc = 1;                         % 1 to apply a Fixed (zero displacement)
                                        % condition (runBands.m:229-331). Set
                                        % to 1 here to CLAMP THE TOP y SURFACE
                                        % of the strip, which is what makes the
                                        % supercell behave like a shield
                                        % anchored to bulk rather than a free
                                        % ribbon.
                                        %
                                        % WHICH boundaries get clamped is NOT
                                        % set here and is NOT an index list.
                                        % buildCrossStrip builds a BoxSelection
                                        % around the y = yTop face
                                        % (buildCrossStrip.m:416) and hands back
                                        % its NAME in
                                        %   P.bndSel.yFixed = 'geom1_yFixedSel'
                                        % runBands hangs the Fixed feature on
                                        % that name with selection.named, the
                                        % same mechanism that has always
                                        % carried the Floquet pair
                                        % ('geom1_xboundaries_bnd') and the z
                                        % symmetry ('geom1_ZsymSel') -- see
                                        % runBands.m:337,:373. A name
                                        % re-resolves every time the geometry
                                        % is rebuilt, so unlike a captured
                                        % index list it cannot go stale when
                                        % P.ncell, P.mbevenz or a fillet radius
                                        % changes.
                                        %
                                        % Measured on COMSOL 6.3 at these
                                        % parameters, 'geom1_yFixedSel'
                                        % resolves to boundary 5, spanning
                                        % x = -457..457 nm, y = 6398 nm,
                                        % z = 0..125 nm -- the top face and
                                        % nothing else -- for all four
                                        % combinations of cutBottomHalfCell and
                                        % mbevenz.
                                        %
                                        % Deliberately the top face ALONE, not
                                        % all four faces and not
                                        % 'geom1_yboundaries_bnd' (which holds
                                        % BOTH y faces, because the
                                        % symmetry/antisymmetry path needs the
                                        % pair):
                                        %   - the x faces must keep their
                                        %     Floquet/periodic condition
                                        %     (runBands.m:334-338); clamping
                                        %     one half of a periodic pair
                                        %     over-constrains it
                                        %   - the y = yLo face is the cut face
                                        %     created by P.cutBottomHalfCell
                                        %     and is meant to stay free
                                        %   - z = 0 keeps its symmetry
                                        %     condition from P.mbevenz = 1
                                        % runBands errors with :fixedSelEmpty
                                        % if the named selection resolves to
                                        % nothing, rather than solving with an
                                        % empty Fixed feature and handing back
                                        % what look like converged
                                        % free-boundary modes.
P.fixed_faces = {'yEnd2'};              % NOT THE ACTIVE PATH FOR THIS SCRIPT.
                                        % P.fixed_faces is runBands' FALLBACK
                                        % selector, used only for celltypes
                                        % whose builder supplies no
                                        % P.bndSel.yFixed (buildCrossStrip
                                        % does, so this value is ignored).
                                        % Kept, set to the correct value, so
                                        % that the fallback would still clamp
                                        % the right face if the named selection
                                        % ever went away. Accepts any nonempty
                                        % subset of
                                        % {'xEnd1','xEnd2','yEnd1','yEnd2'},
                                        % char or cellstr; unset defaults to
                                        % {'yEnd2'}.

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
% BOTH SOLVER BRANCHES NOW EXIST, so 0 no longer stops at the guard below:
%
%   1. solveBands.m:50  -- the 'cross_strip' branch of the fileBase chain,
%      alongside the 'cross' branch at :39. Without it P.fileBase is never
%      assigned and solveBands dies at :134 (fBase = P.fileBase) with an
%      unhelpful "Unrecognized field name 'fileBase'".
%
%   2. runBands.m:139   -- the 'cross_strip' branch of the geometry dispatch,
%      calling buildCrossStrip(model,P). That branch matters MORE than it
%      looks: the if/elseif chain ends in a bare `else` that falls through to
%      buildBoomerangStrip_3D, so an unrecognised celltype does not raise
%      "unknown celltype" -- it quietly tries to build a BOOMERANG strip
%      instead, and whatever it then complains about points at the wrong file.
%      It must therefore stay ABOVE that bare else.
%
%   CreateFileBase.m does NOT need touching: solveBands builds fBase with its
%   own inline chain and never calls it.
%
% The guard below still checks for both branches, so it catches a regression;
% as things stand it passes and a solve is attempted.
%
% CONSIDER RUNNING WITH 1 FIRST, NOW THAT THE FOOTPRINT HAS CHANGED.
% P.cutBottomHalfCell moves the bottom face into the middle of cell 1, which
% (a) makes P.yEnd1 two boundaries instead of one, (b) makes P.r2 actually
% round two corners for the first time, and (c) exposes a 23.5 nm x 125 nm
% ligament face (33.5 nm before that fillet) that the mesh has to resolve. All
% three are worth one look in the GUI before paying for a 70-band, 6-k-point
% solve -- and P.yEnd1 having two entries is exactly what the print block below
% reports.
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

    % Expect two warnings at these parameters -- :filletSelectionBleed (the
    % P.cutBottomHalfCell cut vertices land in disksel2) and :filletR2NoOp (the
    % arm tips do not). Both are the reproduced legacy fillet selections
    % reporting themselves; the first one does change the geometry, by rounding
    % the two cut corners with r2. See the P.r1/P.r2 comments above.
    fprintf('\nbuildCrossStrip: %d cells, Ly = %.1f nm\n', P.ncell, P.Ly*1e9);
    fprintf('  cell centres [nm] : %s\n', num2str(P.ycentres*1e9, '%.1f  '));
    fprintf('  xEnd1 / xEnd2     : [%s] / [%s]\n', num2str(P.xEnd1), num2str(P.xEnd2));
    fprintf('  yLo / yHi    [nm] : %.1f / %.1f\n', P.yLo*1e9, P.yHi*1e9);
    fprintf('  yEnd1 / yEnd2     : [%s] / [%s]\n', num2str(P.yEnd1), num2str(P.yEnd2));
    if P.cutBottomHalfCell
        fprintf(['  bottom face is CUT through cell 1, so yEnd1 should hold ' ...
                 'TWO boundaries,\n  each %.1f nm wide in x (%.1f nm after ' ...
                 'the r2 fillet).\n'], (P.a-P.h)/2*1e9, ((P.a-P.h)/2-P.r2)*1e9);
    end
    % The fixed BC rides on this NAME, not on yEnd2. If yFixedInds is empty the
    % builder has already warned (:yFixedSelEmpty) and runBands would error
    % (:fixedSelEmpty) rather than solve an unconstrained model.
    fprintf('  bndSel.yFixed     : %s -> [%s]\n', ...
        P.bndSel.yFixed, num2str(P.bndSel.yFixedInds));
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
