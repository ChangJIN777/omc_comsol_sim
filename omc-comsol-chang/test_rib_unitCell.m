clear all; clc; close all;
clear P;

%% rib unit cell: mechanical band structure
% One period of the rib cavity (default: the visible design of
% comsol_templates/visible_rib_omc.m),
% built by buildRibUnitCell_LN.m and dispatched from runBands.m for
% P.celltype = 'rib'. A continuous spine (width s) carries one rib per
% cell (width t, reaching |y| = L). The rib is centred in the cell, so the
% Floquet faces at x = +-a/2 cut only the spine. Sibling scripts:
%   test_rib_unitCell_explicitW.m  explicit P.w, fillets, mechanical/optical
%   test_rib_unitCell_LN.m         lithium niobate, 15 deg sidewalls (old LN shape)
%
% Any field left unset takes the default from ribDefaultParams.m, which
% also documents where each number comes from. The values below repeat those
% defaults so they are visible here.
P.xsect = 'rect';
P.beamMat = 'diamond';                  % beam material name
P.celltype = 'rib';             % specify the cell type
P.a  = 300e-9;      % lattice constant = cell width (template mirror pitch)
P.t  = 100e-9;      % rib width along x
P.s  = 100e-9;      % spine width along y (full width, both sides)
P.th = 110e-9;      % slab thickness. 110 nm is what the template EFFECTIVELY
                    % simulates (z-offset defect, see ribDefaultParams);
                    % the nominal design value is 120e-9.
P.r1 = 40e-9;       % concave rib/spine fillet. TODO: not recoverable for the
P.r2 = 40e-9;       % convex rib-tip fillet.    visible set; IR-set values.
P.d  = 0;           % sidewall angle in RADIANS (e.g. 15*pi/180). With d ~= 0:
                    % r1 = r2 = 0 and mbevenz = 0 are required (TwoSymPlanes = 0,
                    % zSymCondition = 0).
% Rib half-length L (beam axis to tip). An explicit P.w (tip-to-tip) or P.L
% wins; otherwise the linear rule L(a) = L_ref + dLda*(a - a_ref) applies:
% P.w = 740e-9;     % e.g. fix the tip-to-tip width directly
P.L_ref = 400e-9;   % L at a = a_ref (template centre cell) ...
P.a_ref = 240e-9;   % ... at a = 240 nm
P.dLda  = -0.4998;  % slope of the template's IR set -> L = 370 nm at a = 300 nm
                    % Set P.dLda = 0 for a constant L = L_ref.
% To get the band structure of taper cell n of the cavity instead of the
% mirror cell, set P.a to that cell's a_n (createRibCavityGeom returns
% them as P.a_cell); L follows from the rule.

P.nperiod = 1;      % no. of periods to simulate for
P.holeatedge = 0;   % read by runBands, unused by this builder
P.run_optical = 0;  % MUST be 0 for runBands (an air domain would become diamond)

P.kpts = 10;                            % no. of k-points, EXCLUDING gamma point
P.nbands = 8;                          % no. of bands to solve for

P.TwoSymPlanes = 1;                     % solve all four y/z symmetry sectors
P.zSymCondition = 1;                    % used only when TwoSymPlanes = 0
P.solveasym = 1;                        % 1 to solve for antisymmetric bands
P.completeBandGaps = 1;                 % 1 to plot complete bandgaps (across all symmetries)
P.plotgeom = 1;                         % 1 to plot the geometry
P.savedat = 1;                          % 1 to save data structures
P.savebndplot = 1;                      % 1 to save bandstructure plot
P.saveplots = 0;                        % 1 to save displacement and strain profiles
P.saveMPH = 0;
P.bandStruct_2D = 0;                    % 0 to simulate 1D band structures

%% mechanical simulation parameters
P.mbeveny = 1;                          % overwritten per sector by solveBands
P.mbevenz = 1;                          % overwritten per sector by solveBands
P.freq = 0;                             % target frequency - 0 for band structures
P.meshSize = 4;                         % mesh quality for mechanical simulations
P.fixed_bc = 0;                         % no flat outer y face exists to clamp

P.anisoMat = 1;
P.rxtal = 0;                            % template's stiffness tensor is unrotated
                                        % (most other scripts use 45)
P.max_dof = 3e6;                        % max # of degrees of freedom

% %% debugging the unit cell geometry (no solve)
% import com.comsol.model.*
% import com.comsol.model.util.*
% ModelUtil.showProgress(true);
% ModelUtil.clear();
% clear model
% model = ModelUtil.create('model');
% [model,P] = buildRibUnitCell_LN(model,P);
% fprintf('L = %.1f nm, xEnd1 = [%s], xEnd2 = [%s], zEnd = [%s]\n', ...
%     P.L*1e9, num2str(P.xEnd1), num2str(P.xEnd2), num2str(P.zEnd));
% mphlaunch(model);

%% Single solve
currentDate = datestr(now,'mmddyyyy');
% test/ribUnitCell is shared with test_rib_unitCell_explicitW.m (the file names
% encode every dimension, so the two never collide). test/rib holds unrelated
% 2024 '2D_ribs' results.
P.datLoc = [fullfile('.','test','ribUnitCell',currentDate),filesep];
bds = solveBands(P);
