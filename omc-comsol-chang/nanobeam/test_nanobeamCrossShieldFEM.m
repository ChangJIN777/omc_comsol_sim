% Front panel script for a nanobeam OMC cavity terminated by a 2D cross-cell
% phononic shield, with a wrap-around mechanical PML for the radiative Q, and
% - when P.solveOpt = 1 - with the optical mode solved in the SAME COMSOL
% model so that the optomechanical coupling g_OM can be computed.
%
% THIS SCRIPT RUNS THE FULL-Y (ASYMMETRIC) GEOMETRY. P.meveny = 0, so there is
% no y = 0 mirror plane: the whole beam width is simulated, the shield sits
% ENTIRELY ON THE +y SIDE with its bottom edge FLUSH with the beam's bottom at
% y = -P.w/2, and the PML follows. In nm, from the parameters below:
%
%     y = -400.0 .......... beam bottom == shield bottom, flush
%     y in [-400, +400] ... beam, for x in [0, 11063.0]
%     y in [-400, 4070.0] . cross shield, for x in [11063.0, 16427.0]
%                           = -400 + P.nShieldY*P.aShield = -400 + 5*894
%     y in [4070, 6070.0] . +y PML arm (+ P.PMLLenY = 2000)
%     nothing below y = -400: bare vacuum for every x
%
% WHY THERE IS NO -y PML AND NO -y SHIELD. A traction-free surface facing
% vacuum carries no elastic energy away - vacuum supports no elastic waves, so
% only material radiates. Below y = -400 nm there is no material at any x, so
% there is nothing to absorb. It is the same argument that already excludes the
% +/-z faces and the beam's own +/-y flanks at x < xS0 (those flanks were free
% vacuum in the old half-y build too, on the +y side). The PML therefore stays
% the SAME three blocks - +x arm, +y arm, corner - and the builder still checks
% that they resolve to three separate domains (PMLDomainCount).
%
% WHAT IT COSTS, HONESTLY. Read before believing a Q:
%   - the clamp is worse. The beam end face at x = xS0 now spans
%     y in [-400, +400] and 416.0 nm of that - 52.00%, against 40.25% in the old
%     mirror-symmetric build - lands on a 15.5 nm free-standing wall fin
%     instead of on a node. The builder warns (beamOverhangsNode) with the real
%     number. P.shieldPadLen > 0 (about aShield/2 = 447 nm) is the fix and is
%     deliberately NOT applied here - see the P.shieldPadLen line below.
%   - the bottom edge is halved. Cutting the lattice on the y = -400 nm
%     cell-corner line bisects the bottom node row and the bottom row of
%     x-running ligaments, leaving 6 free half-ligaments of 15.5 nm right
%     alongside the clamp. Expect edge-localised modes running along that edge
%     towards the +x PML (warning: flushBottomHalvedRow).
%   - the device is different. The old build's physical shield was 8940 nm tall
%     and mirror-symmetric about the beam; this one is 4470 nm and one sided.
%     Half the shield material, and the beam is at its weakest edge rather than
%     at its centre. Do not compare Q against a half-y run and call it a
%     convergence test.
%   - the mode count doubles. With no y parity filter both y-even and y-odd
%     modes come back from one solve, so P.mneigs and P.oneigs may need raising
%     - see the notes at each.
% What it does NOT cost is shield or PML mesh: the shield is still
% P.nShieldY*P.aShield = 4470 nm tall and the +x PML arm still spans that same
% height, so only the BEAM volume doubles. See the DOF BUDGET note.
%
% P.solveOpt (MASTER SWITCH block below, immediately after the clear) is the
% single on/off for the optical half of the run. With it at 0 this is a
% mechanics-only Q simulation and nothing in the optical sections applies;
% results land in a separate datLoc so the two states cannot reuse each
% other's cached solve.
%
% Pipeline:  RunNanobeamFEM -> CreateNanobeamGeom
%                           -> BuildNanobeamCrossShieldFEM   (P.celltype)
%                           -> SetupNanobeamFEM              (reused unchanged)
%                           -> SolveNanobeamFEM              (reused unchanged)
%                           -> CalcGOM                       (needs BOTH studies
%                                                             on one model)
%
% WHY ONE MODEL: CalcGOM.m:25-28 requires ds.ofem AND ds.mfem, and
% CalcGOM.m:124-155 Joins the optical and mechanical datasets, so the two modes
% have to come out of a single geometry. That is the whole reason the builder
% had to learn optics instead of the optical run being a separate script.
%
% THE AIR CYLINDER IS SCOPED TO THE CAVITY (P.solveOpt = 1 only; with the
% switch off no air domain is created at all). The full-length air cylinder of
% BuildNanobeamFEM would enclose the shield and the PML frame and destroy the
% phononic bandgap. Here it stops at P.airCylLen, one free-space wavelength
% short of the shield front face, with 15 of the 18 half-beam holes inside it.
% Read WHY THE AIR CYLINDER STOPS SHORT in BuildNanobeamCrossShieldFEM.m before
% changing P.airrad, P.airCylLen, P.nholes or P.shieldPadLen - the builder
% hard-errors on an overlap, but the convergence question (stage S4b) is yours.
%
% BRING-UP ORDER - do not skip straight to a solve. See the staged block at the
% bottom of this file; S0-S3 need no solver and S1 needs only the geometry.

clear all; close all; clc                                   %#ok<CLALL>

%% MASTER SWITCH - optics on (1) or off (0)
% P.solveOpt is the ONLY toggle for the optical side of this run. Everything
% else is derived from it - P.calcG and P.plotOpt further down, the datLoc
% suffix at the bottom, and every optical addition in
% BuildNanobeamCrossShieldFEM.m (which gates on useOpt at :382). Do NOT add a
% second flag: two sources of truth can disagree, and the builder would then
% have to decide which one it believes.
%
% What = 1 buys, and costs, relative to a mechanics-only run:
%   + the cavity-scoped air cylinder (P.airrad / P.airCylLen below) plus its
%     end-cap imprint, which splits the beam into two domains. With
%     P.meveny = 0 that is a HALF cylinder (all y, z >= 0); it was a quarter
%     (y >= 0, z >= 0) in the half-y build, so the optical mesh is twice the
%     size it was - still comfortable, see the DOF BUDGET note;
%   + the emw physics and its own hmax = lambda/5 mesh pass;
%   + the optical eigenfrequency study and the g_OM calculation in CalcGOM.
% It costs no extra MECHANICAL degrees of freedom - emw and solid are selected
% on disjoint domain sets - so the DOF budget note at P.max_dof binds either
% way. See "WHY ONE MODEL" above for why optics has to live in this model.
%
% What = 0 needs: nothing. The builder's optical field checks, the air
% cylinder, and BOTH consistency guards (oevenz vs mevenz, oeveny vs meveny)
% all sit inside the useOpt branch, so P.airrad / P.airCylLen / P.oeven* below
% may stay set and unused - and cannot fire a spurious error either.
P.solveOpt = 1;

%% geometry parameters - beam and holes
P.xsect = 'rect';                       % v1 of the cross-shield builder is rect only
P.celltype = 'crossShield';             % dispatch key in RunNanobeamFEM
P.beamMat = 'diamond';
P.anisoMat = 1;

P.a = 650e-9;                           % nominal HOLE lattice constant
P.w = 800e-9;                           % beam width
P.theta = 45;                           % etch angle (no effect for rect)
P.th = 250e-9;                          % slab thickness - beam AND shield
P.hx = 343e-9;                          % nominal hole height (along x)
P.hy = 617e-9;                          % nominal hole width  (along y)

P.nholes = 18;                          % # holes in 1/2 beam length
P.ndef = 6;                             % # holes in 1/2 defect region
P.maxdef = 0.16;                      % defect percentage
P.oblong = 1.15;                       % oblong parameter
P.aspect_ratio = 1;

% cavity taper params
P.holeatctr = 1;                        % 1/0 for hole/dielectric in middle
P.taperFunc = 'cubic';

% end waveguide mirror taper params
P.wvgmir = 0;                           % the cross shield is the termination
P.wgmTaper.func = 'cubic';
P.wgmTaper.endtype = 'custom';
P.wgmTaper.a_end = 650e-9;
P.wgmTaper.hx_end = 343e-9;
P.wgmTaper.hy_end = 300e-9;

P.asymCav = 0;                          % not supported by this builder in v1

P.lambda = 1550e-9;
P.nbeam = 2.386;

% Disorder
P.stdDev = [0,0];
P.stdDevPos = 0;
P.asym = 0;
P.extractLocMechModes = 1;

%% cross-cell phononic shield
% ---------------------------------------------------------------------------
% WHAT THESE NUMBERS ACTUALLY ARE. The prose here used to describe a different
% cell (a = 1600, h = 1500, w = 800 nm against a 750 nm beam) and none of its
% conclusions survive the parameters that are really set below. Recomputed:
%   pitch  = aShield = 894 nm. This IS the right pitch for 7 GHz in diamond:
%            v_shear/(2f) = sqrt(c44/rho)/(2f) = sqrt(578e9/3500)/14e9
%            = 918 nm. The old "1.6 um is not that" tension is gone - these
%            numbers are the tight-pitch horn of that dilemma, already chosen.
%   node   = aShield - wShield = 478 nm, against a beam width of P.w = 800 nm.
%            The beam does NOT land entirely on a node, in either build, and
%            the clamp audit does NOT stay quiet at P.shieldPadLen = 0. In the
%            full-y build 416.0 nm of the 800 nm end face - 52.00% - lands on a
%            15.5 nm fin; in the old half-y build it was 40.25%. Expect
%            beamOverhangsNode on every run and read the number it prints.
%   wall   = aShield - hShield = 31 nm. THIS IS BELOW P.shieldMinFeature =
%            50e-9, so isCrossFabricable fails and the builder warns
%            (notFabricable). 31 nm is not a lithography error in the code - it
%            is what these three numbers mean - but it is not fabricable at the
%            stated limit either. LEAVE THE VALUES AS THEY ARE until the S10
%            bandgap check says whether the gap survives widening the wall;
%            changing aShield/hShield without re-running S10 trades a warning
%            for a shield that is outside its own bandgap.
%   solid fraction = 1 - (2*h*w - w^2)/a^2 = 0.318, so there is real mass
%            contrast.
%   shieldHmax default = min((a-h)/3, th/3) = min(10.3, 83.3) = 10.3 nm, NOT
%            the 33.3 nm an older comment claimed. That 10.3 nm is what makes
%            the mesh estimate in the DOF BUDGET block enormous.
%
% Before believing any Q, run stage S10: feed the SAME aShield/hShield/wShield
% and th = P.th through the band structure pipeline (buildCrossUnitCell +
% runBands) and confirm that P.freq actually lands inside the cross lattice's
% bandgap. A shield outside its gap gives a Q that does not grow with
% nShieldX, which is the tell.
% ---------------------------------------------------------------------------
P.aShield = 894e-9;                    % cross lattice pitch, BOTH x and y
P.hShield = 863e-9;                    % cross arm length (must be < aShield)
P.wShield =  416e-9;                    % cross arm width

% Wider-wall alternative if S10 says the gap tolerates it - expect a smaller
% beamOverhangsNode fraction and a 10x cheaper shield mesh:
% P.aShield = 900e-9;  P.hShield = 800e-9;  P.wShield = 400e-9;

P.nShieldX = 6;                         % cells along x
P.nShieldY = 5;                         % cells along y. WITH P.meveny = 0 THIS
                                        % IS THE ABSOLUTE ROW COUNT, not a
                                        % half-model count: 5 physical rows,
                                        % all of them on the +y side of the
                                        % beam, running up from the flush
                                        % bottom at y = -P.w/2 = -400 nm to
                                        % y = -400 + 5*894 = 4070 nm.
                                        % The filename records it as nsyF5.
                                        % (With abs(P.meveny) = 1 the same 5
                                        % would mean 10 physical rows mirrored
                                        % about y = 0, written nsyH5_nsyP10 -
                                        % which is why the two builds cannot
                                        % collide in datLoc.)

P.shieldPadLen = 0;                     % 0 = lattice faithful (default).
                                        % ~aShield/2 = 447e-9 gives a proper
                                        % clamp, bonds away the x = xS0
                                        % free-edge fins, and is the documented
                                        % fix for the 52.00% fin overlap the
                                        % flush bottom causes. It is left at 0
                                        % ON PURPOSE: turning it on changes the
                                        % device, not just the mesh, so that
                                        % is your call and not the builder's.
P.shieldFillet = 'none';                % 'none' | 'analytic' | 'legacy'.
P.r1Shield = 0;                         % ignored while shieldFillet = 'none'
P.r2Shield = 0;
P.shieldMinFeature = 50e-9;             % for the isCrossFabricable audit
% P.shieldHmax = 40e-9;                 % default is min((a-h)/3, th/3)

%% simulation / calculation / plot / save options
P.solveMech = 1;
% P.solveOpt is set ONCE, in the MASTER SWITCH block at the top of the file.
% The two lines below derive from it and must not be hand-edited.
P.calcG = 1*(P.solveMech && P.solveOpt);    % g_OM needs both modes, so this is
                                        % the AND of the two solve flags
P.calcS = 0*P.solveMech;
P.solveMechPML = 1;                     % 1 to extract the radiative mechanical Q

P.plotgeom = 1;                         % mphgeom snapshot into datLoc
P.storeMPH = 1;                         % keep the .mph - you will want to look
P.plotMech = 1*P.solveMech;
P.plotOpt = 1*P.solveOpt;               % follows the master switch
P.plotStrCpl = 1*P.calcS;

%% mechanical simulation parameters
P.mevenx = 1;                           % +/-1 even/odd about x
% P.meveny = 0 IS THE ASYMMETRIC BUILD, and it is the point of this script.
% 0  -> no y = 0 mirror plane. The full beam width is simulated, the shield
%       bottom is flush with the beam bottom at y = -P.w/2, and P.nShieldY
%       counts absolute rows on the +y side only. The builder emits no
%       P.bndSel.beamYsym / PMLYsym / cylYsym at all, which is safe because
%       SetupNanobeamFEM reads them only from inside its meveny/oeveny == +/-1
%       branches (:245, :271, :370, :412) and CalcGOM.m:92 is guarded the same
%       way. CalcGOM.m:48 drops the y factor from symFac automatically.
% +/-1 -> the old half-y build: y = 0 is a symmetry (+1) or antisymmetry (-1)
%       plane, P.nShieldY counts cells in the simulated half, and the device
%       has 2*P.nShieldY mirrored rows. Still fully supported and byte for byte
%       what it was; flip this back and the geometry, selections and filename
%       are identical to before the full-y path existed.
% P.oeveny below MUST agree: 0 here forces 0 there, and the builder hard-errors
% (evenyMismatch) if they disagree while optics is on.
P.meveny = 0;
P.mevenz = 1;                           % +/-1 even/odd about z
P.freq = 7e9;                           % target mechanical frequency
P.mneigs = 30;                          % more than the no-PML case: a PML run
                                        % also returns PML-localised modes.
                                        % WITH P.meveny = 0 the y parity filter
                                        % is gone as well, so both y-even and
                                        % y-odd families come back from one
                                        % solve and the useful modes are
                                        % sparser in the list. If S5 does not
                                        % show a cavity mode near P.freq
                                        % surviving the locRatio filter, raise
                                        % this to 50-60 before suspecting the
                                        % geometry.
P.mMesh = 3;
P.mAdjMesh = 1;

P.rxtal = 0;
P.rxtalInFilename = 1;

%% optical simulation parameters - READ ONLY WHEN P.solveOpt = 1
% Everything in this section is inert with the master switch off: the builder
% reads no optical field outside its useOpt branch, so these may stay set.
%
% Symmetries follow the rect-beam fundamental of test_nanobeamRectFEM.m:99-102,
% the closest sibling, EXCEPT in y where this script has no mirror plane:
%   oevenx = (-1)^holeatctr = -1 for a hole at the centre. SetupNanobeamFEM.m:39
%           then raises it to the power holeatctr again, so the value that
%           reaches the BC is (-1)^1 = -1 -> PEC on the x = 0 plane.
%   oeveny = 0   NOT -1. There is no y = 0 plane in this geometry to hang a
%           PMC or a PEC on, so the only consistent value is 0. WHILE OPTICS IS
%           ON this must match P.meveny: there is one y geometry and two
%           physics, and the builder hard-errors (evenyMismatch) if they
%           disagree. Setting -1 here with P.meveny = 0 would ask
%           SetupNanobeamFEM.m:414 for P.bndSel.cylYsym, a field the full-y
%           builder does not set.
%           WHAT YOU GIVE UP: with oeveny = -1 the solver could only return
%           TE-like modes (E mostly along y, E_y odd about y = 0). With no y
%           BC at all BOTH parities are in the spectrum, so the shift-invert
%           target at P.lambda may land on the y-even mode instead of the
%           TE-like fundamental. That is what P.oneigs below is about.
%   oevenz = +1  even about the slab mid-plane -> PMC on z = 0. WHILE OPTICS IS
%           ON this must have the same |value| as P.mevenz (both 1 here):
%           there is one z = 0 cut in the geometry and the builder
%           hard-errors (evenzMismatch) if the two physics disagree about
%           whether to make it. With P.solveOpt = 0 neither guard is reached
%           and P.mevenz/P.meveny alone decide the geometry, so leftover
%           P.oeven* values cannot fire a spurious error.
P.oevenx = (-1)^(P.holeatctr);
P.oeveny = 0;                           % must be 0 while P.meveny = 0
P.oevenz = 1;

% oneigs = 1 matches every rect sibling (test_nanobeamRectFEM.m:102,
% test_nanobeamRectFEM_withPML.m:98, test_nanobeamTri*.m, and
% test_fabricatednanobeamRectFEM.m:137); only test_1DSnowflakeFEM.m:109 asks
% for 5. Keep 1 here: the shift-invert target is set to the exact design
% wavelength (SolveNanobeamFEM.m:64 shifts to -1i*2*pi*c/P.lambda), and
% RunNanobeamFEM.m:163 passes oModes = 1 to CalcGOM unconditionally - so
% asking for more eigenvalues would not change which one the coupling is
% computed for. If S4b shows the solver landing on a band-edge mode instead of
% the fundamental, raise this to 5 and re-check which solution index is the
% high-Q one (SolveNanobeamFEM.m:199 prints ofem.hiQSol).
%
% THE CASE FOR RAISING IT HERE IS STRONGER THAN IN THE SIBLINGS. Those scripts
% all set oeveny = -1, so the PEC on y = 0 removes the y-even family from the
% spectrum before the solver sees it and the nearest eigenvalue to P.lambda is
% the TE-like fundamental almost by construction. With P.oeveny = 0 there is no
% such filter and both parities are candidates. Run S4 at oneigs = 1 first; if
% PlotEy shows an E_y that is EVEN about y = 0 (a bright lobe on the axis
% rather than a node there), that is the wrong mode - raise this to 5 and take
% the y-odd one. Note RunNanobeamFEM.m:163 passes oModes = 1 to CalcGOM
% unconditionally, so if the fundamental is not solution 1 the g_OM will be
% computed for the wrong optical mode and will look plausible.
P.oneigs = 1;
P.oMesh = 3;
P.oAdjMesh = 1;

% Air cylinder around the CAVITY. 2*lambda + w/2 is the convention in
% test_nanobeamRectFEM.m:105 and five sibling scripts.
%
% WITH P.meveny = 0 THIS IS A HALF CYLINDER, NOT A QUARTER: the builder's
% Revolve sweeps the full 360 degrees (angle1 = -180, angle2 = +180) so the air
% covers ALL y, and the z = 0 cut then removes z < 0, leaving y in
% [-3500, +3500], z in [0, 3500]. In the old half-y build the Revolve swept 180
% degrees (-180 -> 0) and the result was the quarter cylinder y >= 0, z >= 0.
% Twice the air volume, so twice the optical tet count - see the DOF BUDGET.
P.airrad = 2*P.lambda + P.w/2;          % 2*1550 + 400 = 3500 nm
%
% P.airCylLen is the +x extent of that cylinder and is the field that makes
% optics compatible with the shield at all. Leave it unset to take the builder
% default, min(beamLenHalf, xS0 - P.lambda). At the parameters above, verified
% by running LoadMaterialParams + CreateNanobeamGeom headless on R2025b
% (stage S0b - do this again after any change to P.nholes, P.ndef, P.maxdef or
% P.taperFunc, because every number here moves with them):
%
%   beamLenHalf = 11063.0 nm = xpos(18) + a_hole(18)/2 = 10738.0 + 325.0
%   shieldPadLen = 0   ->  xS0 = 11063.0 nm
%   airCylLen   = 11063.0 - 1550 = 9513.0 nm
%   x-gap to the shield front face   = 1550.0 nm  (1.00 lambda)
%   x-gap to the +x PML arm          = 6914.0 nm
%   x-gap to the +y PML arm          = 1550.0 nm, y-gap = 570.0 nm
%     (the y-gap is 4070.0 - 3500.0 now, not 970.0: the +y PML arm starts at
%      yS1 = -400 + 5*894 = 4070 nm, 400 nm lower than the half-y yS1 = 4470)
%   15 of 18 half-beam holes - 6 defect + 9 mirror periods - end before the
%   cap, so the mirror the optical mode sees is a real one.
%
% BUT the default cap at 9513.0 nm lands INSIDE hole 16 (x in [9266.5, 9609.5]),
% so the builder will raise airCylCutsHole: that hole's cross-section sits
% inside the beam bounding box and gets differenced out of P.bndSel.cylXend,
% leaving part of the cap with no scattering BC. Uncomment the override below
% to put the cap in the solid between hole 15 (ends 8959.5) and hole 16, which
% moves the gap to the shield to 1950 nm = 1.26 lambda.
%
% The cylinder reaches |z| = 3500 nm against a 250 nm slab and y = +3500 nm
% against a shield whose top is at 4070 nm, so x is the ONLY separation - in
% BOTH builds, and it is what the builder's hard error guards. On the -y side
% the cylinder now hangs 3100 nm below the flush bottom at y = -400 nm, which
% is harmless: there is no material down there to be enclosed. Shorten the beam
% (P.nholes) or lengthen the pad (P.shieldPadLen) and the x gap moves with it.
% P.airCylLen = 9.113e-6;               % uncomment: mid-solid, hole 15|16

%% SiV strain coupling (unused while calcS = 0)
P.zSiV = {[1 1 1]};
P.xSlc = 0;
P.ySlc = 0;
P.zSlc = 0;

%% mechanical PML
% lambda_shear = sqrt(c44/rho)/f = sqrt(578e9/3500)/7e9 = 1.84 um for diamond
% at 7 GHz, so 2 um is about one wavelength. Confirm with the S8 sweep.
P.PMLLen  = 2e-6;                       % +x frame thickness
P.PMLLenY = 2e-6;                       % +y frame thickness
P.PMLmeshDiv = 8;                       % elements across the absorbing direction
P.PMLScalingType = 'userDefined';       % FREQUENCY DEPENDENT: the PML stretch
                                        % is keyed to P.freq via
                                        % typicalWavelength = v/P.freq.
                                        % Set 'rational' for a wavelength
                                        % independent cross-check in S8.
% Reference speed for typicalWavelength. Unset = c11 (longitudinal, 17.5 km/s,
% lambda = 2.50 um at 7 GHz), which is the legacy choice. Uncomment for the
% shear branch (12.85 km/s, 1.84 um), which is what actually carries radiation
% out of a thin suspended slab.
% P.PMLWaveSpeed = sqrt(P.D(10)/P.rho);

%% simulation settings
% DOF BUDGET FOR THE COMBINED RUN - READ BEFORE THE FIRST SOLVE.
% SolveNanobeamFEM makes ONE mesh node and both studies share it, but the two
% physics are selected on different domains (emw on beam + shield + air, smech
% on beam + shield + PML), so they do NOT share degrees of freedom and adding
% optics does NOT inflate the mechanical count.
%
% HOW MUCH DOES P.meveny = 0 COST? Much less than doubling, because the term
% that dominates is fullY-INVARIANT. The shield is P.nShieldY*P.aShield = 4470
% nm tall in BOTH builds - full-y shifts it down by w/2, it does not resize it
% - and so is the +x PML arm that spans it. Only the BEAM volume doubles, and
% the beam is meshed at hauto, not at P.shieldHmax. Numbers below.
%
% Estimates at these parameters, recomputed headless on R2025b (the builder
% prints its own version of the same line; both use the crude V_tet = h^3/6,
% which understates a real tet fill by roughly 1.4x):
%   optical     hmax = lambda/5 = 310 nm over 183.0 um^3 of HALF cylinder
%               -> 3.69e4 tets, ~7.4e5 DOF.  Comfortable. (It was 91.5 um^3 /
%               1.84e4 tets / 3.7e5 DOF as a quarter cylinder in the half-y
%               build, so this is the one number that really does double.)
%   mechanical  shield hmax = min((a-h)/3, th/3) = min(10.3, 83.3) = 10.3 nm
%               - NOT the 33.3 nm an older copy of this note claimed; the wall
%               a-h is 31 nm at these cross parameters, not 100 nm - over
%               0.954 um^3 of solid shield -> 5.19e6 tets, plus 1.1e3 in the
%               PML at PMLLen/PMLmeshDiv = 250 nm -> ~2.2e7 DOF, and nearer
%               3.0e7 once the h^3/6 optimism is removed.
%               THAT IS ALREADY 4x OVER P.max_dof BELOW, and it was over
%               before P.meveny changed and is over with P.solveOpt = 0. The
%               beam contributes about 3.5e4 more tets in full-y than in half-y
%               (0.754 vs 0.377 um^3 at h ~ 40 nm), i.e. under 1% of the
%               mechanical total. The shield mesh is the problem, not the
%               symmetry.
%
% If the mechanical side runs out of budget, mAdjMesh cannot rescue it: the
% loop raises hauto, but applyMeshOverrides re-applies hmax = P.shieldHmax to
% the shield on every pass, so the dominant term never moves. Levers, in order:
%   1. P.nShieldX = 3, P.nShieldY = 3 for the combined run (the shield tet
%      count scales with nShieldX*nShieldY), then take Q from a mech-only run
%      at the full 6x5 and g_OM from the small one - the optical mode does not
%      touch the shield. NOTE that P.nShieldY = 3 in full-y gives a shield only
%      3*894 = 2682 nm tall from a flush bottom, which still clears the 800 nm
%      beam but leaves very little lattice above it.
%   2. P.shieldHmax = 14e-9 (just over two elements across the 31 nm ligament)
%      cuts the shield tets by 2.5x. Do not go to exactly (a-h)/2 = 15.5e-9:
%      the builder's shieldMeshCoarse warning fires there, and two elements
%      across the feature that carries the entire shield response is already
%      the floor. Even at 14 nm the estimate is ~2.1e6 tets / 8.7e6 DOF, so
%      this lever ALONE does not get under 5e6 - combine it with 1.
%   3. Widen the wall. a-h = 31 nm is the root cause of a 10.3 nm hmax and is
%      also below P.shieldMinFeature, so the notFabricable warning and the mesh
%      cost have the same fix. This changes the lattice, so it needs the S10
%      bandgap re-check first.
%   4. Raise P.max_dof if the machine has the memory - but 3e7 DOF is not a
%      memory setting away, it is a different mesh.
% A fifth fix - separate optical and mechanical mesh nodes - would need an
% edit to SolveNanobeamFEM, which this script deliberately does not require.
P.max_dof = 5e6;

%% single run
% THE OPTICS TAG IN datLoc IS WHAT MAKES P.solveOpt SAFE TO FLIP.
% CreateFileBase.m encodes no optical parameter at all, so a mechanics-only run
% and a combined optical+mechanical run produce the SAME P.fileBase.
% RunNanobeamFEM.m:64 skips the entire solve whenever a matching .mat AND .mph
% already sit in datLoc, takes its elseif branch instead, and there does
% P = ds.P (RunNanobeamFEM.m:157) - which overwrites the caller's P, P.calcG
% included. Run with optics off, flip the master switch on, re-run into the
% same directory and you would silently get the old mechanics-only result:
% no gOM, no error, nothing in the log to say so.
%
% Giving each toggle state its own directory makes that collision impossible
% and keeps the fix inside this script - CreateFileBase.m is shared with the
% band-structure pipeline and every other test script, so it is not the place
% to encode a per-script flag.
%
% THE Y CONVENTION IS A DIFFERENT CASE AND IS HANDLED IN CreateFileBase.m.
% P.solveOpt is a per-script run mode; P.meveny changes the DEVICE, and
% P.nShieldY = 5 means a 5-row one-sided shield here and a 10-row mirrored one
% with abs(P.meveny) = 1. That is a geometry parameter like P.aShield, so it
% belongs in the filename, and CreateFileBase.m now writes nsyF5 for the full-y
% build where it writes nsyH5_nsyP10 for the half-y one. The half-y string is
% unchanged byte for byte, so nothing already cached under test/ is orphaned.
if P.solveOpt
    optTag = 'optMech';                 % emw + solid, gOM computed
else
    optTag = 'mechOnly';                % solid only, no air domain
end
currentDate = datestr(now,'mmddyyyy');                      %#ok<TNOW1,DATST>
datLoc = [fullfile('.','test','1D_OMC_crossShield',currentDate,optTag),filesep];
[ds,model] = RunNanobeamFEM(P,datLoc);

%% ========================================================================
%  STAGED BRING-UP
%  ========================================================================
%
% S0  static analysis, no COMSOL:
%       checkcode('BuildNanobeamCrossShieldFEM.m')
%       checkcode('test_nanobeamCrossShieldFEM.m')
%       checkcode('CreateFileBase.m')
%     The first two should be silent. CreateFileBase has two pre-existing
%     warnings at its last line (the unused fBase) - unrelated, and unchanged.
%
% S0b geometry arithmetic, no COMSOL. P = LoadMaterialParams(P);
%     P = CreateNanobeamGeom(P); then check the numbers the air cylinder
%     depends on. These were re-derived on R2025b at the parameters above;
%     re-derive them again rather than trusting this comment if you change
%     P.nholes, P.ndef, P.maxdef or P.taperFunc:
%       P.beamLenHalf                       expect 11063.0 nm
%       2*P.lambda + P.w/2                   expect 3500.0 nm
%       P.beamLenHalf - P.lambda             expect 9513.0 nm  (= airCylLen)
%       max(P.geomHalf(:,2))/2 vs P.w/2      expect 308.5 vs 400.0 nm
%     The last one is load bearing TWICE OVER in the full-y build: the builder
%     separates the air's z = 0 faces from the beam's by a y bound, and needs
%     max(hy)/2 < w/2 so that the |y| <= max(hy)/2 hole box cannot also catch a
%     solid face. It also guarantees the holes do not break out through the
%     side wall, so each hole stays a hole rather than a notch.
%     Full-y extents worth writing down at the same time:
%       yS0 = -P.w/2                        expect -400.0 nm
%       yS1 = yS0 + P.nShieldY*P.aShield     expect 4070.0 nm
%       yP1 = yS1 + P.PMLLenY                expect 6070.0 nm
%
% S0c the filename, no COMSOL. P = CreateFileBase(P); disp(P.fileBase) and
%     check it ends ..._nsx6_nsyF5_pml2000nm. Then set P.meveny = 1, clear
%     P.fileBase, re-run and check it ends ..._nsx6_nsyH5_nsyP10_pml2000nm.
%     Two different strings is the whole point: RunNanobeamFEM.m:64 skips the
%     solve on a filename match, so if these collided a full-y run would be
%     handed a cached half-y result with nothing in the log to say so.
%
% S1  geometry only, MECHANICS ONLY. Comment out the RunNanobeamFEM call above,
%     set the MASTER SWITCH P.solveOpt = 0, P.nShieldX = P.nShieldY = 2 and
%     P.solveMechPML = 0, then run the block below. What to check in the GUI:
%       - the beam spans y in [-400, +400] nm, with WHOLE holes, not halves.
%       - the shield's bottom edge is at y = -400 nm, coincident with the
%         beam's own bottom face - flush, no step.
%       - nothing at all exists below y = -400 nm. If there is a mirrored
%         shield or a -y PML block down there, the wrong build path ran.
%       - the shield's bottom row is HALVED: the bottom node row is 239 nm deep
%         instead of 478, and the x-running ligaments along y = -400 are 15.5 nm
%         instead of 31. That is expected and the builder warns about it
%         (flushBottomHalvedRow); it is the geometric price of the flush bottom.
%       - the beam unions cleanly onto the shield at x = beamLenHalf.
%       - the console says "FULL-y (no y = 0 mirror) ... 2 physical rows,
%         bottom flush with the beam at y = -400.0 nm".
%       - NO selection warned about being empty. In particular you should NOT
%         see a warning about P.bndSel.beamYsym: the full-y builder does not
%         emit that field at all, and does not audit it either.
%
% S1a the half-y regression, and do it once. Set P.meveny = 1 (and, if optics
%     is on, P.oeveny = -1) and re-run the same block. Everything must look and
%     print exactly as it did before the full-y path existed: shield from y = 0
%     to y = nShieldY*aShield, half beam, "2x2 cross shield (4 physical rows)",
%     P.bndSel.beamYsym and PMLYsym present and non-empty. Every full-y change
%     in the builder is gated on abs(P.meveny) < 1, so this is a real
%     regression test and not a formality.
%
% S1b geometry only, OPTICS ON. Set the MASTER SWITCH back to P.solveOpt = 1
%     and repeat - note this writes into a DIFFERENT datLoc, by design, so S1
%     and S1b never collide. What to check in the GUI, in order of how badly
%     it bites:
%       - the air cylinder is a HALF cylinder (ALL y, z >= 0). Not a quarter -
%         that was the half-y build. If the y < 0 half is missing the Revolve
%         angles did not widen to -180/+180; if the z < 0 half is still there
%         the symZ cut block did not span far enough in z (or far enough in -y,
%         which is the new failure mode) and the z = 0 optical mirror plane is
%         wrong.
%       - there is bare vacuum, not air, between x = P.airCylLen and the
%         shield front face.
%       - the beam is TWO domains, split at x = P.airCylLen by the end-cap
%         imprint. That is expected; see the builder header.
%       - the console lines "Air cylinder clears the cross shield: x by
%         1550.0 nm", "Air cylinder -y reach: y = -3500.0 nm ... (bare vacuum,
%         nothing to overlap)", "Air end cap ... N interior (removed)",
%         "Air cylinder curved surface: N face(s) from the +y and -y 45 degree
%         probes" and "z = 0 plane: ... air (optical only) ... solid
%         (mechanical)". A 0 in any of those is a real failure, not cosmetics.
%       - the z = 0 line should now report MORE air faces than the half-y build
%         did: the outer air at z = 0 is two strips (+y and -y) instead of one,
%         plus one face per enclosed hole.
%       - disp(P.domSel) and disp(P.bndSel): cyl, cylXsym, cylXend, cylZsym,
%         cylCurv all present and non-empty. cylYsym is DELIBERATELY ABSENT -
%         there is no y = 0 plane. Its absence is safe because SetupNanobeamFEM
%         only reads it from inside its oeveny == +/-1 branches.
%
% S2  fillets. P.shieldFillet = 'analytic' with P.r1Shield = P.r2Shield = 10e-9
%     should build. 'legacy' is EXPECTED to warn and probably fail for
%     nx,ny > 1 - run it once so you have seen the failure mode.
%
% S3  mesh only. Re-enable the PML, then in the block below add
%       mesh = model.mesh.create('mesh', P.geomname);
%       mesh.feature('size').set('hauto', P.mMesh); mesh.run;
%       disp(mphxmeshinfo(model).nelements)
%     and tune P.shieldHmax / P.PMLmeshDiv until the count is affordable. See
%     the DOF BUDGET note next to P.max_dof; at 6x5 cells and the default
%     10.3 nm shield hmax this is expected to be roughly 4x OVER budget, with
%     or without optics and in either y build. Full-y is not what does it -
%     compare the element count against S1a at P.meveny = 1 and the difference
%     should be under a per cent, because only the beam doubles and the beam is
%     not what is being meshed at 10.3 nm.
%
% S4  OPTICS ONLY bring-up, and do this before the combined solve. Set
%     P.solveMech = 0... except the builder requires P.solveMech = 1, so
%     instead shrink the mechanics out of the way:
%       P.nShieldX = 1; P.nShieldY = 1; P.solveMechPML = 0; P.mneigs = 3;
%     A 1x1 shield still terminates the beam and still exercises the union,
%     but costs ~1.7e5 shield tets instead of ~5.2e6, so the optical study is
%     the only expensive thing in the model. NOTE that P.nShieldY = 1 in full-y
%     gives a shield only 894 nm tall from a flush bottom, which still clears
%     the 800 nm beam - but only just, and the beamOverhangsNode fraction will
%     be the same 52.00% because it is set by cell 1 alone.
%     Expect from SolveNanobeamFEM.m:210 a single high-Q mode with lambda
%     within a few tens of nm of 1550 and Q > 1e4; PlotEy should show a mode
%     centred on x = 0 and decaying into the mirror. If lambda is far from 1550
%     the taper is off, not the cylinder.
%     WITH P.oeveny = 0 ALSO CHECK THE PARITY, which the half-y build got for
%     free from the PEC on y = 0: PlotEy should show E_y ODD about y = 0, i.e.
%     a node on the beam axis with lobes of opposite sign either side. An E_y
%     with a single bright lobe on the axis is the y-even mode and is the wrong
%     one - raise P.oneigs to 5 and pick the y-odd solution.
%
% S4b optical convergence against the truncated air. Sweep
%     P.airCylLen = [7.5, 8.5, 9.113] um and watch lambda, Q and (once S7 runs)
%     g_OM. lambda and g_OM should be flat - they are set by the field at
%     x ~ 0. Q will NOT be flat and will be OVERESTIMATED at short airCylLen,
%     because the beam beyond the cylinder is clad in the EM interface's
%     default PEC, which is a perfect mirror where the device has vacuum. If
%     lambda or g_OM move by more than a per cent, the air is too short and
%     the cavity is being squeezed - the honest fix is fewer mirror holes
%     between the defect and the shield, not a longer cylinder.
%
% S5  eigenfrequency with P.solveMechPML = 0: real frequencies, cavity mode
%     near P.freq survives, shield modes fail the locRatio filter.
%     WITH P.meveny = 0 expect roughly TWICE as many modes in the same
%     frequency window, because both y-parity families are now in one spectrum.
%     Also expect a NEW family that the half-y build could not produce at all:
%     modes localised on the shield's free bottom edge at y = -400 nm, where
%     the halved 15.5 nm ligaments are. They should fail the locRatio filter
%     (SolveNanobeamFEM.m:395-407 weights |y| < P.w = 800 nm, so an edge mode
%     out at x > xS0 is not "in the centre of the beam"), but look at them
%     once: if the cavity mode HYBRIDISES with one, the flush-bottom design is
%     leaking along that edge and P.shieldPadLen > 0 or a taller shield is the
%     answer.
% S6  P.solveMechPML = 1: complex eigenvalues, real(lambda) > 0, finite QAll.
%     Do NOT compare this Q against a half-y run and call the difference a
%     convergence error - they are different devices. The half-y run's physical
%     shield is 8940 nm tall and mirror-symmetric about the beam; this one is
%     4470 nm and one sided, with the beam clamped at its weakest edge. A lower
%     Q here is the expected physics, not a bug.
% S7  combined solve with P.calcG = 1. CalcGOM Joins the optical and
%     mechanical datasets, so this is the first stage that can fail on the
%     shared mesh node rather than on physics: the optical study solves with
%     hmax = lambda/5 applied globally and no shield refinement, and the
%     mechanical study then re-meshes with hmax = P.shieldHmax on the shield.
%     Sanity checks on the result: ds.cpl.gMax should be of order 100 kHz for
%     a diamond OMC at 7 GHz, and ds.cpl.optWvl should equal the lambda
%     reported in S4. A g_OM that is orders of magnitude off, or a changed
%     optical wavelength, points at the mesh and not at the geometry.
% S8  Q vs P.PMLLen in {1, 2, 3, 4} um and Q vs P.PMLmeshDiv in {4, 8, 16}.
% S9  Q vs P.nShieldX = 1..5 - should climb roughly exponentially inside the
%     bandgap. Immediate saturation means the mode is not in the gap, or a BC
%     is wrong. This is the single most diagnostic test. g_OM must NOT change
%     with nShieldX; if it does, the optical mode is feeling the shield and
%     P.airCylLen is too long.
% S9b Q vs P.nShieldY = 2..6, which matters MORE in the full-y build than it
%     did in the half-y one. With a mirrored shield the beam sat in the middle
%     of 2*nShieldY rows and the y radiation path was nShieldY rows in each
%     direction. One sided, the +y path is still nShieldY rows but the -y
%     direction is a free edge 400 nm from the beam axis, so all the y
%     attenuation now has to come from the +y side. If Q saturates in nShieldX
%     but keeps climbing in nShieldY, the +y path is the limiter.
% S10 independent bandgap check via buildCrossUnitCell + runBands at th = P.th.
%     Do this before S9/S9b, not after: with a-h = 31 nm the cell is well
%     outside P.shieldMinFeature and the pitch was chosen for 7 GHz
%     (v_shear/(2f) = 918 nm vs aShield = 894 nm), so the gap is the one thing
%     here that has never been checked independently.
%
%% test the model building function (geometry only - no solve)
% import com.comsol.model.*
% import com.comsol.model.util.*
% 
% ModelUtil.showProgress(true);
% ModelUtil.clear();
% clear ds model
% ds = []; model = [];
% 
% model = ModelUtil.create('model');
% P = LoadMaterialParams(P);
% P = CreateNanobeamGeom(P);
% [model,P] = BuildNanobeamCrossShieldFEM(model,P);
% 
% disp('--- domain selections ---'); disp(P.domSel);
% disp('--- boundary selections ---'); disp(P.bndSel);
% disp('--- shield extents (m) ---'); disp(P.shield);
% if isfield(P,'airCyl'), disp('--- air cylinder (m) ---'); disp(P.airCyl); end
%
% P.shield gains two fields in the full-y build and only there: yS0 (the flush
% bottom, -P.w/2) and fullY (true). nYPhys is P.nShieldY in full-y and
% 2*P.nShieldY in half-y, so read it rather than assuming either.
%
% mphlaunch(model);
