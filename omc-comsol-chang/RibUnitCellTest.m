classdef RibUnitCellTest < matlab.unittest.TestCase
%RIBUNITCELLTEST  Unit tests for the rib unit-cell builders and their helpers.
%
%   Covers
%     buildRibUnitCell.m      '2D_ribs': slab minus T-shaped voids
%     buildRibUnitCell_LN.m   'rib': the merged rib builder (spine + centred
%                             rib, fillets r1/r2, optional sidewall angle d),
%                             plus its helpers ribDefaultParams, ribLength,
%                             validateRibParams, addRibCell, addRibPrisms
%
%   Named RibUnitCellTest rather than test_*.m: in this repo test_*.m means a
%   simulation front-panel script, and test/ is the data output folder.
%
%   RUNNING (from omc-comsol-chang/, or anywhere: the class adds its folder)
%     results = runtests('RibUnitCellTest');                    % everything
%     results = runtests('RibUnitCellTest','Tag','PureMATLAB'); % no COMSOL
%     results = runtests('RibUnitCellTest','Tag','COMSOL');     % geometry
%     table(results)
%
%   TAGS
%     PureMATLAB  needs only base MATLAB (defaults, rib-length precedence,
%                 validation, file-base names, solver -> builder routing)
%     COMSOL      builds geometry in a fresh ModelUtil model (no mesh, no
%                 solve). Skipped (Incomplete, via an assumption) when LiveLink
%                 is not reachable: ModelUtil must already respond, or mphstart
%                 must connect to an ALREADY-RUNNING server. No server is ever
%                 started.
%     KnownBug    checks of behaviour the existing code gets wrong. They assert
%                 the CORRECT behaviour. While the bug is present the test ends
%                 Incomplete with a "KNOWN BUG (expected failure)" message
%                 (assumeFail), so it neither passes silently nor fails the
%                 suite. Once fixed it passes and logs that the note is stale.
%
%   FILE-BASE TESTS. solveBands / solveBands_noSym / solveOpticalBands build
%   P.fileBase inline and skip the solve when [datLoc fileBase '_bds.mat']
%   exists. The tests write a sentinel file under the EXPECTED name in a
%   temporary folder and call the solver: a cache hit proves the name. The
%   P structs omit the solver fields (TwoSymPlanes, bandStructureDim,
%   beamMat...), so a name MISMATCH errors immediately on a missing field
%   instead of starting a simulation.

    properties (TestParameter)
        % {mode, {field,value,...} overrides, expected error id}
        ValidationErrorCase = struct( ...
            'ribTooWideUnitCell',       {{'unitcell', {'t',300e-9},                  'rib:ribTooWide'}}, ...
            'ribTooWideCavity',         {{'cavity',   {'t',250e-9},                  'rib:ribTooWide'}}, ...
            'ribShorterThanSpine',      {{'unitcell', {'s',800e-9},                  'rib:ribShorterThanSpine'}}, ...
            'tipFilletTooLarge',        {{'unitcell', {'r2',50e-9},                  'rib:tipFilletTooLarge'}}, ...
            'rootFilletTooLargeCell',   {{'unitcell', {'r1',100e-9},                 'rib:rootFilletTooLarge'}}, ...
            'rootFilletTooLargeCavity', {{'cavity',   {'r1',70e-9},                  'rib:rootFilletTooLarge'}}, ...
            'badTaperConvention',       {{'cavity',   {'taperConvention','foo'},     'rib:badTaperConvention'}}, ...
            'badEndBC',                 {{'cavity',   {'mEndBC','pml'},              'rib:badEndBC'}}, ...
            'nonPositiveL',             {{'unitcell', {'dLda',-10},                  'rib:nonPositiveL'}}, ...
            'inconsistentLW',           {{'unitcell', {'L',300e-9,'w',700e-9},       'rib:inconsistentLW'}}, ...
            'airTooSmallCavity',        {{'cavity',   {'airrad',300e-9},             'rib:airTooSmall'}}, ...
            'airTooSmallUnitCell',      {{'unitcell', {'run_optical',1,'airradUC',200e-9}, 'rib:airTooSmall'}}, ...
            'badMode',                  {{'strip',    {},                            'rib:badMode'}}, ...
            'negativeLattice',          {{'unitcell', {'a',-300e-9},                 'MATLAB:validateRibParams:expectedPositive'}}, ...
            'nonIntegerCellCount',      {{'cavity',   {'ncellHalf',13.5},            'MATLAB:validateRibParams:expectedInteger'}}, ...
            'maxdefOutOfRange',         {{'cavity',   {'maxdef',1},                  'MATLAB:validateRibParams:notLess'}}, ...
            'sidewallInDegrees',        {{'unitcell', {'d',15},                      'rib:badSidewallAngle'}}, ...
            'sidewallNegative',         {{'unitcell', {'d',-0.1},                    'rib:badSidewallAngle'}}, ...
            'filletWithSidewall',       {{'unitcell', {'d',10*pi/180,'r1',20e-9},    'rib:filletWithSidewall'}}, ...
            'ribTopWidthNegative',      {{'unitcell', {'d',30*pi/180},               'rib:sidewallTooSteep'}}, ...
            'spineTopWidthNegative',    {{'unitcell', {'d',20*pi/180,'s',60e-9},     'rib:sidewallTooSteep'}}, ...
            'zSymWithSidewall',         {{'unitcell', {'d',5*pi/180,'mbevenz',1},    'rib:zSymWithSidewall'}});

        % {mode, overrides, expected warning id}
        ValidationWarningCase = struct( ...
            'filletsMeet',    {{'unitcell', {'L_ref',120e-9,'dLda',0}, 'rib:filletsMeet'}}, ...
            'stagingOverlap', {{'cavity',   {'dLda',25},               'rib:stagingOverlap'}});

        % {solver/runner file, celltype, expected builder}
        WiringCase = struct( ...
            'runBands_rib',          {{'runBands',          'rib', 'buildRibUnitCell_LN'}}, ...
            'runBandsNoSym_rib',     {{'runBands_noSym',    'rib', 'buildRibUnitCell_LN'}}, ...
            'runOpticalBand1D_rib',  {{'runOpticalBand_1D', 'rib', 'buildRibUnitCell_LN'}}, ...
            'runOpticalBand3D_rib',  {{'runOpticalBand_3D', 'rib', 'buildRibUnitCell_LN'}});

        SlabSymZ  = struct('noZsym', 0, 'zsym', 1);
        SymY      = struct('fullY', 0, 'halfY', 1);
        SymZ      = struct('fullZ', 0, 'halfZ', 1);
        Optical   = struct('noAir', 0, 'air', 1);
    end

    methods (TestClassSetup)
        function addThisFolderToPath(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fileparts(mfilename('fullpath'))));
        end
    end

    %% =====================================================================
    %  (a) PURE MATLAB
    %  =====================================================================
    methods (Test, TestTags = {'PureMATLAB'})

        % ---- ribDefaultParams --------------------------------------------
        function defaultsFillEveryDocumentedField(testCase)
            P = ribDefaultParams(struct());
            expected = {
                'celltype','rib'; 'xsect','rect'; 'beamMat','diamond';
                'anisoMat',1; 'rxtal',0; 'a',300e-9; 't',100e-9; 's',100e-9;
                'L_ref',400e-9; 'a_ref',240e-9; 'dLda',-0.4998; 'd',0;
                'th',110e-9; 'r1',40e-9; 'r2',40e-9; 'ncellHalf',14;
                'ribatctr',1; 'ndef',5; 'maxdef',0.2; 'taperFunc','quadratic';
                'taperConvention','centerspacing'; 'airrad',3.1e-6;
                'airZOffset',60e-9; 'mEndBC','lowreflecting'; 'lbd0',1.55e-6;
                'holeatedge',0; 'plotgeom',0; 'airradUC',1.5e-6};
            for k = 1:size(expected,1)
                testCase.verifyTrue(isfield(P, expected{k,1}), ...
                    sprintf('default field P.%s missing', expected{k,1}));
                testCase.verifyEqual(P.(expected{k,1}), expected{k,2}, ...
                    sprintf('wrong default for P.%s', expected{k,1}));
            end
        end

        function defaultsDoNotPrefillRibLength(testCase)
            % P.w / P.L must stay absent, or they would defeat the L(a) rule
            P = ribDefaultParams(struct());
            testCase.verifyFalse(isfield(P,'w'), 'P.w was pre-filled');
            testCase.verifyFalse(isfield(P,'L'), 'P.L was pre-filled');
        end

        function defaultsKeepCallerValues(testCase)
            in = struct('a',250e-9,'t',80e-9,'th',120e-9,'w',900e-9, ...
                'taperConvention','cellwidth','rxtal',45,'mEndBC','fixed', ...
                'unrelatedField','keep me');
            P = ribDefaultParams(in);
            names = fieldnames(in);
            for k = 1:numel(names)
                testCase.verifyEqual(P.(names{k}), in.(names{k}), ...
                    sprintf('caller value of P.%s was overwritten', names{k}));
            end
            testCase.verifyEqual(P.s, 100e-9, 'unset field not defaulted');
        end

        function defaultsAirradUCFollowsCallerLattice(testCase)
            P = ribDefaultParams(struct('a',250e-9));
            testCase.verifyEqual(P.airradUC, 5*250e-9, 'AbsTol', 1e-20);
            P = ribDefaultParams(struct('a',250e-9,'airradUC',2e-6));
            testCase.verifyEqual(P.airradUC, 2e-6, 'explicit airradUC overwritten');
        end

        function defaultsFilletsZeroWithSidewall(testCase)
            % LN-style call (d ~= 0, no fillet fields) -> the legacy LN shape
            P = ribDefaultParams(struct('d',15*pi/180));
            testCase.verifyEqual([P.r1 P.r2], [0 0]);
            P = ribDefaultParams(struct('d',15*pi/180,'r1',5e-9));
            testCase.verifyEqual(P.r1, 5e-9, 'explicit r1 overwritten');
        end

        function defaultsAcceptNoOrEmptyInput(testCase)
            testCase.verifyEqual(ribDefaultParams(), ribDefaultParams(struct()));
            testCase.verifyEqual(ribDefaultParams([]), ribDefaultParams(struct()));
        end

        % ---- ribLength: precedence P.L > P.w > L(a) rule ------------------
        function ribLengthAtAnchorIsLref(testCase)
            P = ribDefaultParams();
            testCase.verifyEqual(ribLength(P, P.a_ref), P.L_ref);
        end

        function ribLengthOfMirrorCell(testCase)
            P = ribDefaultParams();
            [L, src] = ribLength(P, 300e-9);
            testCase.verifyEqual(L, 370.012e-9, 'AbsTol', 1e-15);
            testCase.verifyEqual(src, 'rule');
        end

        function ribLengthIsElementwiseAndLinear(testCase)
            P = ribDefaultParams();
            a = (240:10:300)*1e-9;
            L = ribLength(P, a);
            testCase.verifySize(L, size(a));
            testCase.verifyEqual(diff(L), P.dLda*diff(a), 'AbsTol', 1e-18);
        end

        function ribLengthConstantForZeroSlope(testCase)
            P = ribDefaultParams(struct('dLda',0));
            testCase.verifyEqual(ribLength(P, [200 300 400]*1e-9), P.L_ref*[1 1 1]);
        end

        function ribLengthRejectsNonPositive(testCase)
            P = ribDefaultParams(struct('dLda',-10));
            testCase.verifyError(@() ribLength(P, 300e-9), 'rib:nonPositiveL');
        end

        function ribLengthExplicitWidthWins(testCase)
            P = ribDefaultParams(struct('w',1747e-9));
            [L, src] = ribLength(P, [240 300]*1e-9);
            testCase.verifyEqual(L, 1747e-9/2*[1 1], 'AbsTol', 1e-18);
            testCase.verifyEqual(src, 'w');
        end

        function ribLengthExplicitHalfLengthWins(testCase)
            P = ribDefaultParams(struct('L',350e-9));
            [L, src] = ribLength(P, 300e-9);
            testCase.verifyEqual(L, 350e-9);
            testCase.verifyEqual(src, 'L');
            P.w = 700e-9;                       % consistent with L: allowed
            testCase.verifyEqual(ribLength(P, 300e-9), 350e-9);
        end

        function ribLengthIgnoresBuilderWrittenValues(testCase)
            % what the builders write back (P.L, P.w + P.ribDerived) must not
            % count as explicit on a second call, e.g. after changing a
            P = ribDefaultParams(struct('a',300e-9));
            P.L = 370.012e-9;  P.w = 2*P.L;
            P.ribDerived.L = P.L;  P.ribDerived.w = P.w;
            [L, src] = ribLength(P, 260e-9);
            testCase.verifyEqual(src, 'rule');
            testCase.verifyEqual(L, P.L_ref + P.dLda*(260e-9 - P.a_ref), 'AbsTol', 1e-18);
            P.w = 800e-9;                       % a NEW value is explicit again
            testCase.verifyEqual(ribLength(P, 260e-9), 400e-9, 'AbsTol', 1e-18);
        end

        % ---- validateRibParams -------------------------------------------
        function validationDefaultsAreWarningFree(testCase)
            P = ribDefaultParams();
            testCase.verifyWarningFree(@() validateRibParams(P,'unitcell'));
            testCase.verifyWarningFree(@() validateRibParams(P,'cavity'));
        end

        function validationAcceptsLegacyLNParams(testCase)
            P = ribDefaultParams(lnRibParams());
            testCase.verifyWarningFree(@() validateRibParams(P,'unitcell'));
        end

        function validationErrors(testCase, ValidationErrorCase)
            [mode, over, id] = ValidationErrorCase{:};
            P = ribDefaultParams(struct(over{:}));
            testCase.verifyError(@() validateRibParams(P, mode), id);
        end

        function validationWarnings(testCase, ValidationWarningCase)
            [mode, over, id] = ValidationWarningCase{:};
            P = ribDefaultParams(struct(over{:}));
            testCase.verifyWarning(@() validateRibParams(P, mode), id);
        end

        function validationIgnoresUnitCellAirWithoutRunOptical(testCase)
            P = ribDefaultParams(struct('run_optical',0,'airradUC',1e-9));
            testCase.verifyWarningFree(@() validateRibParams(P,'unitcell'));
        end

        % ---- file-base names (sentinel cache hit, see class header) -------
        function solveBandsFileBaseRib(testCase)
            testCase.assumeTrue(ispc, ['solveBands appends a literal ''\'' to ' ...
                'P.datLoc, so its cache lookup only works on Windows.']);
            dl = tempDatLoc(testCase);
            P = struct('celltype','rib','datLoc',dl);
            name = 'rib_a_300nm_s_100nm_w_740nm_t_100nm_th_110nm_d_0deg_r1_40nm_r2_40nm_rxt_0_';
            writeSentinel(dl, name, struct('sym',[]));
            ds = solveBands(P);
            testCase.verifyTrue(isfield(ds,'sym') && isempty(ds.sym), ...
                'solveBands did not take the cache branch for the expected rib name');

            P.prefname = 'run7';
            writeSentinel(dl, ['run7_',name], struct('sym',[]));
            ds = solveBands(P);
            testCase.verifyTrue(isfield(ds,'sym') && isempty(ds.sym), 'prefname not prepended');
        end

        function solveBandsFileBaseRibSidewallAndWidth(testCase)
            testCase.assumeTrue(ispc, 'solveBands appends a literal ''\'' to P.datLoc.');
            dl = tempDatLoc(testCase);
            P = lnRibParams();  P.datLoc = dl;
            name = 'rib_a_570nm_s_254nm_w_1747nm_t_391nm_th_283nm_d_15deg_r1_0nm_r2_0nm_rxt_0_';
            writeSentinel(dl, name, struct('sym',[]));
            ds = solveBands(P);
            testCase.verifyTrue(isfield(ds,'sym') && isempty(ds.sym), ...
                'solveBands rib name does not encode the explicit w and d');
        end

        function solveBandsFileBase2DRibs(testCase)
            testCase.assumeTrue(ispc, 'solveBands appends a literal ''\'' to P.datLoc.');
            dl = tempDatLoc(testCase);
            P = ribSlabParams();  P.datLoc = dl;
            % Characterisation: the separator after w is missing in the code
            % ('..._w_3000nmhi_955nm_...'), reproduced here on purpose.
            name = '2D_ribs_a_1500nm_w_3000nmhi_955nm_wi_670nm_th_160nm_ho_975nm_wo_1000nm_ai_150nm_';
            writeSentinel(dl, name, struct('sym',[]));
            ds = solveBands(P);
            testCase.verifyTrue(isfield(ds,'sym') && isempty(ds.sym), ...
                'solveBands did not take the cache branch for the expected 2D_ribs name');
        end

        function solveBandsNoSymFileBaseRibLegacyUnchanged(testCase)
            testCase.assumeTrue(ispc, 'solveBands_noSym appends a literal ''\'' to P.datLoc.');
            dl = tempDatLoc(testCase);
            P = lnRibParams();  P.datLoc = dl;
            % The name of the existing LN runs under test/: prefix 'hole_unitCell'
            % (sic), no separators around it / after s, no fillet suffix.
            name = 'hole_unitCella_570nm_s_254nmw_1747nm_t_391nm_th_283nm_d_15deg';
            writeSentinel(dl, name, struct('full',[]));
            ds = solveBands_noSym(P);
            testCase.verifyTrue(isfield(ds,'full') && isempty(ds.full), ...
                'solveBands_noSym no longer reproduces the legacy LN rib name');
        end

        function solveBandsNoSymFileBaseRibDefaults(testCase)
            % no P.w -> w from the L(a) rule; default fillets -> suffix appended
            testCase.assumeTrue(ispc, 'solveBands_noSym appends a literal ''\'' to P.datLoc.');
            dl = tempDatLoc(testCase);
            P = struct('celltype','rib','datLoc',dl);
            name = 'hole_unitCella_300nm_s_100nmw_740nm_t_100nm_th_110nm_d_0deg_r1_40nm_r2_40nm';
            writeSentinel(dl, name, struct('full',[]));
            ds = solveBands_noSym(P);
            testCase.verifyTrue(isfield(ds,'full') && isempty(ds.full), ...
                'solveBands_noSym rib name with defaults');
        end

        function solveOpticalBandsFileBaseRibLegacyUnchanged(testCase)
            dl = tempDatLoc(testCase);
            P = lnRibParams();  P.datLoc = dl;
            name = 'optical_rib_a_570nm_s_254nmw_1747nm_t_391nm_th_283nm_d_15deg';
            writeSentinel(dl, name, struct('opticalBand','sentinel'));
            ds = solveOpticalBands(P);
            testCase.verifyEqual(ds.opticalBand, 'sentinel');
        end

        function solveOpticalBandsFileBase2DRibs(testCase)
            dl = tempDatLoc(testCase);
            P = ribSlabParams();  P.datLoc = dl;
            name = 'optical_2D_ribs_a_1500nm_w_3000nmhi_955nm_wi_670nm_ho_975nm_wo_1000nm_ai_150nm_';
            writeSentinel(dl, name, struct('opticalBand','sentinel'));
            ds = solveOpticalBands(P);
            testCase.verifyEqual(ds.opticalBand, 'sentinel');
        end

        function createFileBaseRibBothCopiesAgree(testCase)
            P = struct('celltype','rib','mevenx',1,'meveny',1,'mevenz',1);
            expected = ['th110nm_rib_a300nm_t100nm_s100nm_Lr400nm_ar240nm_dL-0.4998_', ...
                'r1_40nm_r2_40nm_sw0deg_nc14_nd5_d0.2000_quadratic_cs_ctr1_az60nm_lr_X1_Y1_Z1'];
            root = fileparts(mfilename('fullpath'));
            % root copy: call it with root as the current folder so that it,
            % not nanobeam/CreateFileBase.m, wins the lookup
            testCase.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(root));
            Pr = CreateFileBase(P);
            testCase.verifyEqual(Pr.fileBase, expected, 'root CreateFileBase');
            Pexp = P;  Pexp.w = 800e-9;
            Pw = CreateFileBase(Pexp);
            testCase.verifySubstring(Pw.fileBase, '_Lc400nm_', 'explicit P.w not encoded');
            testCase.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture( ...
                fullfile(root,'nanobeam')));
            Pn = CreateFileBase(P);
            testCase.verifyEqual(Pn.fileBase, expected, 'nanobeam/CreateFileBase');
        end

        % ---- solver -> builder routing (source inspection) ---------------
        function solverDispatchesRibCelltype(testCase, WiringCase)
            [file, celltype, builder] = WiringCase{:};
            testCase.verifyEqual(dispatchedBuilder(file, celltype), builder, ...
                sprintf('%s.m does not route celltype ''%s'' to %s', file, celltype, builder));
        end

        function noVisibleRibCelltypeLeft(testCase)
            % the device was renamed 'visible_rib' -> 'rib'; no dispatcher may
            % still branch on the old name
            files = {'runBands','runBands_noSym','solveBands','solveBands_noSym', ...
                'solveOpticalBands','runOpticalBand_1D','runOpticalBand_3D','CreateFileBase'};
            root = fileparts(mfilename('fullpath'));
            nb = {'RunNanobeamFEM','SetupNanobeamFEM','PlotDefectCells','PlotEy','CreateFileBase'};
            paths = [cellfun(@(f) fullfile(root,[f '.m']), files, 'UniformOutput', false), ...
                     cellfun(@(f) fullfile(root,'nanobeam',[f '.m']), nb, 'UniformOutput', false)];
            for k = 1:numel(paths)
                testCase.verifyFalse(contains(fileread(paths{k}), '''visible_rib'''), ...
                    sprintf('%s still references ''visible_rib''', paths{k}));
            end
        end
    end

    methods (Test, TestTags = {'PureMATLAB','KnownBug'})
        function knownBug_runBandsHas2DRibsBranch(testCase)
            % solveBands names '2D_ribs' files and test_nanobeamRibGeometry.m
            % calls solveBands with celltype '2D_ribs', but runBands.m has no
            % '2D_ribs' branch, so the final else builds a BOOMERANG STRIP
            % (buildBoomerangStrip_3D). buildRibUnitCell's only caller is
            % RunNanobeamBands.m:161, and that call is unconditional.
            expectedFailure(testCase, ...
                strcmp(dispatchedBuilder('runBands','2D_ribs'), 'buildRibUnitCell'), ...
                ['runBands.m has no ''2D_ribs'' branch; solveBands(P) with celltype ' ...
                 '''2D_ribs'' silently builds buildBoomerangStrip_3D.']);
        end

        function knownBug_optical2DRibsFileBaseEncodesThickness(testCase)
            % solveOpticalBands' '2D_ribs' name omits th, so two designs that
            % differ only in thickness share a file and the second silently
            % reloads the first's bands.
            dl = tempDatLoc(testCase);
            P = ribSlabParams();  P.datLoc = dl;
            name = 'optical_2D_ribs_a_1500nm_w_3000nmhi_955nm_wi_670nm_ho_975nm_wo_1000nm_ai_150nm_';
            writeSentinel(dl, name, struct('opticalBand','th160'));
            P.th = 200e-9;
            collided = false;
            try
                ds = solveOpticalBands(P);
                collided = isfield(ds,'opticalBand') && isequal(ds.opticalBand,'th160');
            catch
                % no cache hit -> it tried to solve and stopped on a missing
                % solver field: the names differ, which is the correct behaviour
            end
            expectedFailure(testCase, ~collided, ...
                ['solveOpticalBands ''2D_ribs'' file base omits th: runs differing ' ...
                 'only in thickness reuse each other''s cached bands.']);
        end
    end

    %% =====================================================================
    %  (b) COMSOL GEOMETRY (no mesh, no solve)
    %  =====================================================================
    methods (Test, TestTags = {'COMSOL'})

        % ---- buildRibUnitCell ('2D_ribs'), values of test_nanobeamRibGeometry.m
        function ribSlabBuildsWithAnalyticShape(testCase, SlabSymZ)
            [model, P] = buildRibSlab(testCase, SlabSymZ);
            S = geomStats(model);
            testCase.verifyGreaterThanOrEqual(S.nDom, 1);
            testCase.verifyTrue(S.connected, 'solid is not connected');
            zLo = -P.th/2*(~SlabSymZ);
            verifyBox(testCase, S, [-P.a/2 P.a/2; 0 P.w/2; zLo P.th/2]);
            % slab a x w/2 minus the inverted-U void (wo x ho minus the wi x hi stub)
            area = P.a*P.w/2 - (P.wo*P.ho - P.wi*P.hi);
            verifyVolume(testCase, model, S, area*(P.th/2 + P.th/2*(~SlabSymZ)));
        end

        % ---- buildRibUnitCell_LN ('rib'), vertical walls, defaults ---------
        function ribBuildsWithAnalyticShape(testCase, SymY, SymZ, Optical)
            P = ribDefaultParams(struct('mbeveny',SymY,'mbevenz',SymZ,'run_optical',Optical));
            [model, P] = buildRib(testCase, P);
            S = geomStats(model);
            testCase.verifyEqual(S.nDom >= 1 + Optical, true, sprintf('%d domain(s)', S.nDom));
            testCase.verifyTrue(S.connected, 'geometry is not connected');
            L = ribLength(ribDefaultParams(), P.a);
            testCase.verifyEqual(P.L, L, 'AbsTol', 1e-18);
            testCase.verifyEqual(P.w, 2*L, 'AbsTol', 1e-18);
            halfY = SymY || Optical;            % optical runs always build y >= 0
            if Optical
                R = P.airradUC;
                verifyBox(testCase, S, [-P.a/2 P.a/2; 0 R; -R*(~SymZ) R]);
                vol = pi*R^2*P.a/2 * (1 - SymZ/2);
            else
                verifyBox(testCase, S, [-P.a/2 P.a/2; -L*(~halfY) L; -P.th/2*(~SymZ) P.th/2]);
                % per y-half: half spine + rib beyond it, minus 2 convex tip
                % fillets, plus 2 concave root fillets
                halfArea = P.a*P.s/2 + P.t*(L - P.s/2) + 2*(1 - pi/4)*(P.r1^2 - P.r2^2);
                vol = halfArea*(2 - halfY) * P.th*(1 - SymZ/2);
            end
            verifyVolume(testCase, model, S, vol);
        end

        function ribProvidesSolverSelections(testCase, SymY, SymZ, Optical)
            P = ribDefaultParams(struct('mbeveny',SymY,'mbevenz',SymZ,'run_optical',Optical));
            [model, P] = buildRib(testCase, P);
            xb = resolveSel(model,'geom1_xboundaries_bnd',2);
            testCase.verifyNotEmpty(xb, 'geom1_xboundaries_bnd empty');
            testCase.verifyNotEmpty(P.xEnd1, 'P.xEnd1 empty');
            testCase.verifyNotEmpty(P.xEnd2, 'P.xEnd2 empty');
            testCase.verifyEqual(sort(xb), sort([P.xEnd1(:)' P.xEnd2(:)']), ...
                'xboundaries disagree with the bndindex lookups');
            testCase.verifyNotEmpty(P.zEnd, 'P.zEnd empty');
            testCase.verifyNotEmpty(P.domSel.beam, 'P.domSel.beam empty');
            if SymY || Optical
                testCase.verifyNotEmpty(resolveSel(model,'geom1_yboundaries_bnd',2), 'yboundaries empty');
                testCase.verifyNotEmpty(P.yEnd1, 'P.yEnd1 empty');
            end
            if SymZ
                testCase.verifyNotEmpty(resolveSel(model,'geom1_ZsymSel',2), 'geom1_ZsymSel empty');
                testCase.verifyNotEmpty(P.bndSel.Zsym, 'P.bndSel.Zsym empty');
            end
            if Optical
                testCase.verifyNotEmpty(P.domSel.air, 'P.domSel.air empty');
                testCase.verifyEmpty(intersect(P.domSel.air, P.domSel.beam), ...
                    'air and beam domain selections overlap');
                testCase.verifyNotEmpty(P.bndSel.airOuter, 'geom1_airOuter_bnd empty');
                testCase.verifyEqual(P.zEnd2, P.bndSel.airOuter, ...
                    'P.zEnd2 (runOpticalBand_3D scattering) is not the outer air surface');
            end
        end

        % ---- sidewall angle: the original LN cell ----------------------------
        function ribSidewallTrapezoidVolume(testCase, SymY, Optical)
            % LN values (explicit w, d = 15 deg, no fillets). Cross-section at
            % height u: a*s(u) + t(u)*(2L - s), linear in u, so
            % V = th*[a*(s - th*tan d) + (t - th*tan d)*(2L - s)].
            P = lnRibParams();  P.mbeveny = SymY;  P.run_optical = Optical;
            [model, P] = buildRib(testCase, P);
            S = geomStats(model);
            testCase.verifyEqual(S.nDom >= 1 + Optical, true, sprintf('%d domain(s)', S.nDom));
            testCase.verifyTrue(S.connected, 'geometry is not connected');
            testCase.verifyEqual(P.L, P.w/2, 'explicit P.w did not set L');
            halfY = SymY || Optical;
            if Optical
                R = P.airradUC;
                verifyBox(testCase, S, [-P.a/2 P.a/2; 0 R; -R R]);
                vol = pi*R^2*P.a/2;
            else
                verifyBox(testCase, S, [-P.a/2 P.a/2; -P.L*(~halfY) P.L; -P.th/2 P.th/2]);
                k = tan(P.d);
                vol = P.th*(P.a*(P.s - k*P.th) + (P.t - k*P.th)*(2*P.L - P.s)) / (1 + halfY);
            end
            verifyVolume(testCase, model, S, vol);
            testCase.verifyNotEmpty(resolveSel(model,'geom1_xboundaries_bnd',2), 'xboundaries empty');
        end

        function ribExplicitWidthWinsInGeometry(testCase)
            P = ribDefaultParams(struct('w',600e-9,'mbeveny',1,'mbevenz',0));
            [model, P] = buildRib(testCase, P);
            S = geomStats(model);
            testCase.verifyEqual(P.L, 300e-9, 'AbsTol', 1e-18);
            testCase.verifyEqual(S.box(2,2), 300e-9, 'AbsTol', 1e-10, 'rib tip not at w/2');
        end
    end

    methods (Test, TestTags = {'COMSOL','KnownBug'})
        % buildRibUnitCell draws the slab over x in [-a/2, a/2]
        % (pos [-beamLen/2 0]) but places several selections as if the cell
        % ran from x = 0 to x = nperiod*a, the nanobeam convention. It also ends
        % in runAll and reads selections after runCurrent, and both leave the
        % geometry unfinalized. Each consequence gets its own test, so the
        % report lists them separately.
        function knownBug_ribSlabXsymSelAtCellEdge(testCase)
            [~, P] = buildRibSlab(testCase, 1);
            expectedFailure(testCase, ~isempty(P.bndSel.Xsym), ...
                ['buildRibUnitCell XsymSel is a box at x = 0, the middle of the ' ...
                 'cell, where no face exists; P.bndSel.Xsym (a periodic face for ' ...
                 'RunNanobeamBands) is empty.']);
        end

        function knownBug_ribSlabXendSelAtCellEdge(testCase)
            [~, P] = buildRibSlab(testCase, 1);
            expectedFailure(testCase, ~isempty(P.bndSel.Xend), ...
                ['buildRibUnitCell XendSel is at x = nperiod*a while the cell ends ' ...
                 'at x = a/2; P.bndSel.Xend is empty.']);
        end

        function knownBug_ribSlabYsymSelCoversWholeFace(testCase)
            [~, P] = buildRibSlab(testCase, 1);
            expectedFailure(testCase, ~isempty(P.bndSel.Ysym), ...
                ['buildRibUnitCell YsymSel box spans x in [0, nperiod*a] but the ' ...
                 'y = 0 face spans [-a/2, a/2]; with condition allvertices ' ...
                 'P.bndSel.Ysym is empty.']);
        end

        function knownBug_ribSlabZsymSelBoxCoversBeamWidth(testCase)
            [model, P] = buildRibSlab(testCase, 1);
            ok = ~isempty(P.bndSel.Zsym) && ...
                 ~isempty(resolveSel(model,'geom1_ZsymSel',2));
            expectedFailure(testCase, ok, ...
                ['buildRibUnitCell ZsymSel box limits y to +-nperiod*a/2 (750 nm) ' ...
                 'but the z = 0 face reaches y = w/2 (1500 nm); with allvertices ' ...
                 'it selects nothing. P.bndSel.Zsym is also read right after ' ...
                 'runCurrent, before the geometry is finalized.']);
        end

        function knownBug_ribSlabIndexLookupsAfterFinalize(testCase)
            [~, P] = buildRibSlab(testCase, 1);
            ok = ~isempty(P.xEnd1) && ~isempty(P.xEnd2) && ~isempty(P.yEnd1);
            expectedFailure(testCase, ok, ...
                ['buildRibUnitCell ends in runAll (never runs fin), so its ' ...
                 'bndindex lookups P.xEnd1/xEnd2/yEnd1 see 0 vertices and return [].']);
        end
    end
end

%% ========================================================================
%  Local helpers
%  ========================================================================

function P = ribSlabParams()
% buildRibUnitCell values from test_nanobeamRibGeometry.m
P = struct('celltype','2D_ribs','xsect','rect','beamMat','diamond', ...
    'a',1500e-9,'w',3000e-9,'wi',670e-9,'ho',975e-9,'hi',955e-9, ...
    'wo',1000e-9,'th',160e-9,'ai',150e-9,'nperiod',1,'holeatedge',0, ...
    'mbevenz',1,'mbeveny',1);
end

function P = lnRibParams()
% values of the original LN rib design (old test_LN_unitCell.m, now
% test_rib_unitCell_explicitW.m / test_rib_unitCell_LN.m): explicit w,
% sidewall angle 15 deg in radians, no fillet fields
P = struct('celltype','rib','xsect','rect', ...
    'a',570e-9,'s',254e-9,'w',1747e-9,'th',283e-9,'t',391e-9, ...
    'd',15*pi/180,'run_optical',0,'mbevenz',0,'mbeveny',0);
end

function [model, P] = buildRibSlab(testCase, mbevenz)
P = ribSlabParams();
P.mbevenz = mbevenz;
model = newComsolModel(testCase);
[model, P] = buildRibUnitCell(model, P);
end

function [model, P] = buildRib(testCase, P)
model = newComsolModel(testCase);
[model, P] = buildRibUnitCell_LN(model, P);
end

% ---- COMSOL session ------------------------------------------------------

function tf = livelinkAvailable()
% True when ModelUtil answers, connecting with mphstart to an ALREADY-RUNNING
% server if needed. Never starts a server. Cached for the MATLAB session.
persistent cached
if ~isempty(cached)
    tf = cached;
    return
end
tf = false;
try
    com.comsol.model.util.ModelUtil.tags();
    tf = true;
catch
    if exist('mphstart','file') == 2
        try
            mphstart();
            com.comsol.model.util.ModelUtil.tags();
            tf = true;
        catch
        end
    end
end
cached = tf;
end

function model = newComsolModel(testCase)
testCase.assumeTrue(livelinkAvailable(), ...
    ['COMSOL LiveLink not reachable (ModelUtil unavailable and no running ' ...
     'server to connect to); COMSOL geometry test skipped.']);
figsBefore = findall(groot, 'Type', 'figure');
testCase.addTeardown(@() closeNewFigures(figsBefore));   % builders plot unconditionally
com.comsol.model.util.ModelUtil.clear();
model = com.comsol.model.util.ModelUtil.create('model');
testCase.addTeardown(@() com.comsol.model.util.ModelUtil.clear());
end

function closeNewFigures(figsBefore)
figsNow = findall(groot, 'Type', 'figure');
close(setdiff(figsNow, figsBefore));
end

% ---- geometry measurements ------------------------------------------------

function S = geomStats(model)
% Finalize (the solvers call geom.run themselves before adding physics; the
% builder's own returned P fields are checked BEFORE this, in the tests),
% then count domains, test connectivity and take the bounding box.
g = model.geom('geom1');
g.run();
S.nDom = double(g.getNDomains());
ud = double(g.getUpDown());                 % 2 x nBoundaries, 0 = outside
S.connected = domainsConnected(ud, S.nDom);
[S.box, S.boxExact] = boundingBox(g);
end

function tf = domainsConnected(ud, nDom)
if nDom <= 1
    tf = nDom == 1;
    return
end
pairs = ud(:, all(ud > 0, 1) & ud(1,:) ~= ud(2,:));
G = graph(pairs(1,:), pairs(2,:), [], nDom);
tf = all(conncomp(G) == 1);
end

function [box, exact] = boundingBox(g)
% GeomSequence.getBoundingBox where the API offers it; otherwise the vertex
% hull, which is exact only for planar faces. A curved face (air cylinder)
% need not carry a vertex at its extreme, so verifyBox relaxes those checks.
try
    bb = double(g.getBoundingBox());
    box = reshape(bb(:), 2, []).';
    exact = true;
catch
    vc = double(g.getVertexCoord());
    box = [min(vc,[],2) max(vc,[],2)];
    exact = false;
end
end

function verifyBox(testCase, S, expected)
tol = 1e-10;                                % 0.1 nm
if S.boxExact
    testCase.verifyEqual(S.box, expected, 'AbsTol', tol, 'bounding box');
else
    % vertex hull: never larger than the true box, and planar extremes exact
    testCase.verifyLessThanOrEqual(S.box(:,2), expected(:,2) + tol, 'vertex hull exceeds box');
    testCase.verifyGreaterThanOrEqual(S.box(:,1), expected(:,1) - tol, 'vertex hull exceeds box');
    testCase.verifyEqual(S.box(1,:), expected(1,:), 'AbsTol', tol, 'x extent');
end
end

function verifyVolume(testCase, model, S, expected)
testCase.assumeTrue(exist('mphmeasure','file') == 2, 'mphmeasure (LiveLink) not on the path');
v = 0;
for d = 1:S.nDom
    v = v + mphmeasure(model, 'geom1', 'volume', 'selection', d);
end
testCase.verifyEqual(v, expected, 'RelTol', 1e-6, 'geometry volume vs analytic');
end

function inds = resolveSel(model, tag, dim)
try
    inds = double(model.selection(tag).entities(dim));
catch
    try
        inds = double(model.selection(tag).inputEntities());
    catch
        inds = [];
    end
end
inds = inds(:)';
end

% ---- pure-MATLAB helpers ---------------------------------------------------

function expectedFailure(testCase, correctBehaviour, bugText)
% KnownBug semantics: Incomplete while the bug is present, pass once fixed.
if correctBehaviour
    testCase.log(1, ['Known bug no longer reproduces; update RibUnitCellTest: ', bugText]);
else
    testCase.assumeFail(['KNOWN BUG (expected failure): ', bugText]);
end
end

function name = dispatchedBuilder(file, celltype)
% Builder that FILE.m calls in its strcmp(P.celltype,'<celltype>') branch, ''
% if there is no such branch. Comment lines are stripped first, so a
% commented-out branch does not count. The branch text runs to the next
% elseif/else at the start of a line.
src = regexprep(fileread(which(file)), '(?m)^[ \t]*%[^\n]*', '');
idx = regexp(src, ['strcmp\(P\.celltype\s*,\s*''', celltype, '''\)'], 'once');
name = '';
if isempty(idx)
    return
end
rest = src(idx:end);
stop = regexp(rest, '\n[ \t]*(elseif|else)\>', 'once');
if ~isempty(stop)
    rest = rest(1:stop);
end
tok = regexp(rest, '\[\s*model\s*,\s*P\s*\]\s*=\s*(\w+)\s*\(', 'tokens', 'once');
if ~isempty(tok)
    name = tok{1};
end
end

function dl = tempDatLoc(testCase)
f = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
dl = [f.Folder, filesep];
end

function writeSentinel(datLoc, fileBase, ds)  % ds is saved by name
save([datLoc, fileBase, '_bds.mat'], 'ds');
end
