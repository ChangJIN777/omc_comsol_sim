classdef BandTrackingTest < matlab.unittest.TestCase
%BANDTRACKINGTEST  Tests for the rib band-tracking tool.
%
%   results = runtests('BandTrackingTest','Tag','PureMATLAB');   % no COMSOL
%   results = runtests('BandTrackingTest','Tag','COMSOL');       % 3 small solves
%
%   Covers ribModeGrid, bandOverlap, modeFingerprints, connectBandsAlongK,
%   bandTrackCandidate, matchBandModes, surrogateRibBandsTracked,
%   trackRibBands (dry run, bridging, pre-solved inputs, outputs) and
%   plotRibBandTrack.

    methods (TestClassSetup)
        function addThisFolderToPath(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fileparts(mfilename('fullpath'))));
        end
    end

    methods (Test, TestTags = {'PureMATLAB'})

        % ---------------- building blocks ----------------
        function gridIsDesignIndependent(testCase)
            A = ribModeGrid(surrP(100e-9));
            B = ribModeGrid(struct('t',140e-9,'s',70e-9,'L',420e-9,'a',320e-9, ...
                                   'mbeveny',1,'mbevenz',1));
            testCase.verifyEqual(A.q, B.q, 'AbsTol', 1e-12);
            testCase.verifyEqual(A.region, B.region);
            testCase.verifyNotEqual(A.xyz, B.xyz);
            testCase.verifyEqual(sum(A.w), 1, 'AbsTol', 1e-12);
            n = size(A.q, 2);
            testCase.verifyEqual(A.mirrorX(A.mirrorX), 1:n, 'mirrorX is not an involution');
            testCase.verifyEmpty(A.mirrorY);
            testCase.verifyEmpty(A.mirrorZ);
            testCase.verifyLessThan(max(abs(A.xyz(2,:))), 370e-9);
            testCase.verifyGreaterThan(min(A.xyz(2,:)), 0);
        end

        function gridFullCellHasAllMirrors(testCase)
            G = ribModeGrid(struct('mbeveny',0,'mbevenz',0));
            n = size(G.q, 2);
            testCase.verifyEqual(n, 4*size(ribModeGrid(surrP(100e-9)).q, 2));
            testCase.verifyEqual(G.mirrorY(G.mirrorY), 1:n);
            testCase.verifyEqual(G.mirrorZ(G.mirrorZ), 1:n);
        end

        function overlapIgnoresSignPhaseAndScale(testCase)
            rs = RandStream('mt19937ar','Seed',1);
            U = rs.randn(30, 4) + 1i*rs.randn(30, 4);
            V = U.*[-1, exp(0.7i), 3, -2i];
            M = bandOverlap(U, V, rs.rand(10,1) + 0.1);
            testCase.verifyEqual(diag(M), ones(4,1), 'AbsTol', 1e-12);
            testCase.verifyLessThan(max(M(~eye(4))), 0.9);
        end

        function overlapSkipsPointsMissingInEitherDesign(testCase)
            rs = RandStream('mt19937ar','Seed',2);
            U = rs.randn(30, 3);  V = U;
            V(4,:) = NaN;                       % point 4, component x
            U(27,2) = NaN;                      % point 7, component z, one band
            [M, nValid] = bandOverlap(U, V);
            testCase.verifyEqual(nValid, 8);
            testCase.verifyEqual(diag(M), ones(3,1), 'AbsTol', 1e-12);
        end

        function subspaceOverlapSeesThroughDegenerateRotation(testCase)
            [Q, ~] = qr(randn(RandStream('mt19937ar','Seed',3), 30, 3), 0);
            Rot = [1 1; -1 1]/sqrt(2);
            B = Q(:,1:2)*Rot;
            M = bandOverlap(Q(:,1), B);
            testCase.verifyEqual(M, [0.5 0.5], 'AbsTol', 1e-12);
            testCase.verifyEqual(bandOverlap(Q(:,1), B, [], 'subspace'), 1, 'AbsTol', 1e-12);
            testCase.verifyEqual(bandOverlap(Q(:,3), B, [], 'subspace'), 0, 'AbsTol', 1e-12);
        end

        function fingerprintsOfSurrogateModes(testCase)
            ds = surrogateRibBandsTracked(surrP(100e-9), struct('coupling',0));
            md = ds.sym.modes;
            FP = modeFingerprints(md.U(:,:,end), md.grid);
            testCase.verifyEqual(sum(FP.pol, 1), ones(1,10), 'AbsTol', 1e-5);
            testCase.verifyEqual(sum(FP.reg, 1), ones(1,10), 'AbsTol', 1e-5);
            % even-numbered characters are x-odd, the rest x-even
            expected = ones(1,10);  expected(2:2:end) = -1;
            testCase.verifyEqual(FP.px, expected(ds.sym.truth(end,:)), 'AbsTol', 1e-4);
            testCase.verifyTrue(all(isnan(FP.py)) && all(isnan(FP.pz)));
            testCase.verifyTrue(all(isnan(FP.dil)), 'no divU was supplied');
        end

        function surrogateHybridOverlapIsCosSquared(testCase)
            % The premise of the bridging tests: characters 1 and 3 mix by an
            % angle theta, and the weighted MAC between two designs is
            % cos^2(dTheta) even though their grids and weights differ.
            [sV, ~] = crossingS();
            g = 0.1e9;
            A = bandTrackCandidate(surrogateRibBandsTracked(surrP(sV - 5e-9)).sym);
            B = bandTrackCandidate(surrogateRibBandsTracked(surrP(sV)).sym);
            thA = 0.5*atan2(2*g, detuneX(sV - 5e-9));
            M = bandOverlap(A.U, B.U, {A.w, B.w});
            [~, iA] = min(abs(A.f - lowerBranch(sV - 5e-9, g)));
            [~, iB] = min(abs(B.f - lowerBranch(sV, g)));
            testCase.verifyEqual(M(iA, iB), cos(pi/4 - thA)^2, 'AbsTol', 1e-4);
        end

        function connectRecoversWindowSlide(testCase)
            kn = (0:10)'/10;
            Fall = (1:14).*(1 + 0.15*kn)*1e9;           % 14 non-crossing bands
            sh = [0 0 0 1 1 1 1 2 2 2 2]';              % window slides upward
            F = zeros(11, 10);
            for r = 1:11, F(r,:) = Fall(r, sh(r) + (1:10)); end
            S = connectBandsAlongK(F);                  % frequency fallback
            testCase.verifyEqual(S.off, sh);
            testCase.verifyEqual(S.Flev(:,3:10), Fall(:,3:10));
            testCase.verifyEqual(S.complete, [false false true(1,8) false false]);
        end

        function connectFollowsCharacterThroughCrossingAlongK(testCase)
            % s = 190 nm, no coupling: the flat intruder cuts through the
            % valence band along k
            ds = surrogateRibBandsTracked(surrP(190e-9), struct('coupling',0));
            md = ds.sym.modes;
            S = connectBandsAlongK(ds.sym.F, md.U, md.grid.w);
            tr = ds.sym.truth;
            testCase.verifyGreaterThan(numel(unique(sum((tr == 1).*(1:10), 2))), 1, ...
                'test premise: character 1 should change column along k');
            for j = 1:size(tr, 2)
                for r = 1:size(tr, 1)
                    testCase.verifyEqual(tr(r, S.branch(r,:) == j), tr(end, j), ...
                        sprintf('branch %d lost its character at row %d', j, r));
                end
            end
            testCase.verifyEqual(S.off, zeros(11,1));
            testCase.verifyGreaterThan(min(S.conf(1:end-1)), 0.8);
        end

        function candidateNeedsModes(testCase)
            testCase.verifyError(@() bandTrackCandidate(struct('F', ones(3))), ...
                'bandTrackCandidate:noModes');
            testCase.verifyError(@() bandTrackCandidate([]), 'bandTrackCandidate:noModes');
        end

        % ---------------- matchBandModes ----------------
        function matchIgnoresSignFlipsAndPhases(testCase)
            [ref, t] = refCand();
            c = ref;
            c.U = c.U.*[-1, 1i, exp(2.1i), -1, 1, -1i, 3, -0.5, 1, -1];
            R = matchBandModes(ref, t, c);
            testCase.verifyEqual(R.idx, t);
            testCase.verifyEqual(R.status, {'confident','confident'});
            testCase.verifyEqual(R.mac, [1 1], 'AbsTol', 1e-6);
            testCase.verifyTrue(R.ok);
        end

        function matchSeesThroughDegenerateRotation(testCase)
            % The upper band and its neighbour become degenerate and the
            % solver returns an arbitrary rotation of the pair.
            [ref, t] = refCand();
            c = ref;
            j = t(2);  j2 = j + 1;
            c.f(j2) = c.f(j)*(1 + 1e-4);
            c.U(:,[j j2]) = c.U(:,[j j2])*[1 1; -1 1]/sqrt(2);
            R = matchBandModes(ref, t, c);
            testCase.verifyTrue(R.degenerate(2));
            testCase.verifyEqual(R.groupMac(2), 1, 'AbsTol', 1e-6);
            testCase.verifyLessThan(R.mac(2), 0.8);
            testCase.verifyEqual(R.status{2}, 'confident');
            testCase.verifyEqual(R.status{1}, 'confident');
        end

        function matchVetoesParityFlipUnlessDisabled(testCase)
            [ref, t] = refCand();
            c = ref;
            c.px(t(1)) = -c.px(t(1));           % same shape, flipped parity flag
            R = matchBandModes(ref, t, c);
            testCase.verifyEqual(R.status{1}, 'lost');
            testCase.verifyTrue(isnan(R.idx(1)));
            testCase.verifyEqual(R.status{2}, 'confident');
            R = matchBandModes(ref, t, c, struct('useParity', false));
            testCase.verifyEqual(R.status, {'confident','confident'});
        end

        function matchRejectsBadTargets(testCase)
            ref = refCand();
            testCase.verifyError(@() matchBandModes(ref, [1 NaN], ref), ...
                'matchBandModes:badTargets');
        end

        % ---------------- trackRibBands (dry run) ----------------
        function trackFollowsBothBandsThroughExactCrossings(testCase)
            % No coupling: the intruder sweeps down through the conduction
            % band, the gap and the valence band. The band ORDER changes,
            % the two tracked modes do not.
            f = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            Ps = sweepS((103:10:193)*1e-9);
            res = trackRibBands(Ps, struct('dryRun',true, 'targetFreq',11.2e9, ...
                'rootDir',fullfile(f.Folder,'run'), 'surrogate',struct('coupling',0)));
            T = res.table;
            testCase.verifyEqual(res.vars, {'s'});
            testCase.verifyEqual(T.status{1}, 'init');
            testCase.verifyTrue(all(strcmp(T.status(2:end), 'confident')));
            testCase.verifyEqual(T.truth_lo', ones(1,10), 'the lower band was swapped');
            testCase.verifyEqual(T.truth_up', 2*ones(1,10), 'the upper band was swapped');
            testCase.verifyGreaterThan(numel(unique(T.idx_up)), 1, ...
                'test premise: the band order should change in this sweep');
            testCase.verifyEqual(res.nBridgeSolves, 0);
            testCase.verifyEqual(T.gapX, T.f_up - T.f_lo, 'AbsTol', 1);
            testCase.verifyTrue(any(T.nBetween > 0), ...
                'test premise: the intruder should sit between the bands somewhere');
            testCase.verifyGreaterThan(min(T.macGold_lo), 0.99);
        end

        function bridgingResolvesAvoidedCrossing(testCase)
            % 100 MHz coupling. The step across the anticrossing rotates the
            % hybrid by about 48 deg (MAC 0.45 / 0.55): ambiguous. One
            % midpoint solve splits it into two steps of about 24 deg (MAC
            % 0.84), which are confident.
            f = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            [sV, fEdge] = crossingS();
            Ps = sweepS(sV + (-25:10:25)*1e-9);
            o = struct('dryRun',true, 'initFreq',fEdge, ...
                'rootDir',fullfile(f.Folder,'bridged'));
            res = trackRibBands(Ps, o);
            T = res.table;
            testCase.verifyEqual(T.status', ...
                {'init','confident','confident','bridged','confident','confident'});
            testCase.verifyEqual(res.nBridgeSolves, 1);
            testCase.verifyEqual(T.nBridge', [0 0 0 1 0 0]);
            testCase.verifyEqual(res.bridges(1).design, 4);
            testCase.verifyEqual(res.bridges(1).status, 'confident');
            testCase.verifyEqual(res.bridges(1).x, sV, 'RelTol', 1e-9);
            % the tracker stays on the continuous (lower) branch, whose
            % character turns from valence (1) into intruder (3)
            testCase.verifyEqual(T.truth_lo([1 2 3 5 6])', [1 1 1 3 3]);
            testCase.verifyEqual(T.truth_up', 2*ones(1,6));
            testCase.verifyTrue(all(diff(T.f_lo) < 0), 'the lower branch should fall');
            testCase.verifyLessThan(T.macGold_lo(6), 0.1);

            % without bridging the same step is reported, not resolved
            o.maxBridge = 0;  o.rootDir = fullfile(f.Folder,'unbridged');
            res0 = trackRibBands(Ps, o);
            testCase.verifyEqual(res0.table.status{4}, 'ambiguous');
            testCase.verifyEqual(res0.nBridgeSolves, 0);
            testCase.verifyEqual(height(res0.table), 6);
            testCase.verifyFalse(isnan(res0.table.f_lo(4)), 'best guess should be reported');
        end

        function maxBridgeExhaustionMarksAmbiguousAndContinues(testCase)
            % 30 MHz coupling: the hybrid region is about 2 nm wide and the
            % second design sits in its centre, so halving the step twice is
            % not enough (30.7 deg left); a third level is.
            f = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            [sV, fEdge] = crossingS();
            o = struct('dryRun',true, 'initFreq',fEdge, 'maxBridge',2, ...
                'surrogate',struct('coupling',0.03e9), ...
                'rootDir',fullfile(f.Folder,'two'));
            res = trackRibBands(sweepS(sV + [-10 0 20]*1e-9), o);
            T = res.table;
            testCase.verifyEqual(T.status', {'init','ambiguous','confident'});
            testCase.verifyEqual(res.nBridgeSolves, 2);
            testCase.verifyEqual(T.nBridge', [0 2 0]);
            testCase.verifyEqual(T.status_up{2}, 'confident');
            % design 3 is matched from design 1 (the last trusted design):
            % the narrow anticrossing is jumped and the valence character kept
            testCase.verifyEqual(T.truth_lo(3), 1);
            testCase.verifyEqual(T.step(3), 1, 'AbsTol', 1e-9);

            o.maxBridge = 3;  o.rootDir = fullfile(f.Folder,'three');
            res = trackRibBands(sweepS(sV + [-10 0]*1e-9), o);
            testCase.verifyEqual(res.table.status', {'init','bridged'});
            testCase.verifyEqual(res.nBridgeSolves, 3);
            testCase.verifyEqual([res.bridges.depth], [1 2 3]);
        end

        function presolvedStructsAndFilesGiveTheSameTrack(testCase)
            f = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            [sV, fEdge] = crossingS();
            sList = sV + (-25:10:25)*1e-9;
            o = struct('dryRun',true, 'initFreq',fEdge, ...
                'rootDir',fullfile(f.Folder,'fromP'));
            resP = trackRibBands(sweepS(sList), o);

            dsList = cell(1, numel(sList));
            files  = cell(1, numel(sList));
            for k = 1:numel(sList)
                ds = surrogateRibBandsTracked(surrP(sList(k)));
                dsList{k} = ds;
                files{k} = fullfile(f.Folder, sprintf('d%d_bds.mat', k));
                save(files{k}, 'ds');
            end
            o.rootDir = fullfile(f.Folder,'fromDs');
            resD = trackRibBands(dsList, o);       % bridge solves: surrogate
            o.rootDir = fullfile(f.Folder,'fromFiles');
            resF = trackRibBands(files, o);
            testCase.verifyEqual(resD.table.status, resP.table.status);
            testCase.verifyEqual(resF.table.status, resP.table.status);
            testCase.verifyEqual(resD.table.f_lo, resP.table.f_lo);
            testCase.verifyEqual(resF.table.f_up, resP.table.f_up);
            testCase.verifyTrue(all(strcmp(resD.table.source, 'given')));
            testCase.verifyTrue(all(strcmp(resF.table.source, 'file')));
            testCase.verifyEqual(resF.nBridgeSolves, 1);
        end

        function failedDesignIsReportedAndSkipped(testCase)
            f = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            good = surrogateRibBandsTracked(surrP(103e-9), struct('coupling',0));
            bad  = surrogateRibBandsTracked(surrP(113e-9), struct('coupling',0));
            bad.sym = rmfield(bad.sym, 'modes');
            last = surrogateRibBandsTracked(surrP(123e-9), struct('coupling',0));
            o = struct('dryRun',true, 'targetFreq',11.2e9, ...
                'rootDir',fullfile(f.Folder,'run'), 'surrogate',struct('coupling',0));
            res = testCase.verifyWarning(@() trackRibBands({good, bad, last}, o), ...
                'trackRibBands:designFailed');
            testCase.verifyEqual(res.table.status', {'init','failed','confident'});
            testCase.verifyNotEmpty(res.table.error{2});
            testCase.verifyTrue(isnan(res.table.f_lo(2)));
        end

        function outputsAreWritten(testCase)
            f = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            res = trackRibBands(sweepS((103:10:133)*1e-9), struct('dryRun',true, ...
                'targetFreq',11.2e9, 'rootDir',fullfile(f.Folder,'run'), ...
                'surrogate',struct('coupling',0)));
            testCase.verifyTrue(isfile(res.csvFile));
            testCase.verifyTrue(isfile(res.matFile));
            Tc = readtable(res.csvFile, 'TextType', 'string');
            testCase.verifyEqual(height(Tc), 4);
            testCase.verifyTrue(all(ismember({'design','s_nm','status','f_lo','f_up', ...
                'gapX','mac_lo','mac_up','nBridge'}, Tc.Properties.VariableNames)));
            testCase.verifyEqual(Tc.s_nm', 103:10:133, 'AbsTol', 1e-6);
            Ld = load(res.matFile);
            testCase.verifyEqual(Ld.res.table.f_lo, res.table.f_lo);
            testCase.verifyEqual(numel(Ld.res.designs), 4);
        end

        function initialBandChoiceIsValidated(testCase)
            f = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            o = struct('dryRun',true, 'rootDir',fullfile(f.Folder,'run'));
            Ps = sweepS([103 113]*1e-9);
            testCase.verifyError(@() trackRibBands(Ps, o), 'trackRibBands:noTarget');
            o.initIndex = [3 3];
            testCase.verifyError(@() trackRibBands(Ps, o), 'trackRibBands:badInitIndex');
            testCase.verifyError(@() trackRibBands({}, o), 'trackRibBands:badDesigns');
        end

        function summaryPlotDraws(testCase)
            f = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            [sV, fEdge] = crossingS();
            res = trackRibBands(sweepS(sV + (-25:10:25)*1e-9), struct('dryRun',true, ...
                'initFreq',fEdge, 'rootDir',fullfile(f.Folder,'run')));
            fig = figure('Visible', 'off');
            testCase.addTeardown(@() close(fig));
            ax = axes(fig);
            plotRibBandTrack(res, ax);
            names = get(findobj(ax, 'Type', 'line'), 'DisplayName');
            testCase.verifyTrue(any(contains(names, 'bridged')));
            testCase.verifyTrue(any(contains(names, 'below the gap')));
        end
    end

    methods (Test, TestTags = {'COMSOL'})
        function smallRealSweep(testCase)
            % Three small solves (t = 100, 105, 110 nm; 3 k-points, 6 bands,
            % coarse mesh), then the same track again from the saved files.
            testCase.assumeTrue(livelinkReachable(), ...
                'COMSOL LiveLink not reachable; real-sweep test skipped.');
            f = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            figsBefore = findall(groot, 'Type', 'figure');
            testCase.addTeardown(@() close(setdiff(findall(groot,'Type','figure'), figsBefore)));
            testCase.addTeardown(@() com.comsol.model.util.ModelUtil.clear());

            Ps = arrayfun(@comsolP, [100 105 110]*1e-9);
            o = struct('targetFreq',11.2e9, 'maxBridge',1, ...
                'rootDir',fullfile(f.Folder,'solve'));
            res = trackRibBands(Ps, o);
            T = res.table;
            testCase.verifyEqual(T.status{1}, 'init');
            testCase.verifyTrue(all(ismember(T.status(2:3), {'confident','bridged'})), ...
                sprintf('statuses: %s', strjoin(T.status', ', ')));
            testCase.verifyGreaterThan(min([T.mac_lo(2:3); T.mac_up(2:3)]), 0.8);
            testCase.verifyTrue(all(T.f_up > T.f_lo));

            % the exported modes of the first solve
            d = dir(fullfile(res.designs(1).folder, '*_bds.mat'));
            testCase.assertNumElements(d, 1);
            Ld = load(fullfile(d(1).folder, d(1).name));
            md = Ld.ds.sym.modes;
            nPts = size(md.grid.xyz, 2);
            testCase.verifySize(md.U, [3*nPts, 6, 3]);
            testCase.verifyLessThan(max(md.freqCheck), 1e-6, ...
                'sampled modes are not in the column order of ds.F');
            testCase.verifyGreaterThan(mean(all(isfinite(md.U(:,:,end)), 2)), 0.9, ...
                'too many grid points fall outside the solid');
            M = bandOverlap(md.U(:,:,end), md.U(:,:,end), md.grid.w);
            testCase.verifyLessThan(max(M(~eye(6))), 0.3, ...
                'distinct X-point modes should be nearly orthogonal on the grid');
            FP = modeFingerprints(md.U(:,:,end), md.grid, md.divU(:,:,end));
            testCase.verifyGreaterThan(min(abs(FP.px)), 0.8, ...
                'X-point modes should have a definite x parity at rxtal = 0');
            if ~ispc
                testCase.verifyFalse(isfolder([res.designs(1).folder, filesep, '\']), ...
                    'stray backslash folder was not removed');
            end

            % same track from the saved results, no new solves
            files = cell(1, 3);
            for k = 1:3
                dk = dir(fullfile(res.designs(k).folder, '*_bds.mat'));
                files{k} = fullfile(dk(1).folder, dk(1).name);
            end
            o.rootDir = fullfile(f.Folder,'fromFiles');  o.maxBridge = 0;
            res2 = trackRibBands(files, o);
            testCase.verifyEqual(res2.nSolves, 0);
            testCase.verifyEqual(res2.table.f_lo(1), T.f_lo(1));
            testCase.verifyEqual(res2.table.idx_up(1), T.idx_up(1));
        end
    end
end

% =========================================================================

function P = surrP(s)
% Quarter-model surrogate design; only s is varied by the sweeps.
P = struct('t',100e-9, 's',s, 'L',370e-9, 'kpts',10, 'nbands',10, ...
           'mbeveny',1, 'mbevenz',1);
end

function Ps = sweepS(sList)
Ps = arrayfun(@surrP, sList);
end

function d = detuneX(s)
% Uncoupled intruder-minus-valence frequency at X (Hz).
ds = surrogateRibBandsTracked(surrP(s), struct('coupling',0));
d = ds.sym.detune(end);
end

function fl = lowerBranch(s, g)
% X-point frequency of the lower branch of the coupled valence/intruder pair.
ds = surrogateRibBandsTracked(surrP(s), struct('coupling',0));
f1 = ds.sym.F(end, ds.sym.truth(end,:) == 1);
f3 = ds.sym.F(end, ds.sym.truth(end,:) == 3);
fl = (f1 + f3)/2 - sqrt(((f1 - f3)/2)^2 + g^2);
end

function [sV, fEdge] = crossingS()
% Spine width at which the intruder meets the valence band at X, and the
% X-point frequencies [valence conduction] of the two gap-edge bands.
sV = fzero(@detuneX, [150e-9 199e-9]);
ds = surrogateRibBandsTracked(surrP(150e-9), struct('coupling',0));
fEdge = [ds.sym.F(end, ds.sym.truth(end,:) == 1), ...
         ds.sym.F(end, ds.sym.truth(end,:) == 2)];
end

function [ref, t] = refCand()
% X-point modes of one uncoupled surrogate design, and the sorted columns of
% its valence (1) and conduction (2) characters.
ds = surrogateRibBandsTracked(surrP(103e-9), struct('coupling',0));
ref = bandTrackCandidate(ds.sym);
tr = ds.sym.truth(end, ref.order);
t = [find(tr == 1), find(tr == 2)];
end

function P = comsolP(t)
% Small real solve: y-even / z-even sector on a quarter model.
P = ribDefaultParams();
P.t = t;
P.kpts = 2;  P.nbands = 6;  P.freq = 11.2e9;
P.TwoSymPlanes = 0;  P.zSymCondition = 1;  P.solveasym = 0;
P.mbeveny = 1;  P.mbevenz = 0;
P.completeBandGaps = 1;  P.bandStruct_2D = 0;  P.run_optical = 0;
P.plotgeom = 0;  P.savedat = 1;  P.savebndplot = 0;  P.saveplots = 0;  P.saveMPH = 0;
P.meshSize = 5;  P.fixed_bc = 0;  P.max_dof = 3e6;
P.anisoMat = 1;  P.rxtal = 0;
end

function tf = livelinkReachable()
% As RibUnitCellTest's livelinkAvailable: never starts a server.
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
end
