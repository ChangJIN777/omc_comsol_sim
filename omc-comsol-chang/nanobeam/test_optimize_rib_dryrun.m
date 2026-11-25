function test_optimize_rib_dryrun()
%TEST_OPTIMIZE_RIB_DRYRUN  Run optimize_rib_oblong_maxdef end to end, no COMSOL.
%
%   The rib counterpart of test_bayesopt_wiring / test_OptimStore_resume.
%   Instead of re-implementing the optimizer's adapters, it calls the REAL
%   script with OPT.dryRun = 1, which swaps RunNanobeamFEM for an analytic
%   stand-in and writes to tempdir. Everything else is the production code
%   path: isFabricableRib as XConstraintFcn, the bayesopt and Nelder-Mead
%   adapters, OptimStore, seedBayesFromStore, the CSV and the plots.
%
%   What it asserts:
%     1. A cold Bayesian study evaluates the start design first, spends
%        exactly its budget, and stores fitness == fitnessRaw (gate off)
%        together with a finite Q_mech estimate.
%     2. Re-running with the same runId and a larger budget solves only the
%        difference, and never reports a worse incumbent.
%     3. A Nelder-Mead session on the same runId reuses the stored start
%        point (cross-optimizer resume).
%     4. With OPT.useQmechGate = 1, gated rows have fitness 0, a positive
%        ungated fitness, Q_mech below the threshold and status Qmech_gated.
%     5. A mistyped option name is an error.
%     6. The coupling breakdown (gMB, gPE) travels with gOM into the store,
%        the CSV and stats.bestCpl, on cold, resumed and cached paths.
%
%   Requires the Statistics and Machine Learning Toolbox (for bayesopt).
%
%   Usage:  test_optimize_rib_dryrun

if exist('bayesopt', 'file') == 0 && exist('bayesopt', 'builtin') == 0
    error('test_optimize_rib_dryrun:noBayesopt', ...
        'bayesopt is unavailable, so this test cannot run.');
end

rid = ['dryrun_test_', char(datetime('now', 'Format', 'yyyyMMddHHmmssSSS'))];
fprintf('=== optimize_rib_oblong_maxdef dry-run test ===\nrun id: %s\n\n', rid);

%% ---------- 1. cold Bayesian study ----------
[ob1, md1, fit1, ds1, O1, s1] = optimize_rib_oblong_maxdef( ...
    struct('dryRun', true, 'runId', rid, 'bo', struct('maxEvals', 10)));
store = OptimStore(O1.dbPath, O1.storeBackend);
T1 = store.loadEvals(rid);
assert(height(T1) == 10 && s1.nSolved == 10 && s1.nCacheHits == 0, ...
    'a cold study must solve exactly its budget');
assert(s1.nInfeasible == 0, 'XConstraintFcn let an unfabricable candidate through');
assert(abs(T1.oblong(1) - (-2)) < 1e-9 && abs(T1.maxdef(1) - 0.2) < 1e-12, ...
    'the start design must be the first evaluation');
assert(all(strcmp(T1.status, 'ok')), 'every dry-run evaluation should be ok');
assert(all(abs(T1.fitness - T1.fitnessRaw) <= 1e-12*T1.fitnessRaw), ...
    'with the gate off, fitness and fitnessRaw must agree');
assert(all(isfinite(T1.Qmech) & T1.Qmech > 0), ...
    'the Q_mech estimate must be recorded even with the gate off');
assert(isfinite(fit1) && abs(fit1 - max(T1.fitness)) <= 1e-9*fit1, ...
    'the reported best fitness must be the best stored one');
assert(isfinite(ob1) && md1 >= O1.maxdef_min && md1 <= O1.maxdef_max && isstruct(ds1), ...
    'best point / dsBest malformed');
assert(exist([O1.rootLoc, 'optimize_log.csv'], 'file') == 2 && ...
       exist([O1.rootLoc, 'optimize_convergence.png'], 'file') == 2, ...
    'CSV log or convergence figure missing');
% Coupling breakdown: the stub has gMB = 1.4 g and gPE = -0.4 g.
assert(all(abs(T1.gMB + T1.gPE - T1.gOM) <= 1e-9*T1.gOM) && all(T1.gPE < 0), ...
    'stored gMB + gPE must equal the stored net gOM, parts signed');
[~, k1] = max(T1.fitness);
c1 = s1.bestCpl;
assert(abs(c1.gOM - T1.gOM(k1)) <= 1e-9*c1.gOM && ...
       abs(c1.gMB - T1.gMB(k1)) <= 1e-9*c1.gOM && ...
       abs(c1.gPE - T1.gPE(k1)) <= 1e-9*c1.gOM, ...
    'stats.bestCpl must be the breakdown of the best stored row');
assert(abs(c1.cancelRatio - 1.8) < 1e-9, 'cancelRatio must be (|gMB|+|gPE|)/gOM');
assert(abs(c1.rerun.gMB - c1.gMB) <= 1e-9*c1.gOM && ...
       abs(c1.rerun.gPE - c1.gPE) <= 1e-9*c1.gOM, ...
    'the BEST re-run breakdown must reproduce the search-time one in a dry run');
hit1 = store.lookup(rid, [T1.dAR(k1), T1.maxdef(k1)]);
assert(abs(hit1.gMB - T1.gMB(k1)) <= 1e-9*c1.gOM && ...
       abs(hit1.gPE - T1.gPE(k1)) <= 1e-9*c1.gOM, ...
    'lookup must return the stored breakdown');
csvHead = strtrim(fgetl_first([O1.rootLoc, 'optimize_log.csv']));
assert(contains(csvHead, 'gOM_kHz,gMB_kHz,gPE_kHz'), 'CSV lacks the breakdown columns');
fprintf('\n[1] cold study : 10 solved, best fitness %.4g\n\n', fit1);

%% ---------- 2. resume with a larger budget ----------
[~, ~, fit2, ~, ~, s2] = optimize_rib_oblong_maxdef( ...
    struct('dryRun', true, 'runId', rid, 'bo', struct('maxEvals', 14)));
T2 = store.loadEvals(rid);
assert(s2.nSolved + s2.nCacheHits == 4, ...
    'a resumed study must spend only the budget the store does not cover');
assert(height(T2) == 10 + s2.nSolved, 'store row count disagrees with the solves');
assert(fit2 >= fit1*(1 - 1e-12), 'a resumed study reported a worse incumbent');
assert(isfinite(s2.bestCpl.gMB) && isfinite(s2.bestCpl.gPE) && ...
       abs(s2.bestCpl.gMB + s2.bestCpl.gPE - s2.bestCpl.gOM) <= 1e-9*s2.bestCpl.gOM, ...
    'a resumed study lost the coupling breakdown of its incumbent');
fprintf('\n[2] resume     : %d new solve(s), best fitness %.4g\n\n', s2.nSolved, fit2);

%% ---------- 3. cross-optimizer: Nelder-Mead on the same run ----------
[~, ~, fit3, ~, ~, s3] = optimize_rib_oblong_maxdef( ...
    struct('dryRun', true, 'runId', rid, 'optimizer', 'neldermead', ...
           'nm', struct('MaxFunEvals', 25)));
assert(s3.nCacheHits >= 1, ...
    'Nelder-Mead should have been served its start point from the store');
assert(s3.nSolved >= 1, 'Nelder-Mead made no new evaluation');
assert(fit3 >= fit2*(1 - 1e-12), 'the Nelder-Mead session reported a worse incumbent');
assert(isfinite(s3.bestCpl.gMB) && isfinite(s3.bestCpl.gPE), ...
    'the Nelder-Mead session lost the coupling breakdown of its incumbent');
fprintf('\n[3] neldermead : %d cached, %d solved, best fitness %.4g\n\n', ...
    s3.nCacheHits, s3.nSolved, fit3);

%% ---------- 4. the optional Q_mech gate ----------
rid4 = [rid, '_gate'];
[~, ~, ~, ~, O4] = optimize_rib_oblong_maxdef( ...
    struct('dryRun', true, 'runId', rid4, 'useQmechGate', 1, ...
           'bo', struct('maxEvals', 14)));
T4 = store.loadEvals(rid4);
gated = T4.fitness == 0;
fprintf('\n[4] gate on    : %d of %d design(s) gated\n\n', nnz(gated), height(T4));
if any(gated)
    assert(all(T4.Qmech(gated) < O4.QmechMin), ...
        'a row was gated despite passing the Q_mech threshold');
    assert(all(isfinite(T4.fitnessRaw(gated)) & T4.fitnessRaw(gated) > 0), ...
        'gated rows lost their ungated fitness');
    assert(all(strcmp(T4.status(gated), 'Qmech_gated')), ...
        'gated rows must carry the status Qmech_gated');
end
passed = ~gated & isfinite(T4.fitness);
assert(all(T4.Qmech(passed) >= O4.QmechMin), ...
    'a row passed the gate despite failing the Q_mech threshold');

%% ---------- 5. a mistyped option is an error ----------
threw = false;
try
    optimize_rib_oblong_maxdef(struct('dryRun', true, 'optimiser', 'bayesopt'));
catch ME
    threw = strcmp(ME.identifier, 'optimize_rib_oblong_maxdef:unknownOption');
end
assert(threw, 'an unknown option name must raise unknownOption');
fprintf('[5] unknown option rejected\n');

fprintf('\nALL RIB DRY-RUN ASSERTIONS PASSED\n');
end

function line = fgetl_first(fileName)
%FGETL_FIRST  First line of a text file.
fidIn = fopen(fileName, 'r');
closer = onCleanup(@() fclose(fidIn));
line = fgetl(fidIn);
end
