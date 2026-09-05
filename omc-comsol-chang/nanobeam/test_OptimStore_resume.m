function test_OptimStore_resume()
%TEST_OPTIMSTORE_RESUME  Verify stop/resume works, with no COMSOL needed.
%
%   Exercises exactly the plumbing that optimize_oblong_maxdef.m uses for
%   persistence -- OptimStore + a memoized fminsearch objective -- but with a
%   cheap analytic fitness standing in for RunNanobeamFEM, so the whole
%   stop/resume cycle runs in a second instead of days.
%
%   What it asserts:
%     1. A first session persists every evaluation it makes.
%     2. A second session with the same runId replays that path entirely from
%        the database -- zero re-solves -- and then extends the search.
%     3. The replayed path is IDENTICAL to the original, which is what makes
%        the fast-forward sound: fminsearch is deterministic given x0.
%     4. A fresh runId in the same database shares nothing with the first.
%     5. The incumbent survives across sessions.
%
%   Usage:  test_OptimStore_resume

close all;

dbPath = fullfile(tempdir, ['optimstore_resume_', ...
    char(datetime('now','Format','yyyyMMddHHmmssSSS')), '.sqlite3']);
fprintf('=== OptimStore resume test ===\ndb: %s\n\n', dbPath);

% Analytic stand-in for the FEM fitness: a smooth peak inside the same kind of
% box optimize_oblong_maxdef searches, with the same "higher is better,
% objective is -log10(fitness)" convention.
fitnessFcn = @(x) 1e16 * exp(-((x(1)-1.18)/0.25)^2 - ((x(2)-0.207)/0.05)^2);

x0    = [1.1135, 0.22];
runA  = 'resume_test_A';
runB  = 'resume_test_B';

%% ---------- session 1: a short run that we then "interrupt" ----------
[nSolve1, nCache1, path1, best1] = runSession(dbPath, runA, x0, fitnessFcn, 12);
fprintf('session 1 : %d solved, %d cached, best fitness %.6g\n', ...
    nSolve1, nCache1, best1);
assert(nSolve1 > 0, 'session 1 solved nothing');
assert(nCache1 == 0, 'session 1 should have had an empty database');

%% ---------- session 2: same runId, resume, extend the budget ----------
[nSolve2, nCache2, path2, best2] = runSession(dbPath, runA, x0, fitnessFcn, 30);
fprintf('session 2 : %d solved, %d cached, best fitness %.6g\n', ...
    nSolve2, nCache2, best2);

% The replay must cover everything session 1 paid for.
assert(nCache2 >= nSolve1, ...
    'resume re-solved points it should have served from the database (%d cached < %d stored)', ...
    nCache2, nSolve1);

% And the retraced prefix must match step for step -- if it does not, the
% fast-forward is serving answers for a path the optimizer is not on.
nCommon = min(size(path1,1), size(path2,1));
assert(nCommon >= size(path1,1), 'session 2 was shorter than session 1');
assert(max(max(abs(path1 - path2(1:size(path1,1), :)))) < 1e-12, ...
    'resumed path diverged from the original path');

% A bigger budget must not make the answer worse.
assert(best2 >= best1 - 1e-6*abs(best1), ...
    'resuming lost ground: best went from %.6g to %.6g', best1, best2);

fprintf('           replayed %d point(s) identically, then extended\n', nCache2);

%% ---------- session 3: a different runId shares nothing ----------
[nSolve3, nCache3] = runSession(dbPath, runB, x0, fitnessFcn, 8);
fprintf('session 3 : %d solved, %d cached (fresh runId)\n', nSolve3, nCache3);
assert(nCache3 == 0, 'a fresh runId reused another run''s evaluations');

%% ---------- the database agrees with what the sessions reported ----------
store = OptimStore(dbPath);
R = store.listRuns();
fprintf('\nruns in database:\n');
disp(R(:, {'run_id','optimizer','status','nEvals'}));
assert(height(R) == 2, 'expected exactly 2 runs, found %d', height(R));

TA = store.loadEvals(runA);
assert(height(TA) == store.evalCount(runA), 'loadEvals/evalCount disagree');
assert(all(diff(TA.eval) > 0), 'stored eval indices are not strictly increasing');
assert(~any(isnan(TA.fitness)), 'analytic fitness should never be NaN');

[xb, fb] = store.best(runA);
assert(abs(fb - best2) < 1e-6*abs(best2), 'store.best disagrees with the session');
fprintf('\nbest stored: dAR=%.6f maxdef=%.6f fitness=%.6g\n', xb(1), xb(2), fb);

fprintf('\nALL RESUME ASSERTIONS PASSED\n');
end

% =========================================================================
function [nSolve, nCache, pathXY, bestFit] = runSession(dbPath, runId, x0, ...
    fitnessFcn, maxFunEvals)
%RUNSESSION  One optimizer session against the store, mirroring the real script.

Jpenalty = 1e3;
bounds   = struct('dARmin', 0.5567, 'dARmax', 1.6703, ...
                  'mdMin',  0.11,   'mdMax',  0.33);

store  = OptimStore(dbPath);
store.startRun(runId, 'test_OptimStore_resume', 'neldermead', ...
    struct('maxFunEvals', maxFunEvals));

nSolve   = 0;
nCache   = 0;
pathXY   = zeros(0, 2);
bestFit  = -Inf;
evalCount = store.evalCount(runId);

opts = optimset('Display','off','TolX',1e-6,'TolFun',1e-6, ...
    'MaxIter',1000,'MaxFunEvals',maxFunEvals);
fminsearch(@objective, x0, opts);

store.finishRun(runId, 'complete');

    function J = objective(x)
        % Same sloped barrier as the real script.
        dARrange = bounds.dARmax - bounds.dARmin;
        mdRange  = bounds.mdMax  - bounds.mdMin;
        viol = max([(bounds.dARmin - x(1))/dARrange, ...
                    (x(1) - bounds.dARmax)/dARrange, ...
                    (bounds.mdMin - x(2))/mdRange, ...
                    (x(2) - bounds.mdMax)/mdRange, 0]);
        if viol > 0
            J = Jpenalty * (2 + viol);
            return;
        end

        pathXY(end+1, :) = x(:).';

        % --- the resume fast-forward under test ---
        hit = store.lookup(runId, x);
        if ~isempty(hit)
            J      = hit.J;
            nCache = nCache + 1;
            if hit.fitness > bestFit; bestFit = hit.fitness; end
            evalCount = max(evalCount, hit.eval_idx);
            return;
        end

        % --- the "expensive solve" ---
        evalCount = evalCount + 1;
        nSolve    = nSolve + 1;
        fit = fitnessFcn(x);
        J   = -log10(fit);
        if fit > bestFit; bestFit = fit; end

        store.logEval(runId, evalCount, struct('x', x, ...
            'oblong', log(x(1)/(578/196))/(2*log(1-x(2))), ...
            'fitness', fit, 'J', J, 'status', 'ok'));
    end
end
