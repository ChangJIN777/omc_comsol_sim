function test_bayesopt_wiring()
%TEST_BAYESOPT_WIRING  Verify the bayesopt backend plumbing, with no COMSOL.
%
%   Exercises exactly the machinery optimize_oblong_maxdef.m uses when
%   OPT.optimizer = 'bayesopt' -- the table-in/objective+constraints-out
%   adapter, the Q_mech coupled constraint, OptimStore persistence of BOTH the
%   gated and ungated fitness, and seeding a resumed study from the store --
%   but with a cheap analytic fitness standing in for RunNanobeamFEM.
%
%   What it asserts:
%     1. bayesopt runs to completion against the constrained adapter.
%     2. Every evaluation is persisted with both `fitness` and `fitnessRaw`.
%     3. Gated designs (Q_mech < QmechMin) get fitness == 0 but keep a finite
%        ungated fitness -- which is the whole point of storing both.
%     4. The Q_mech constraint is positive exactly on the gated designs.
%     5. seedBayesFromStore returns well-formed InitialX/Objective/Constraints
%        that bayesopt accepts, and a seeded study starts from them.
%     6. A Nelder-Mead-style cache read can re-derive its own penalty-ladder J
%        from rows a Bayesian run wrote -- the cross-optimizer resume path.
%
%   Requires the Statistics and Machine Learning Toolbox (for bayesopt).
%
%   Usage:  test_bayesopt_wiring

close all;

if exist('bayesopt', 'file') == 0 && exist('bayesopt', 'builtin') == 0
    error('test_bayesopt_wiring:noBayesopt', ...
        ['bayesopt is unavailable, so this test cannot run.\n' ...
         '  licence permits it (1=yes) : %d\n' ...
         '  installed on disk (1=yes)  : %d'], ...
        license('test','Statistics_Toolbox'), ...
        isfolder(fullfile(matlabroot,'toolbox','stats')));
end

dbPath = fullfile(tempdir, ['bayeswiring_', ...
    char(datetime('now','Format','yyyyMMddHHmmssSSS')), '.sqlite3']);
fprintf('=== bayesopt wiring test ===\ndb: %s\n\n', dbPath);

OPT = struct();
OPT.defectAspectRatio_min = 0.5567;
OPT.defectAspectRatio_max = 1.6703;
OPT.maxdef_min = 0.11;
OPT.maxdef_max = 0.33;
OPT.QmechMin   = 1e7;
OPT.Jpenalty   = 1e3;
OPT.useStore   = 1;
OPT.resume     = 1;
OPT.runId      = 'bo_wiring_A';

%% ---------- session 1: a cold Bayesian study ----------
[nEval1, res1] = runBO(dbPath, OPT, 14, 0);
fprintf('session 1 : %d evaluation(s), %d objective row(s) in the study\n', ...
    nEval1, height(res1.XTrace));
assert(nEval1 > 0, 'session 1 evaluated nothing');

store = OptimStore(dbPath);
T = store.loadEvals(OPT.runId);
assert(height(T) == nEval1, 'store row count disagrees with the objective calls');

%% ---------- both fitness flavours are persisted ----------
assert(any(isfinite(T.fitnessRaw)), 'no ungated fitness was stored');
gated = T.fitness == 0;
fprintf('           %d of %d design(s) gated by Q_mech\n', nnz(gated), height(T));
if any(gated)
    % This is the property the whole design rests on: a gated design is
    % worthless to Nelder-Mead but still informative to the GP.
    assert(all(isfinite(T.fitnessRaw(gated)) & T.fitnessRaw(gated) > 0), ...
        'gated rows lost their ungated fitness');
    assert(all(T.Qmech(gated) < OPT.QmechMin), ...
        'a row was gated despite passing the Q_mech threshold');
end
ungated = ~gated & isfinite(T.fitness);
if any(ungated)
    assert(all(T.Qmech(ungated) >= OPT.QmechMin), ...
        'a row passed the gate despite failing the Q_mech threshold');
    assert(all(abs(T.fitness(ungated) - T.fitnessRaw(ungated)) < ...
               1e-6*abs(T.fitnessRaw(ungated))), ...
        'above the gate the two fitness values must agree');
end

%% ---------- the seed helper produces what bayesopt expects ----------
[initX, initObj, initCon] = seedBayesFromStore(store, OPT);
fprintf('           seed: %d point(s), %d objective(s), %dx%d constraint(s)\n', ...
    height(initX), numel(initObj), size(initCon,1), size(initCon,2));
assert(height(initX) == numel(initObj), 'seed X/objective length mismatch');
assert(size(initCon,1) == height(initX) && size(initCon,2) == 2, ...
    'seed constraint matrix has the wrong shape');
assert(all(strcmp(initX.Properties.VariableNames, {'dAR','maxdef'})), ...
    'seed table variable names must match the optimizableVariables');
assert(all(isfinite(initObj)), 'seed objectives must all be finite');
% Constraint 1 is positive exactly where the gate is violated.
qmSeed = T.Qmech(isfinite(T.fitnessRaw) & T.fitnessRaw > 0);
expectCon1 = 1 - qmSeed./OPT.QmechMin;
assert(max(abs(initCon(:,1) - expectCon1)) < 1e-9, ...
    'seed Q_mech constraint does not match 1 - Qmech/QmechMin');
assert(all(initCon(:,2) == -1), 'usable-solve constraint should be satisfied');

%% ---------- session 2: resumed, seeded from the store ----------
% bayesopt counts seeded points against MaxObjectiveEvaluations, so a
% resumed study's budget must exceed what the store already holds.
[nEval2, res2] = runBO(dbPath, OPT, height(initX) + 6, height(initX));
fprintf('session 2 : %d new evaluation(s), study holds %d point(s)\n', ...
    nEval2, height(res2.XTrace));
assert(height(res2.XTrace) > height(initX), ...
    'the seeded study did not extend beyond its seed points');

%% ---------- cross-optimizer: Nelder-Mead can read these rows ----------
% The Bayesian run wrote objective_j on ITS scale; a Nelder-Mead resume must
% re-derive its own penalty-ladder J from the stored physics instead.
hit = store.lookup(OPT.runId, [T.dAR(1), T.maxdef(1)]);
assert(~isempty(hit), 'lookup missed a row the Bayesian run wrote');
assert(isfield(hit,'fitnessRaw') && isfield(hit,'Qmech'), ...
    'lookup must expose fitnessRaw and Qmech for cross-optimizer resume');
if isfinite(hit.fitness) && hit.fitness > 0
    Jnm = -log10(hit.fitness);
elseif isfinite(hit.fitness)
    Jnm = OPT.Jpenalty;
else
    Jnm = 2*OPT.Jpenalty;
end
fprintf('           NM re-derives J = %.4g from a BO-written row\n', Jnm);
assert(isfinite(Jnm), 'could not re-derive a Nelder-Mead objective');

fprintf('\nALL BAYESOPT WIRING ASSERTIONS PASSED\n');
end

% =========================================================================
function [nEval, results] = runBO(dbPath, OPT, maxEvals, nSeeded)
%RUNBO  One Bayesian session against the store, mirroring the real dispatch.

store = OptimStore(dbPath);
store.startRun(OPT.runId, 'test_bayesopt_wiring', 'bayesopt', OPT);

nEval     = 0;
evalCount = store.evalCount(OPT.runId);

optVars = [ ...
    optimizableVariable('dAR', ...
        [OPT.defectAspectRatio_min, OPT.defectAspectRatio_max], 'Type','real'), ...
    optimizableVariable('maxdef', ...
        [OPT.maxdef_min, OPT.maxdef_max], 'Type','real')];

[initX, initObj, initCon] = seedBayesFromStore(store, OPT);

boArgs = {'MaxObjectiveEvaluations', maxEvals, ...
          'NumSeedPoints', max(1, 4 - nSeeded), ...
          'AcquisitionFunctionName', 'expected-improvement-plus', ...
          'IsObjectiveDeterministic', false, ...
          'NumCoupledConstraints', 2, ...
          'AreCoupledConstraintsDeterministic', [false false], ...
          'PlotFcn', [], 'Verbose', 0};
if height(initX) > 0
    boArgs = [boArgs, {'InitialX', initX, 'InitialObjective', initObj, ...
                       'InitialConstraintViolations', initCon}];
end

results = bayesopt(@boObjective, optVars, boArgs{:});
store.finishRun(OPT.runId, 'complete');

    function [obj, con] = boObjective(t)
        x = [t.dAR, t.maxdef];

        % --- serve from the store if already solved (the resume path) ---
        hit = store.lookup(OPT.runId, x);
        if ~isempty(hit)
            fitRaw = hit.fitnessRaw;
            qm     = hit.Qmech;
        else
            % --- the "expensive solve", analytic here ---
            evalCount = evalCount + 1;
            nEval     = nEval + 1;

            % Smooth peak, plus a Q_mech field that fails the gate in one
            % corner so the constraint actually has something to learn.
            fitRaw = 1e16 * exp(-((x(1)-1.18)/0.25)^2 - ((x(2)-0.207)/0.05)^2);
            qm     = 3e7 * exp(-((x(2)-0.19)/0.09)^2);

            if qm < OPT.QmechMin
                fit = 0;
                status = 'Qmech_gated';
            else
                fit = fitRaw;
                status = 'ok';
            end

            store.logEval(OPT.runId, evalCount, struct('x', x, ...
                'oblong', log(x(1)/(578/196))/(2*log(1-x(2))), ...
                'Qmech', qm, 'fitness', fit, 'fitnessRaw', fitRaw, ...
                'J', -log10(max(fitRaw, realmin)), 'status', status));
        end

        % --- the adapter under test ---
        if isfinite(fitRaw) && fitRaw > 0
            obj = -log10(fitRaw);
        else
            obj = NaN;
        end
        if isfinite(qm)
            con1 = 1 - qm/OPT.QmechMin;
        else
            con1 = 1;
        end
        con2 = -1;
        if ~(isfinite(fitRaw) && fitRaw > 0)
            con2 = 1;
        end
        con = [con1, con2];
    end
end
