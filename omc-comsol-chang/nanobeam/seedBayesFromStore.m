function [initX, initObj, initCon] = seedBayesFromStore(store, OPT)
%SEEDBAYESFROMSTORE  Turn stored evaluations into bayesopt seed data.
%
%   Returns a table of design points plus the matching objective and
%   coupled-constraint values, in the form bayesopt's 'InitialX' /
%   'InitialObjective' / 'InitialConstraintViolations' expect.
%
%   This is strictly better than a .mat checkpoint of the study object: it
%   survives a MATLAB version change, it is queryable with plain SQL, and it
%   lets a Nelder-Mead run seed a Bayesian run -- which is the comparison
%   worth doing, since both write the same physics to the same table.
initX   = table(zeros(0,1), zeros(0,1), 'VariableNames', {'dAR','maxdef'});
initObj = zeros(0,1);
initCon = zeros(0,2);

if isempty(store) || ~OPT.useStore || ~OPT.resume
    return;
end

T = store.loadEvals(OPT.runId);
if isempty(T) || height(T) == 0
    return;
end

% Rows written before fitness_raw existed cannot seed the GP: their ungated
% score is unknown, and inventing one would poison the surrogate.  They stay
% in the database and still serve the Nelder-Mead cache.
usable = isfinite(T.fitnessRaw) & T.fitnessRaw > 0;
if ~any(usable)
    return;
end

initX   = table(T.dAR(usable), T.maxdef(usable), ...
                'VariableNames', {'dAR','maxdef'});
initObj = -log10(T.fitnessRaw(usable));

qm   = T.Qmech(usable);
con1 = 1 - qm./OPT.QmechMin;
con1(~isfinite(con1)) = 1;
con2 = -ones(size(con1));   % all of these rows have a usable fitness
initCon = [con1, con2];
end
