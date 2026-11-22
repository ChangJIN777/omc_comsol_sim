%% BAYESOPT_RIB
% Bayesian optimization of the 'rib' unit cell for a mechanical bandgap.
% bayesopt counterpart of optimize_rib_unitCell.m (fminsearch); the loop is
% bayesoptUnitCell.m, the rib-specific parts are in unitCellAdapter_rib.m.
%
% USAGE
%   1. cfg.solverBackend = 'surrogate';  run. The whole loop executes in
%      seconds with no COMSOL; every number is FABRICATED.
%   2. cfg.solverBackend = 'comsol';     run the real study (COMSOL with
%      LiveLink running; one sector solve per evaluation).
%   3. Extend:  S = load(<state file>);
%               results = resume(S.results, 'MaxObjectiveEvaluations', 20);
%      NOT across a change of bounds, searched variables or fabrication
%      limits: start a fresh study (optionally cfg.seedFromState).
%
% Run from omc-comsol-chang/. Needs the Statistics and Machine Learning Toolbox.
% See also BAYESOPTUNITCELL, UNITCELLADAPTER_RIB, OPTIMIZE_RIB_UNITCELL.

clear all; clc; close all; %#ok<CLALL>

A   = unitCellAdapter_rib();
cfg = A.defaultCfg();           % every default and its source is in the adapter

%% ===================== CONFIGURATION =====================
cfg.solverBackend = 'surrogate';
% cfg.solverBackend = 'comsol';

% --- design variables: bounds in INTEGER nm (d: integer degrees) -----------
% Sources: optimize_rib_unitCell.m:47-60.
cfg.vars.t.search = true;    cfg.vars.t.bounds = [ 90 200];
cfg.vars.s.search = true;    cfg.vars.s.bounds = [ 50 200];
cfg.vars.L.search = true;    cfg.vars.L.bounds = [250 500];
cfg.vars.a.search = false;   cfg.vars.a.fixed  = 300e-9;    % TODO(confirm): search a? [240 360]
cfg.vars.d.search = false;   cfg.vars.d.fixed  = 0;         % TODO(confirm): search d? [0 15] deg
cfg.vars.th.search = false;  cfg.vars.th.fixed = 110e-9;    % effective; nominal 120e-9

cfg.fixFootprint = false;    % true: L derived from the footprint area; set
                             % cfg.vars.L.search = false as well

% --- objective ------------------------------------------------------------
cfg.targetFreq = 11.2e9;     % Hz, optimize_rib_unitCell.m:72
cfg.sigma      = 1.0e9;      % Hz, TODO(confirm)
cfg.gapField   = 'sym';      % y-even & z-even sector (y-even only with a sidewall)
cfg.gapSelect  = 'bestScore';% TODO(confirm): 'closest' = optimize_rib rule
cfg.rxtal      = 0;          % TODO(confirm): 0 (<100>) or 45 (<110>)

% --- budget -----------------------------------------------------------------
cfg.maxEvaluations = 60;     % time ONE solve (test_rib_unitCell.m) first
cfg.numSeedPoints  = 12;
cfg.runTag = 'bayesopt_rib_run1';

%% ===================== RUN =====================
results = bayesoptUnitCell(cfg, A);
