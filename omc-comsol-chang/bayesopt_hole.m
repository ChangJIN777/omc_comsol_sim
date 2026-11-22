%% BAYESOPT_HOLE
% Bayesian optimization of the 'hole' unit cell for a mechanical bandgap.
% bayesopt counterpart of optimize_hole_unitCell.m / sweep_hole_unitCell.m; the
% loop is bayesoptUnitCell.m, the hole-specific parts are in
% unitCellAdapter_hole.m.
%
% USAGE: as bayesopt_rib.m -- 'surrogate' first, then 'comsol'.
%
% Run from omc-comsol-chang/. Needs the Statistics and Machine Learning Toolbox.
% See also BAYESOPTUNITCELL, UNITCELLADAPTER_HOLE, OPTIMIZE_HOLE_UNITCELL.

clear all; clc; close all; %#ok<CLALL>

A   = unitCellAdapter_hole();
cfg = A.defaultCfg();           % every default and its source is in the adapter

%% ===================== CONFIGURATION =====================
cfg.solverBackend = 'surrogate';
% cfg.solverBackend = 'comsol';

% --- design variables: bounds in INTEGER nm ---------------------------------
cfg.vars.hx.search = true;   cfg.vars.hx.bounds = [100 500];   % TODO(confirm) box
cfg.vars.hy.search = true;   cfg.vars.hy.bounds = [150 600];   % TODO(confirm) box
cfg.vars.a.search  = false;  cfg.vars.a.fixed   = 618e-9;      % TODO(confirm) baseline
cfg.vars.bw.search = false;  cfg.vars.bw.fixed  = 600e-9;      %   (optimize_hole_unitCell.m:65-71
cfg.vars.th.search = false;  cfg.vars.th.fixed  = 380e-9;      %    vs sweep/test: 528/800/500)

% Fixed ellipse area, as optimize_hole_unitCell.m: 1-D search over hx with
% hy = A0/hx. To use it: cfg.fixArea = true; cfg.vars.hy.search = false;
% and optionally cfg.A0 = 270e-9^2; cfg.aspectRange = [0.3 1.5];
cfg.fixArea = false;

% --- objective ------------------------------------------------------------
cfg.targetFreq = 8e9;        % Hz
cfg.sigma      = 0.5e9;      % Hz, detuning tolerance
% fitness = gapSize[GHz] * exp(-((targetFreq - midGap)/sigma)^2), absolute gap
cfg.gapField   = 'sym';      % y-even sector
cfg.gapSelect  = 'bestScore';% gap with the highest fitness ('closest' = optimize_hole rule)
cfg.rxtal      = 0;          % 0 for now (optimize_hole_unitCell.m:97 uses 45)

% --- budget -----------------------------------------------------------------
cfg.maxEvaluations = 60;     % time ONE solve (test_hole_unitCell.m) first
cfg.numSeedPoints  = 12;
cfg.runTag = 'bayesopt_hole_run1';

%% ===================== RUN =====================
results = bayesoptUnitCell(cfg, A);
