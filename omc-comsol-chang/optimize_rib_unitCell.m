function [xBest, fitnessBest, dsBest, PBest, OPT] = optimize_rib_unitCell(optOverrides, pOverrides)
%OPTIMIZE_RIB_UNITCELL Optimize the 'rib' unit-cell geometry for a target mechanical bandgap.
%   Nelder-Mead (fminsearch) search over a configurable set of rib design
%   variables (default t, s, L; a and the sidewall angle d are optional) that
%   maximizes a user-defined fitness on the band structure returned by
%   solveBands. Each evaluation writes to a unique datLoc folder, so solveBands
%   never reuses a cached result. Candidates that fail validateRibParams are
%   scored with a penalty, and COMSOL is never started for them.
%
%   [xBest, fitnessBest, dsBest, PBest, OPT] = optimize_rib_unitCell()
%   [...] = optimize_rib_unitCell(optOverrides)              % e.g. struct('dryRun',true)
%   [...] = optimize_rib_unitCell(optOverrides, pOverrides)  % also override P fields
%
%   xBest       struct of the best design variables (meters; d in radians)
%   fitnessBest best fitness (units of the chosen fitness function)
%   dsBest      band-structure struct from a final plotted re-run (in BEST_...)
%   PBest       the P struct of that re-run
%   OPT         the options actually used (after overrides)
%
%   optOverrides / pOverrides are structs merged over the knobs / P defaults
%   below. Nested structs (OPT.init, OPT.lb, OPT.ub) merge field by field.
%   This is how RibUnitCellTest runs the dry run in a temporary folder.
%
%   Based on optimize_hole_unitCell.m. The rib device and its parameters are
%   described in rib_README.md (section "Optimization").
if nargin < 1 || isempty(optOverrides), optOverrides = struct(); end
if nargin < 2 || isempty(pOverrides),   pOverrides   = struct(); end

%% ===================== TUNABLE KNOBS =====================
% --- dry run: analytic surrogate instead of solveBands (no COMSOL) ---
% Exercises the plumbing only (normalization, clamping, feasibility penalty,
% CSV log, BEST re-run path). The surrogate is NOT physics.
OPT.dryRun = false;

% --- design variables -----------------------------------------------------
% Any subset of {'t','s','L','a','d'}, in meters (d in radians):
%   t  rib width along x (bottom face)
%   s  spine width along y (bottom face)
%   L  rib half-length, beam axis to tip. Set explicitly as P.L, which wins
%      over the L(a) rule (see ribLength.m)
%   a  lattice constant
%   d  sidewall angle. Needs P.r1 = P.r2 = 0 (fillets cannot be combined
%      with d ~= 0) and switches the solve to y-symmetry only (see below)
% Variables that are not listed keep their P value.
OPT.vars = {'t','s','L'};

% Initial values and [lb, ub] box per variable. Sources:
%   t, s : 100 nm = ribDefaultParams.m, t0/s0 in rib_optimize_sweep_diamond.m,
%          and every logged run in test/diamond_ribUnitCell*/ (t, s in 100-105 nm)
%   L    : 370 nm = the L(a) rule at a = 300 nm (ribDefaultParams). The logged
%          diamond runs used w = 2L = 600-840 nm (L = 300-420 nm)
%   a    : 300 nm = ribDefaultParams (template mirror pitch). The logged runs
%          used 250-262 nm
%   d    : 0 = vertical walls (visible design); LN rib scripts use 15 deg
% Bounds: t >= 90 nm keeps 2*r2 < t with the default 40 nm fillets, and
% t <= 200 nm keeps r1 < (a - t)/2 at a = 300 nm. s >= 50 nm is the
% fabrication floor used in rib_optimize_sweep_diamond.m (min(a-t, s, th) >= 50 nm).
OPT.init = struct('t',100e-9, 's',100e-9, 'L',370e-9, 'a',300e-9, 'd',0);
OPT.lb   = struct('t', 90e-9, 's', 50e-9, 'L',250e-9, 'a',240e-9, 'd',0);
OPT.ub   = struct('t',200e-9, 's',200e-9, 'L',500e-9, 'a',360e-9, 'd',15*pi/180);

% --- optional fixed-quantity constraint (off by default) ---------------------
% Analogue of the hole optimizer's fixed area. Hold the cell's footprint
% area A = a*s + t*(2L - s) (full y) at its initial value by solving for one
% dependent variable, which must NOT be in OPT.vars.
OPT.fixFootprint = false;
OPT.fixDependent = 'L';           % currently only 'L' is supported

% --- target mechanical frequency (Hz) the gap should be centred on ---
% Visible design: the template's mechanical eigen-shift is 11.2 GHz
% (comsol_templates/visible_rib_omc.m:23488, eigen-unit GHz at :1262).
OPT.targetFreq = 11.2e9;

% --- penalty (in the MINIMIZED objective) for no usable gap / infeasible ---
% Must exceed any realistic fitness. Infeasible candidates get
% noGapPenalty*(1 + violation), so the penalty grows with the distance
% outside the feasible region and fminsearch is pushed back.
% Out-of-box proposals are clamped and additionally penalized by
% boxPenaltySlope per unit of normalized distance outside the box.
OPT.noGapPenalty    = 1e9;
OPT.boxPenaltySlope = 1e9;

% --- output root: unique per run (time-stamped), one sub-folder per evaluation ---
% The time stamp matters: if two runs on the same day shared a root, the
% second run's eval_0001... folders would hit the first run's cached
% results, and solveBands would return empty bands.
OPT.rootDir = fullfile('.','test','ribUnitCell_opt', ['opt_', datestr(now,'mmddyyyy_HHMMSS')]);

% --- fminsearch options (TolX on the normalized variables, TolFun in fitness units) ---
OPT.fminOpts = optimset('Display','iter','TolX',1e-3,'TolFun',1e6, ...
    'MaxIter',80,'MaxFunEvals',160);

% --- FITNESS FUNCTION (user-configurable) ---------------------------------
% Maps the solveBands output struct ds -> scalar score (HIGHER is better).
% Return NaN/[] when there is no usable gap (the objective then applies
% noGapPenalty). Options (local helpers at the end of this file):
%   fitness_closestGap_sym       (default) size of the SYMMETRIC gap nearest
%                                targetFreq; uses ds.sym. Needs
%                                P.TwoSymPlanes = 0 (see P defaults)
%   fitness_closestGapRatio_sym  gap-to-midgap ratio of that symmetric gap
%                                (dimensionless; also set TolFun ~1e-3)
%   fitness_closestGap           size of the COMPLETE gap nearest targetFreq;
%                                uses ds.full, which is complete only over
%                                the sectors solved: set P.TwoSymPlanes = 1
%                                (all four y/z sectors, d = 0 only) or
%                                P.solveasym = 1
% A handle @(ds, targetFreq) or one of the names 'closestGap_sym',
% 'closestGapRatio_sym', 'closestGap' (names also work through optOverrides).
OPT.fitnessFcn = 'closestGap_sym';

OPT = mergeStruct(OPT, optOverrides);
if ischar(OPT.fitnessFcn) || isstring(OPT.fitnessFcn)
    switch char(OPT.fitnessFcn)
        case 'closestGap_sym',      OPT.fitnessFcn = @fitness_closestGap_sym;
        case 'closestGapRatio_sym', OPT.fitnessFcn = @fitness_closestGapRatio_sym;
        case 'closestGap',          OPT.fitnessFcn = @fitness_closestGap;
        otherwise
            error('optimize_rib_unitCell:badFitness', ...
                'Unknown OPT.fitnessFcn ''%s''.', char(OPT.fitnessFcn));
    end
end

%% ===================== P STRUCT DEFAULTS =====================
% Geometry and material from ribDefaultParams (the single defaults source);
% only the solver / IO fields are set here, following test_rib_unitCell.m.
P = ribDefaultParams();
P.celltype = 'rib';
% TODO(confirm): crystal orientation. 0 = beam along <100> (template and
% cavity); 45 = beam along <110>. The gap shifts by several percent between
% the two (in-plane shear velocity -9 %, longitudinal +5 %), so match this to
% the fabricated devices and to the cavity runs.
P.rxtal = 0;
P.anisoMat = 1;
P.d = 0;                    % vertical walls: fillets r1 = r2 = 40 nm allowed
P.kpts = 10;                % k-points, EXCLUDING gamma
P.nbands = 10;              % bands per symmetry sector
% Symmetry: with TwoSymPlanes = 0 solveBands solves ONE sector into ds.sym
% (plus ds.asym if solveasym = 1):
%   zSymCondition = 1 -> mbevenz = +1, mbeveny kept from P. With mbeveny = 1
%                        ds.sym is the y-even AND z-even sector (the OMC
%                        breathing modes) on a quarter model. d = 0 only.
%   zSymCondition = 0 -> mbeveny = +1, mbevenz kept (0): y-even, both z
%                        parities. Required when d ~= 0 (set automatically).
P.TwoSymPlanes  = 0;
P.zSymCondition = 1;
P.mbeveny = 1;
P.mbevenz = 0;
P.solveasym = 0;            % only the symmetric sector (one FEM solve per point)
P.completeBandGaps = 1;
P.plotgeom = 0;             % off during optimization
P.savedat = 1;
P.savebndplot = 0;          % off during optimization
P.saveplots = 0;
P.saveMPH = 0;
P.bandStruct_2D = 0;
P.run_optical = 0;          % mechanics only (runBands refuses the air domain)
P.freq = OPT.targetFreq;
P.meshSize = 3;
P.fixed_bc = 0;             % no flat outer y face exists to clamp
P.max_dof = 3e6;

P = mergeStruct(P, pOverrides);

% --- consistency of the knobs ---
nv = numel(OPT.vars);
allowed = {'t','s','L','a','d'};
if nv == 0 || ~all(ismember(OPT.vars, allowed))
    error('optimize_rib_unitCell:badVars', ...
        'OPT.vars must be a nonempty subset of {%s}.', strjoin(allowed, ', '));
end
usesSidewall = ismember('d', OPT.vars) || P.d > 0;
if usesSidewall
    if P.r1 > 0 || P.r2 > 0
        error('optimize_rib_unitCell:filletWithSidewall', ...
            ['A sidewall angle is optimized or set (d ~= 0), but P.r1/P.r2 are ' ...
             'nonzero. Set pOverrides.r1 = pOverrides.r2 = 0.']);
    end
    if P.zSymCondition || P.TwoSymPlanes
        fprintf(['Note: d ~= 0 is not mirror-symmetric in z; switching to ' ...
                 'zSymCondition = 0, TwoSymPlanes = 0, mbevenz = 0 (y-even sector).\n']);
        P.zSymCondition = 0;  P.TwoSymPlanes = 0;  P.mbevenz = 0;
    end
end
if OPT.fixFootprint
    if ~strcmp(OPT.fixDependent, 'L') || ismember(OPT.fixDependent, OPT.vars)
        error('optimize_rib_unitCell:badFixFootprint', ...
            'OPT.fixDependent must be ''L'' and must not be one of OPT.vars.');
    end
    x0 = OPT.init;
    OPT.A0 = x0.a*x0.s + x0.t*(2*x0.L - x0.s);   % fixed footprint area (m^2)
end
if ~OPT.dryRun
    close all;
end

%% ===================== OPTIMIZATION =====================
evalCounter = 0;        % shared with the nested objective below

% Normalized coordinates v in [1, 2] map linearly onto [lb, ub]. Offsetting by
% 1 (instead of using [0, 1]) makes fminsearch's default initial simplex,
% which perturbs each nonzero component by 5 %, a uniform 5-10 % step of
% every variable's range, so nm-scale variables with different ranges stay
% equally well conditioned.
lb = cellfun(@(f) OPT.lb.(f),   OPT.vars);
ub = cellfun(@(f) OPT.ub.(f),   OPT.vars);
x0 = cellfun(@(f) OPT.init.(f), OPT.vars);
if any(ub <= lb) || any(x0 < lb | x0 > ub)
    error('optimize_rib_unitCell:badBounds', ...
        'Each OPT.init must lie in [OPT.lb, OPT.ub] and lb < ub.');
end
toX = @(v) lb + (v - 1).*(ub - lb);
toV = @(x) 1 + (x - lb)./(ub - lb);

% CSV log of the trajectory, created fresh for each run.
if ~exist(OPT.rootDir, 'dir'); mkdir(OPT.rootDir); end
OPT.logFile = fullfile(OPT.rootDir, 'opt_log.csv');
colNames = cellfun(@varColName, OPT.vars, 'UniformOutput', false);
fid = fopen(OPT.logFile, 'w');
fprintf(fid, 'eval,%s,%s,feasible,errorID,boxDist,fitness,J\n', ...
    strjoin(colNames, ','), strjoin(strcat('v_', OPT.vars), ','));
fclose(fid);

[vBest, Jbest] = fminsearch(@objective, toV(x0), OPT.fminOpts);

% Clamp the returned optimum to the box (fminsearch is unconstrained).
xb = toX(min(max(vBest, 1), 2));
xBest = struct();
for k = 1:nv
    xBest.(OPT.vars{k}) = xb(k);
end
fitnessBest = -Jbest;   % negate: fminsearch minimizes, so Jbest = -(best fitness)
noFeasible = Jbest >= OPT.noGapPenalty;
if noFeasible
    fitnessBest = NaN;
    warning('optimize_rib_unitCell:noUsableDesign', ...
        ['No feasible design with a usable gap was found (best J = %.3g). ' ...
         'Check OPT.init / OPT.lb / OPT.ub and the log %s.'], Jbest, OPT.logFile);
end

fprintf('\n==================== OPTIMIZATION COMPLETE ====================\n');
for k = 1:nv
    fprintf('Best %-2s = %s\n', OPT.vars{k}, varString(OPT.vars{k}, xb(k)));
end
fprintf('Best fitness = %.4g\n', fitnessBest);
fprintf('Evaluations  = %d\n', evalCounter);
fprintf('Log          = %s\n', OPT.logFile);

%% ===================== FINAL RE-RUN WITH PLOTS =====================
PBest = applyDesign(P, OPT, xb);
if noFeasible || ~checkFeasible(PBest)
    % An infeasible geometry would only make the builder error out.
    fprintf('\nBest point is infeasible or gap-free; skipping the final re-run.\n');
    dsBest = [];
    return
end
PBest.plotgeom = ~OPT.dryRun;
PBest.savebndplot = ~OPT.dryRun;
PBest.datLoc = [fullfile(OPT.rootDir, 'BEST'), filesep];
if isfield(PBest,'fileBase'); PBest = rmfield(PBest,'fileBase'); end
fprintf('\nRe-running best design -> %s\n', PBest.datLoc);
if OPT.dryRun
    dsBest = surrogateBands(PBest);
else
    dsBest = solveBands(PBest);
    removeStrayBackslashDir(PBest.datLoc);
end

%% ===================== NESTED OBJECTIVE =====================
    function J = objective(v)
        % v = normalized design vector. Clamp to the box, map to P, check
        % feasibility (no COMSOL if infeasible), solve in a unique folder,
        % score with the user fitness, return -fitness (fminsearch minimizes).
        % evalCounter lives in the enclosing function's workspace; nested
        % functions share it directly, so no globals are needed.
        evalCounter = evalCounter + 1;

        % fminsearch can propose points outside the box. Clamp onto the box,
        % and penalize the distance so the simplex is pushed back inside.
        vc = min(max(v(:)', 1), 2);
        boxDist = norm(v(:)' - vc);
        x = toX(vc);
        Pe = applyDesign(P, OPT, x);
        Pe.plotgeom    = 0;
        Pe.savebndplot = 0;
        Pe.saveplots   = 0;

        % Feasibility before COMSOL: the same checks the builder runs.
        [feasible, errId, viol] = checkFeasible(Pe);
        fit = NaN;
        if ~feasible
            J = OPT.noGapPenalty*(1 + viol);
        else
            % Unique folder per evaluation: solveBands loads cached results and
            % skips the COMSOL solve if <datLoc>/<fileBase>_bds.mat exists.
            Pe.datLoc = [fullfile(OPT.rootDir, sprintf('eval_%04d', evalCounter)), filesep];
            % Strip fileBase so it is regenerated from the current geometry.
            if isfield(Pe,'fileBase'); Pe = rmfield(Pe,'fileBase'); end
            try
                if OPT.dryRun
                    ds = surrogateBands(Pe);
                else
                    ds = solveBands(Pe);
                    removeStrayBackslashDir(Pe.datLoc);
                end
                fit = OPT.fitnessFcn(ds, OPT.targetFreq);
            catch ME
                warning('optimize_rib_unitCell:solveFailed', ...
                    'solveBands failed at eval %d: %s', evalCounter, ME.message);
                errId = ME.identifier;
                fit = NaN;
            end
            if isempty(fit) || ~isfinite(fit)
                J = OPT.noGapPenalty;
                fit = NaN;
            else
                J = -fit;
            end
        end
        J = J + OPT.boxPenaltySlope*boxDist;

        fprintf('  eval %4d: %s  feasible=%d  fitness=%10.4g  J=%10.4g\n', evalCounter, ...
            strjoin(arrayfun(@(k) sprintf('%s=%s', OPT.vars{k}, ...
                varString(OPT.vars{k}, x(k))), 1:nv, 'UniformOutput', false), ' '), ...
            feasible, fit, J);

        % Append to the CSV log so the trajectory survives a mid-run crash.
        fid = fopen(OPT.logFile, 'a');
        fprintf(fid, '%d,%s,%s,%d,%s,%.6g,%.6g,%.6g\n', evalCounter, ...
            strjoin(arrayfun(@(k) sprintf('%.4f', varScale(OPT.vars{k})*x(k)), 1:nv, ...
                'UniformOutput', false), ','), ...
            strjoin(arrayfun(@(k) sprintf('%.6f', vc(k)), 1:nv, 'UniformOutput', false), ','), ...
            feasible, errId, boxDist, fit, J);
        fclose(fid);
    end
end

%% ===================== LOCAL HELPERS =====================
function Pe = applyDesign(P, OPT, x)
% Write the design vector x (order of OPT.vars) into a copy of P.
Pe = P;
for k = 1:numel(OPT.vars)
    Pe.(OPT.vars{k}) = x(k);
end
if ismember('L', OPT.vars)
    % Explicit half-length: drop any width so ribLength sees no P.w/P.L conflict.
    if isfield(Pe,'w'); Pe = rmfield(Pe,'w'); end
end
if OPT.fixFootprint
    % A0 = a*s + t*(2L - s)  ->  L = ((A0 - a*s)/t + s)/2
    Pe.L = ((OPT.A0 - Pe.a*Pe.s)/Pe.t + Pe.s)/2;
    if isfield(Pe,'w'); Pe = rmfield(Pe,'w'); end
end
if ismember('a', OPT.vars) && isfield(Pe,'airradUC')
    Pe.airradUC = 5*Pe.a;       % keep the unit-cell air radius on its default rule
end
end

% -------------------------------------------------------------------------

function [ok, errId, viol] = checkFeasible(Pe)
% Same checks the builder runs, plus a continuous violation measure.
%   ok     true if validateRibParams accepts the candidate
%   errId  the validation error identifier ('' when feasible)
%   viol   sum of relative constraint violations (0 when feasible), so the
%          infeasible penalty grows with the distance outside the region
Pv = ribDefaultParams(Pe);
% solveBands sets mbevenz = 1 when zSymCondition = 1 (TwoSymPlanes = 0), so
% check the candidate with the symmetry it will actually be solved with.
if isfield(Pv,'TwoSymPlanes') && ~Pv.TwoSymPlanes && isfield(Pv,'zSymCondition') && Pv.zSymCondition
    Pv.mbevenz = 1;
end
ok = true;  errId = '';
try
    ribLength(Pv, Pv.a);
    validateRibParams(Pv, 'unitcell');
catch ME
    ok = false;
    errId = ME.identifier;
end
viol = 0;
if ~ok
    L = Pv.L_ref + Pv.dLda*(Pv.a - Pv.a_ref);
    if isfield(Pv,'L') && ~isempty(Pv.L), L = Pv.L; end
    c = [Pv.t/Pv.a - 1, Pv.s/(2*max(L,eps)) - 1];
    if Pv.d == 0
        % (t + 2*r1)/a - 1 rather than r1/((a - t)/2) - 1: the latter blows up
        % as t -> a and vanishes for t > a, which would make the penalty fall
        % again past the infeasible boundary. This form grows monotonically.
        c = [c, 2*Pv.r2/Pv.t - 1, (Pv.t + 2*Pv.r1)/Pv.a - 1];
    else
        k = Pv.th*tan(Pv.d);
        c = [c, 2*k/Pv.t - 1, 2*k/Pv.s - 1];
    end
    viol = sum(max(0, c));
end
end

% -------------------------------------------------------------------------

function ds = surrogateBands(P)
% Cheap analytic stand-in for solveBands (OPT.dryRun). NOT physics: a smooth
% bump in (t, s, L) with an interior optimum, and a gap centre that scales
% like a band frequency (~ 1/a, softened by L and t). It only returns the
% fields the fitness helpers read.
L = ribLength(ribDefaultParams(P), P.a);
fc = 11.2e9*(300e-9/P.a)*sqrt(370e-9/L)*(P.t/100e-9)^0.2;
gap = 2e9*exp(-((P.t-130e-9)/60e-9)^2 - ((P.s-80e-9)/50e-9)^2 - ((L-400e-9)/120e-9)^2);
ds.sym.midGap  = [0.6*fc, fc];
ds.sym.gapSize = [0.05*gap, gap];
ds.full.midGap  = fc;
ds.full.gapSize = 0.5*gap;
ds.P = P;
end

% -------------------------------------------------------------------------

function removeStrayBackslashDir(datLoc)
% solveBands appends a literal '\' to datLoc. On Windows that is the
% separator; on macOS/Linux it creates an empty sub-folder named '\' inside
% the evaluation folder (results are still written to P.datLoc itself, and
% its cache check never hits there, which is harmless with one fresh folder
% per evaluation). Remove that empty stray folder.
if ispc, return; end
stray = [datLoc, '\'];
if exist(stray, 'dir')
    try
        rmdir(stray);   % fails, harmlessly, if it is not empty
    catch
    end
end
end

% -------------------------------------------------------------------------

function s = varColName(f)
if strcmp(f, 'd'), s = 'd_deg'; else, s = [f, '_nm']; end
end

function k = varScale(f)
if strcmp(f, 'd'), k = 180/pi; else, k = 1e9; end
end

function s = varString(f, x)
if strcmp(f, 'd')
    s = sprintf('%.2fdeg', x*180/pi);
else
    s = sprintf('%.1fnm', x*1e9);
end
end

% -------------------------------------------------------------------------

function out = mergeStruct(base, over)
% Copy fields of 'over' into 'base'; nested structs merge field by field.
out = base;
fn = fieldnames(over);
for k = 1:numel(fn)
    f = fn{k};
    if isfield(out, f) && isstruct(out.(f)) && isstruct(over.(f)) && isscalar(over.(f))
        out.(f) = mergeStruct(out.(f), over.(f));
    else
        out.(f) = over.(f);
    end
end
end

%% ===================== FITNESS HELPERS =====================
function fit = fitness_closestGap(ds, targetFreq)
% Size (Hz) of the COMPLETE bandgap whose midgap is nearest targetFreq.
% NaN if ds has no usable complete bandgap.
fit = NaN;
if ~isfield(ds,'full') || ~isfield(ds.full,'midGap') || isempty(ds.full.midGap)
    return;
end
[fit, ~] = closestGap(ds.full.midGap, ds.full.gapSize, targetFreq);
end

function fit = fitness_closestGap_sym(ds, targetFreq)
% Size (Hz) of the SYMMETRIC bandgap (ds.sym) whose midgap is nearest
% targetFreq. NaN if there is no usable symmetric bandgap.
fit = NaN;
if ~isfield(ds,'sym') || ~isfield(ds.sym,'midGap') || isempty(ds.sym.midGap)
    return;
end
[fit, ~] = closestGap(ds.sym.midGap, ds.sym.gapSize, targetFreq);
end

function fit = fitness_closestGapRatio_sym(ds, targetFreq)
% Gap-to-midgap ratio (dimensionless) of the SYMMETRIC bandgap whose midgap
% is nearest targetFreq. NaN if there is no usable symmetric bandgap.
fit = NaN;
if ~isfield(ds,'sym') || ~isfield(ds.sym,'midGap') || isempty(ds.sym.midGap)
    return;
end
[gap, mid] = closestGap(ds.sym.midGap, ds.sym.gapSize, targetFreq);
if isfinite(gap) && mid > 0
    fit = gap/mid;
end
end

function [gap, mid] = closestGap(midGap, gapSize, targetFreq)
% The gap closest to targetFreq, ignoring zero-size entries (findGaps can
% return them at band-touching points). Picking the CLOSEST gap rather than
% the widest keeps the optimizer at the operating frequency.
gap = NaN;  mid = NaN;
valid = gapSize > 0;
if ~any(valid), return; end
midGap = midGap(valid);  gapSize = gapSize(valid);
[~, k] = min(abs(midGap - targetFreq));
gap = gapSize(k);  mid = midGap(k);
end
