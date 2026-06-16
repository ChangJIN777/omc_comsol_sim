function [hxBest, hyBest, fitnessBest, dsBest] = optimize_hole_unitCell()
%OPTIMIZE_HOLE_UNITCELL Optimize hole unit-cell hx/hy for a target mechanical bandgap.
%   Nelder-Mead (fminsearch) search over [hx, hy] that maximizes a user-defined
%   fitness on the band structure returned by solveBands. Each evaluation writes
%   to a unique datLoc subfolder so solveBands never skips a solve.
%
%   Returns the best hx, hy (meters), the best fitness, and the band-structure
%   struct from a final plotted re-run.
%
%   Based on test_hole_unitCell.m.
close all;

%% ===================== TUNABLE KNOBS =====================
% --- search bounds (meters) ---
OPT.hx_min = 150e-9;   OPT.hx_max = 320e-9;
OPT.hy_min = 250e-9;   OPT.hy_max = 450e-9;

% --- initial guess (meters) ---
OPT.hx0 = 220e-9;
OPT.hy0 = 330e-9;

% --- target mechanical frequency (Hz) the gap should be centered on ---
OPT.targetFreq = 3e9;

% --- penalty returned (in the MINIMIZED objective) when no usable gap exists ---
% Large positive value (J > 0) steers the optimizer away from gap-free designs
% without aborting the run.  Must exceed any realistic fitness value (in Hz).
OPT.noGapPenalty = 1e9;

% --- root location for the per-evaluation output subfolders ---
currentDate = datestr(now,'mmddyyyy');
OPT.rootLoc = ['.\test\hole_DiamondMechanical\opt_',currentDate,'\'];

% --- fminsearch options (TolX in meters, TolFun in Hz) ---
OPT.fminOpts = optimset('Display','iter','TolX',2e-9,'TolFun',1e6, ...
    'MaxIter',60,'MaxFunEvals',120);

% --- FITNESS FUNCTION (user-configurable) ---------------------------------
% Maps the solveBands output struct ds -> scalar score (HIGHER is better).
% Default: size of the complete bandgap whose midgap is closest to targetFreq.
% Return NaN/[] if there is no usable complete bandgap (objective then applies
% noGapPenalty).
% Swap in any @(ds, targetFreq) handle returning a scalar (higher = better),
% e.g. to target a specific gap-to-midgap ratio instead of absolute gap size.
OPT.fitnessFcn = @(ds, targetFreq) fitness_closestGap(ds, targetFreq);

%% ===================== P STRUCT DEFAULTS =====================
% Identical to test_hole_unitCell.m except plot/save toggles, which the
% objective forces off per-evaluation and the final re-run turns back on.
P.xsect = 'isoFit';
P.beamMat = 'diamond';
P.celltype = 'hole';
P.unitcell = 'rectrangular';
P.maxdef = 0.15;
P.oblong = 0.7265;
P.a = 550e-9;
P.hx = OPT.hx0;
P.hy = OPT.hy0;
P.beam_width = 600e-9;
P.d_in = 0;
P.d_out = 0;
P.th = 400e-9;
P.nperiod = 1;
P.holeatedge = 0;
P.mbevenz = 0;
P.kpts = 10;
P.nbands = 5;
P.TwoSymPlanes = 0;
P.zSymCondition = 0;
P.solveasym = 1;
P.completeBandGaps = 1;
P.plotgeom = 0;          % off during optimization
P.savedat = 1;
P.savebndplot = 0;       % off during optimization
P.saveplots = 0;
P.saveMPH = 0;
P.bandStruct_2D = 0;
P.bandStructureDim = 1;
P.airDiskH = 1000e-9;
P.run_optical = 0;
P.optical_freq = 200;
P.nbeam = 2.386;
P.mbeveny = 0;
P.freq = OPT.targetFreq;
P.meshSize = 3;
P.fixed_bc = 0;
P.anisoMat = 1;
P.rxtal = 45;
P.max_dof = 3e6;

%% ===================== OPTIMIZATION =====================
evalCounter = 0;        % shared with the nested objective below

x0 = [OPT.hx0, OPT.hy0];
[xBest, fBest] = fminsearch(@objective, x0, OPT.fminOpts);

% clamp the returned optimum to the box (fminsearch is unconstrained)
hxBest = min(max(xBest(1), OPT.hx_min), OPT.hx_max);
hyBest = min(max(xBest(2), OPT.hy_min), OPT.hy_max);
fitnessBest = -fBest;   % negate: fminsearch minimizes, so fBest = -(best fitness)

fprintf('\n==================== OPTIMIZATION COMPLETE ====================\n');
fprintf('Best hx       = %.1f nm\n', hxBest*1e9);
fprintf('Best hy       = %.1f nm\n', hyBest*1e9);
fprintf('Best fitness  = %.4g Hz  = %.4f GHz (gap size)\n', fitnessBest, fitnessBest*1e-9);
fprintf('Evaluations   = %d\n', evalCounter);

%% ===================== FINAL RE-RUN WITH PLOTS =====================
Pfinal = P;
Pfinal.hx = hxBest;
Pfinal.hy = hyBest;
Pfinal.plotgeom = 1;
Pfinal.savebndplot = 1;
Pfinal.datLoc = [OPT.rootLoc, 'BEST\'];
if isfield(Pfinal,'fileBase'); Pfinal = rmfield(Pfinal,'fileBase'); end
fprintf('\nRe-running best design with plots -> %s\n', Pfinal.datLoc);
dsBest = solveBands(Pfinal);

%% ===================== NESTED OBJECTIVE =====================
    function J = objective(x)
        % x = [hx, hy] (meters). Clamp to box, solve in a unique folder,
        % score with the user fitness, return -fitness (fminsearch minimizes).
        % evalCounter lives in the enclosing function's workspace; MATLAB nested
        % functions share that workspace directly (read + write), so no globals
        % or persistent variables are needed to track evaluation count.
        evalCounter = evalCounter + 1;

        % fminsearch is unconstrained and can propose x outside [hx_min, hx_max].
        % Clamping maps the out-of-bounds proposal onto the nearest feasible point
        % rather than passing an illegal geometry to solveBands and crashing.
        hx = min(max(x(1), OPT.hx_min), OPT.hx_max);
        hy = min(max(x(2), OPT.hy_min), OPT.hy_max);

        Pe = P;
        Pe.hx = hx;
        Pe.hy = hy;
        Pe.plotgeom    = 0;
        Pe.savebndplot = 0;
        Pe.saveplots   = 0;
        % Unique folder per evaluation: solveBands loads cached results and skips
        % the COMSOL solve if <datLoc>/<fileBase>_bds.mat already exists there.
        % A fresh folder for every call guarantees a full solve every time.
        Pe.datLoc = sprintf('%seval_%04d\\', OPT.rootLoc, evalCounter);
        % Strip fileBase so CreateFileBase regenerates it from the current hx/hy;
        % a stale value copied from P would point to the wrong cache file.
        if isfield(Pe,'fileBase'); Pe = rmfield(Pe,'fileBase'); end

        try
            ds = solveBands(Pe);
        catch ME
            warning('solveBands failed at hx=%.1fnm hy=%.1fnm: %s', ...
                hx*1e9, hy*1e9, ME.message);
            J = OPT.noGapPenalty;
            return;
        end

        fit = OPT.fitnessFcn(ds, OPT.targetFreq);
        if isempty(fit) || ~isfinite(fit)
            J = OPT.noGapPenalty;
            fit = NaN;
        else
            J = -fit;   % negate: fminsearch minimizes, we want to maximize fitness
        end

        fprintf('  eval %4d: hx=%6.1fnm hy=%6.1fnm  fitness=%10.4g  J=%10.4g\n', ...
            evalCounter, hx*1e9, hy*1e9, fit, J);
    end
end

%% ===================== LOCAL FITNESS HELPER =====================
function fit = fitness_closestGap(ds, targetFreq)
% Size (Hz) of the complete bandgap whose midgap is nearest targetFreq.
% NaN if ds has no usable complete bandgap.
fit = NaN;
if ~isfield(ds,'full') || ~isfield(ds.full,'midGap') || isempty(ds.full.midGap)
    return;
end
midGap  = ds.full.midGap;
gapSize = ds.full.gapSize;
% findGaps can return zero-size entries at band-touching (degenerate) points
% that are not true bandgaps; discard them before searching for a usable gap.
valid = gapSize > 0;
if ~any(valid)
    return;
end
midGap  = midGap(valid);
gapSize = gapSize(valid);
% Select the gap closest to targetFreq (not the widest gap) so the optimizer
% stays focused on the desired operating frequency rather than drifting to a
% large but irrelevant gap elsewhere in the spectrum.
[~, k] = min(abs(midGap - targetFreq));
fit = gapSize(k);
end
