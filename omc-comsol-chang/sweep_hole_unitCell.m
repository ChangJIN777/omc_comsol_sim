function [rGrid, fitnessGrid, hxBest, hyBest, fitnessBest, dsBest] = sweep_hole_unitCell()
%SWEEP_HOLE_UNITCELL Coarse 1-D sweep of hole aspect ratio for a target mechanical bandgap.
%   Replaces the Nelder-Mead optimizer in optimize_hole_unitCell.m with a coarse
%   linspace sweep over the aspect ratio r = hx/hy. The ellipse area
%   A0 = hx0*hy0 is held FIXED, so each candidate maps r -> hx = sqrt(A0*r),
%   hy = sqrt(A0/r). Each evaluation writes to a unique datLoc subfolder so
%   solveBands never skips a solve.
%
%   Unlike the optimizer there is no objective/penalty: a failed solve simply
%   stores NaN at that grid point so the sweep continues. The best (highest-
%   fitness) point is reported, plotted, and re-run with plots enabled.
%
%   Returns the r grid, the fitness at each grid point, the best hx, hy
%   (meters), the best fitness, and the band-structure struct from the final
%   plotted re-run.
%
%   Based on optimize_hole_unitCell.m.
close all;

%% ===================== TUNABLE KNOBS =====================
% --- initial design point (meters) ---
% Sets the FIXED ellipse area A0 = hx0*hy0 that is held constant throughout the
% sweep. Only the aspect ratio r = hx/hy is varied (see below).
OPT.hx0 = 270e-9;
OPT.hy0 = 270e-9;

% --- fixed area + aspect-ratio sweep bounds (dimensionless) ---
% Area is locked to A0; the sweep walks r = hx/hy with the area-preserving map
%   hx = sqrt(A0*r),  hy = sqrt(A0/r).
%   r = 1 -> hx == hy; r < 1 -> taller/narrower hole (hy>hx);
%   r > 1 -> shorter/wider hole (hx>hy).
OPT.A0    = OPT.hx0 * OPT.hy0;   % fixed ellipse area (m^2)
OPT.r_min = 0.3;                 % smallest hx/hy (tall, narrow holes)
OPT.r_max = 1.5;                 % largest  hx/hy (short, wide holes)

% --- number of r points to evaluate (linspace from r_min to r_max) ---
OPT.nSweep = 10;

% --- target mechanical frequency (Hz) the gap should be centered on ---
OPT.targetFreq = 5e9;

% --- penalty (unused by the sweep, kept for parity with the fitness helpers) ---
% A sweep records NaN at failed/gap-free points rather than applying a penalty,
% so this is retained only so the knob block matches optimize_hole_unitCell.m.
OPT.noGapPenalty = 1e9;

% --- root location for the per-evaluation output subfolders ---
currentDate = datestr(now,'mmddyyyy');
OPT.rootLoc = ['.\test\hole_DiamondMechanical\sweep_',currentDate,'\'];

% --- FITNESS FUNCTION (user-configurable) ---------------------------------
% Maps the solveBands output struct ds -> scalar score (HIGHER is better).
% Default: size of the symmetric bandgap whose midgap is closest to targetFreq.
% Return NaN/[] if there is no usable bandgap (the sweep then stores NaN).
% Swap in any @(ds, targetFreq) handle returning a scalar (higher = better),
% e.g. to target a specific gap-to-midgap ratio instead of absolute gap size.
OPT.fitnessFcn = @(ds, targetFreq) fitness_closestGap_sym(ds, targetFreq);

%% ===================== P STRUCT DEFAULTS =====================
% Identical to test_hole_unitCell.m except plot/save toggles, which the
% sweep forces off per-evaluation and the final re-run turns back on.
P.xsect = 'isoFit';
P.beamMat = 'diamond';
P.celltype = 'hole';
P.unitcell = 'rectrangular';
P.maxdef = 0.15;
P.oblong = 0.7265;
P.a = 618e-9;
P.hx = OPT.hx0;
P.hy = OPT.hy0;
P.beam_width = 600e-9;
P.d_in = 0;
P.d_out = 0;
P.th = 380e-9;
P.nperiod = 1;
P.holeatedge = 0;
P.mbevenz = 0;
P.kpts = 10;
P.nbands = 5;
P.TwoSymPlanes = 0;
P.zSymCondition = 0;
P.solveasym = 0;    % only solve symmetric bands (skips asym FEM solve in solveBands)
P.completeBandGaps = 1;
P.plotgeom = 0;          % off during sweep
P.savedat = 1;
P.savebndplot = 0;       % off during sweep
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

%% ===================== SWEEP =====================
rGrid       = linspace(OPT.r_min, OPT.r_max, OPT.nSweep);
fitnessGrid = nan(1, OPT.nSweep);

% Open CSV log for the sweep trajectory (create fresh each run).
OPT.logFile = [OPT.rootLoc, 'sweep_log.csv'];
if ~exist(OPT.rootLoc, 'dir'); mkdir(OPT.rootLoc); end
fid = fopen(OPT.logFile, 'w');
fprintf(fid, 'eval,r,hx_nm,hy_nm,fitness_Hz\n');
fclose(fid);

for i = 1:OPT.nSweep
    r = rGrid(i);
    % Area invariant: hx*hy == A0 for every r (hx=sqrt(A0*r), hy=sqrt(A0/r)).
    hx = sqrt(OPT.A0 * r);
    hy = sqrt(OPT.A0 / r);

    Pe = P;
    Pe.hx = hx;
    Pe.hy = hy;
    Pe.plotgeom    = 0;
    Pe.savebndplot = 0;
    Pe.saveplots   = 0;
    % Unique folder per evaluation: solveBands loads cached results and skips
    % the COMSOL solve if <datLoc>/<fileBase>_bds.mat already exists there.
    % A fresh folder for every call guarantees a full solve every time.
    Pe.datLoc = sprintf('%ssweep_%04d\\', OPT.rootLoc, i);
    % Strip fileBase so CreateFileBase regenerates it from the current hx/hy;
    % a stale value copied from P would point to the wrong cache file.
    if isfield(Pe,'fileBase'); Pe = rmfield(Pe,'fileBase'); end

    try
        ds  = solveBands(Pe);
        fit = OPT.fitnessFcn(ds, OPT.targetFreq);
        if isempty(fit) || ~isfinite(fit)
            fit = NaN;
        end
    catch ME
        warning('solveBands failed at hx=%.1fnm hy=%.1fnm: %s', ...
            hx*1e9, hy*1e9, ME.message);
        fit = NaN;   % sweep, not an optimizer: record NaN and keep going
    end

    fitnessGrid(i) = fit;

    fprintf('  eval %4d/%d: r=%6.4f  hx=%6.1fnm hy=%6.1fnm  fitness=%10.4g\n', ...
        i, OPT.nSweep, r, hx*1e9, hy*1e9, fit);

    % Append row to CSV log so trajectory survives a mid-run crash.
    fid = fopen(OPT.logFile, 'a');
    fprintf(fid, '%d,%.6f,%.4f,%.4f,%.6g\n', i, r, hx*1e9, hy*1e9, fit);
    fclose(fid);
end

%% ===================== BEST POINT =====================
[fitnessBest, iBest] = max(fitnessGrid, [], 'omitnan');
if isnan(fitnessBest)
    error('sweep_hole_unitCell:noUsableGap', ...
        'No usable bandgap found at any sweep point.');
end
rBest  = rGrid(iBest);
hxBest = sqrt(OPT.A0 * rBest);
hyBest = sqrt(OPT.A0 / rBest);

fprintf('\n==================== SWEEP COMPLETE ====================\n');
fprintf('Best r        = %.4f (hx/hy)\n', rBest);
fprintf('Best hx       = %.1f nm\n', hxBest*1e9);
fprintf('Best hy       = %.1f nm\n', hyBest*1e9);
fprintf('Best fitness  = %.4g Hz  = %.4f GHz (gap size)\n', fitnessBest, fitnessBest*1e-9);
fprintf('Sweep points  = %d\n', OPT.nSweep);

%% ===================== SWEEP PLOT =====================
fig = figure('Name','aspect-ratio sweep');
plot(rGrid, fitnessGrid*1e-9, '-o', 'LineWidth', 1.2); hold on;
plot(rBest, fitnessBest*1e-9, 'ro', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', 'best');
xlabel('r = hx/hy');
ylabel('gap size (GHz)');
title(sprintf('Aspect-ratio sweep (A0 fixed, target %.2f GHz)', OPT.targetFreq*1e-9));
grid on;
saveas(fig, [OPT.rootLoc, 'sweep_fitness.png']);
savefig(fig, [OPT.rootLoc, 'sweep_fitness.fig']);

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
% Select the gap closest to targetFreq (not the widest gap) so the sweep
% stays focused on the desired operating frequency rather than drifting to a
% large but irrelevant gap elsewhere in the spectrum.
[~, k] = min(abs(midGap - targetFreq));
fit = gapSize(k);
end

%% ===================== SYMMETRIC BANDGAP FITNESS HELPER =====================
function fit = fitness_closestGap_sym(ds, targetFreq)
% Size (Hz) of the SYMMETRIC bandgap whose midgap is nearest targetFreq.
% Uses ds.sym (y-symmetric bands, mbeveny=1) rather than ds.full (complete
% bandgap across all symmetries). NaN if no usable symmetric bandgap exists.
fit = NaN;
if ~isfield(ds,'sym') || ~isfield(ds.sym,'midGap') || isempty(ds.sym.midGap)
    return;
end
midGap  = ds.sym.midGap;
gapSize = ds.sym.gapSize;
% findGaps can return zero-size entries at band-touching (degenerate) points
% that are not true bandgaps; discard them before searching for a usable gap.
valid = gapSize > 0;
if ~any(valid); return; end
midGap  = midGap(valid);
gapSize = gapSize(valid);
% Select the gap closest to targetFreq so the sweep stays focused on the
% desired operating frequency rather than drifting to an irrelevant gap.
[~, k] = min(abs(midGap - targetFreq));
fit = gapSize(k);
end
