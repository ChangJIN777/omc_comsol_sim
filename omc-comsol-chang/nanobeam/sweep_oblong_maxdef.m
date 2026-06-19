function [oVec, mVec, fitnessMap, best_oblong, best_maxdef, fitnessBest, dsBest] = ...
        sweep_oblong_maxdef()
%SWEEP_OBLONG_MAXDEF  2-D coarse sweep of nanobeam cavity taper parameters.
%
%   Sweeps P.oblong x P.maxdef on an (nOblong x nMaxdef) grid using the
%   nanobeam FEM pipeline (RunNanobeamFEM). Each grid point is evaluated in
%   a unique datLoc subfolder so RunNanobeamFEM never reuses a cached result.
%
%   P.maxdef  – maximum defect ratio (how much the lattice constant shrinks
%               at the cavity center, dimensionless, e.g. 0.15 = 15%).
%   P.oblong  – relative hx/hy taper exponent.  hy scales as
%               (1-maxdef)^(1+oblong), hx as (1-maxdef)^(1-oblong).
%               oblong=0 -> both dimensions scale equally;
%               oblong=1 -> constant hole height, only width changes.
%
%   Unlike an optimizer, failed evals (fabrication-tolerance error, COMSOL
%   crash) store NaN so the sweep can continue.  The best grid point (highest
%   fitness) is re-run with plots enabled.
%
%   Returns:
%     oVec        – (nOblong x 1) oblong sweep vector
%     mVec        – (nMaxdef x 1) maxdef sweep vector
%     fitnessMap  – (nOblong x nMaxdef) fitness at each grid point
%     best_oblong – oblong value at the best grid point
%     best_maxdef – maxdef value at the best grid point
%     fitnessBest – fitness at the best grid point
%     dsBest      – ds struct from the final re-run of the best point
%
%   Based on sweep_hole_unitCell.m; adapted for the nanobeam FEM pipeline.

close all;

%% ===================== TUNABLE KNOBS =====================

% --- oblong sweep bounds (dimensionless) ---
% Controls how hx vs hy scale differently into the defect taper.
% hy ~ (1-maxdef)^(1+oblong), hx ~ (1-maxdef)^(1-oblong).
% oblong=0 -> equal scaling; oblong=1 -> hx constant, only hy shrinks.
% For the fabricated OMC geometry (isoFit, w=800 nm) oblong ≈ 6.53 is typical.
OPT.oblong_min = 1.0;
OPT.oblong_max = 10.0;
OPT.nOblong    = 5;     % number of oblong grid points (linspace)

% --- maxdef sweep bounds (dimensionless, 0–<1) ---
% Fraction by which the mirror lattice constant is reduced at cavity centre.
% Typical range 0.02 – 0.20; larger = deeper defect = stronger localisation.
OPT.maxdef_min = 0.07;
OPT.maxdef_max = 0.2;
OPT.nMaxdef    = 5;     % number of maxdef grid points (linspace)

% --- target mechanical frequency (used for wM heatmap title only) ---
OPT.targetFreq = 10e9;  % Hz

% --- optical-wavelength fitness target and tolerance ---
OPT.targetLambda = 1550;   % nm  — target optical resonance wavelength
OPT.lambdaTol    = 50;     % nm  — Gaussian decay scale for wavelength mismatch

% --- root folder for per-evaluation subfolders ---
currentDate = datestr(now, 'mmddyyyy');
OPT.rootLoc = ['.\test\1D_OMC_hole\sweep_oblongMaxdef_', currentDate, '\'];

% --- FITNESS FUNCTION (user-configurable) ---------------------------------
% fitness = Q_opt * Q_mech * exp(-|lambda_opt - targetLambda| / lambdaTol)
% Graceful degradation: if optics are off / no high-Q optical mode found,
% the fitness falls back to Q_mech alone. Returns NaN if Q_mech is absent.
OPT.fitnessFcn = @(ds) fitness_QxQxLambda(ds, OPT.targetLambda, OPT.lambdaTol);

%% ===================== P STRUCT DEFAULTS =====================
% Geometry matches test_nanobeamRectFEM_withPML.m (fabricated OMC device).
% Adjust a, w, th, hx, hy as needed; oblong/maxdef are swept below.
P.xsect     = 'isoFit';
P.beamMat   = 'diamond';
P.celltype  = 'hole';
P.anisoMat  = 1;
P.rxtal     = 0;
P.rxtalInFilename = 1;

P.a   = 529e-9;     % nominal lattice constant (m)
P.w   = 800e-9;     % beam width (m)
P.theta = 45;       % etch angle (degrees; no effect for isoFit/rect)
P.th  = 500e-9;     % beam thickness (m)
P.hx  = 196e-9;     % nominal mirror hole height (m)
P.hy  = 651e-9;     % nominal mirror hole width  (m)

P.nholes    = 18;   % number of holes in half-beam
P.ndef      = 8;    % number of holes in half-defect region
P.maxdef    = 1 - 490/529;  % initial value — will be swept (~7.4 %)
P.oblong    = 6.53;         % initial value — will be swept
P.taperFunc = 'cubic';
P.holeatctr = 1;    % hole at cavity centre (matches fabricated device)

P.useManualDefect = 0;

P.wvgmir          = 0;
P.wgmTaper.func   = 'cubic';
P.wgmTaper.endtype = 'custom';
P.wgmTaper.a_end   = 650e-9;
P.wgmTaper.hx_end  = 343e-9;
P.wgmTaper.hy_end  = 300e-9;

P.asymCav = 0;      % disable asymmetric cavity for speed during sweep

P.lambda = 1550e-9;
P.nbeam  = 2.386;

P.stdDev    = [0, 0];
P.stdDevPos = 0;
P.asym      = 0;

% --- physics toggles (mechanical + PML; enable solveOpt/calcG/calcS if needed) ---
P.solveMech    = 1;
P.solveOpt     = 0;
P.calcG        = 1 * (P.solveMech && P.solveOpt);
P.calcS        = 0;
P.solveMechPML = 1;     % mechanical PML for radiative Q estimation

% Filter out non-localised modes; RunNanobeamFEM records warning when none pass.
P.extractLocMechModes = 1;

% --- plot / save: mode profiles saved to each per-eval subfolder ---
% Set P.plotMech = 0 / P.saveplots = 0 to skip plots and run faster.
P.plotgeom   = 0;
P.storeMPH   = 0;
P.plotMech   = 1 * P.solveMech;  % displacement/strain plot per eval
P.plotOpt    = 1 * P.solveOpt;
P.plotStrCpl = 1 * P.calcS;
P.saveplots  = 1;   % write .png/.fig into the per-eval evLoc subfolder
P.savedat    = 1;

% --- mechanical solver ---
P.mevenx  = 1;
P.meveny  = 1;
P.mevenz  = 1;
P.freq    = OPT.targetFreq;
P.mneigs  = 20;
P.mMesh   = 1;     % fine mesh (matches withPML test script)
P.mAdjMesh = 1;
P.max_dof = 5e6;

% --- optical solver (used only if solveOpt = 1) ---
P.oevenx  = (-1)^P.holeatctr;  % -1 for fundamental mode with hole at centre
P.oeveny  = -1;
P.oevenz  = 1;
P.oneigs  = 1;
P.oMesh   = 3;
P.oAdjMesh = 1;
P.airrad  = 2 * P.lambda + P.w / 2;

P.g0min = 80e3;

% SiV strain coupling axes
P.zSiV  = {[1 1 1]};
P.xSlc  = 0;
P.ySlc  = 0;
P.zSlc  = 0;
P.LStats.xmin = [];
P.LStats.xmax = [];
P.LStats.ymin = 0;
P.LStats.ymax = 60e-9;
P.LStats.zmin = P.th / 2 - 80e-9;
P.LStats.zmax = P.th / 2;

% PML settings match test_nanobeamRectFEM_withPML.m (PMLmesh not set here;
% the adaptive mesh path in SolveNanobeamFEM controls PML element size).
P.PMLmeshDiv = 20;
P.PMLLen     = 10e-6;
P.PMLstr     = 0.008;

%% ===================== BUILD SWEEP GRIDS =====================
oVec = linspace(OPT.oblong_min, OPT.oblong_max, OPT.nOblong)';
mVec = linspace(OPT.maxdef_min, OPT.maxdef_max, OPT.nMaxdef)';

fitnessMap = nan(OPT.nOblong, OPT.nMaxdef);
wM_map     = nan(OPT.nOblong, OPT.nMaxdef);
gOM_map    = nan(OPT.nOblong, OPT.nMaxdef);
LSiV_map   = nan(OPT.nOblong, OPT.nMaxdef);
Qmech_map  = nan(OPT.nOblong, OPT.nMaxdef);
Qopt_map   = nan(OPT.nOblong, OPT.nMaxdef);
lambda_map = nan(OPT.nOblong, OPT.nMaxdef);  % optical wavelength (nm)

%% ===================== OPEN CSV LOG =====================
if ~exist(OPT.rootLoc, 'dir'); mkdir(OPT.rootLoc); end
logFile = [OPT.rootLoc, 'sweep_log.csv'];
fid = fopen(logFile, 'w');
fprintf(fid, ['eval,oblong_idx,maxdef_idx,oblong,maxdef,wM_GHz,gOM_kHz,LSiV_MHz,' ...
              'Q_mech,Q_opt,lambda_opt_nm,fitness,status\n']);
fclose(fid);

evalCount = 0;
nTotal    = OPT.nOblong * OPT.nMaxdef;

%% ===================== 2-D SWEEP =====================
for io = 1:OPT.nOblong
    for im = 1:OPT.nMaxdef
        evalCount = evalCount + 1;
        oblong_i = oVec(io);
        maxdef_i = mVec(im);

        Pe = P;
        Pe.oblong = oblong_i;
        Pe.maxdef = maxdef_i;

        % Unique folder per (oblong, maxdef) — RunNanobeamFEM skips the
        % solve if the output file already exists in datLoc.
        evLoc = sprintf('%seval_%04d_ob%.4f_md%.4f%s', ...
            OPT.rootLoc, evalCount, oblong_i, maxdef_i, filesep);

        % Strip any stale fileBase so CreateFileBase regenerates it.
        if isfield(Pe, 'fileBase'); Pe = rmfield(Pe, 'fileBase'); end

        fprintf('  eval %3d/%d: oblong=%.4f  maxdef=%.4f  -> %s\n', ...
            evalCount, nTotal, oblong_i, maxdef_i, evLoc);

        wM_i     = NaN;
        gOM_i    = NaN;
        LSiV_i   = NaN;
        Qmech_i  = NaN;
        Qopt_i   = NaN;
        lambda_i = NaN;
        fit      = NaN;
        status_i = 'failed';

        try
            [ds, ~] = RunNanobeamFEM(Pe, evLoc);

            % Treat a missing or empty ds as a silent failure.
            if isempty(ds) || ~isstruct(ds)
                error('sweep_oblong_maxdef:emptyDs', 'RunNanobeamFEM returned empty ds');
            end

            % --- check for localised mechanical modes ---
            hasMech = false;
            if isfield(ds, 'mfem') && isfield(ds.mfem, 'freqs') && ...
                    ~isempty(ds.mfem.freqs)
                hasMech = true;
            end
            if ~hasMech
                warning('sweep_oblong_maxdef:noMechModes', ...
                    'eval %d (oblong=%.4f, maxdef=%.4f): no localised mechanical modes found', ...
                    evalCount, oblong_i, maxdef_i);
                status_i = 'no_mech_modes';
            else
                status_i = 'ok';
            end

            % --- extract mechanical frequency (prefer coupled-mode result) ---
            if isfield(ds, 'cpl') && isfield(ds.cpl, 'mechFreq') && ...
                    ~isempty(ds.cpl.mechFreq)
                wM_i = max(real(ds.cpl.mechFreq));
            elseif hasMech
                wM_i = max(real(ds.mfem.freqs));
            end

            % --- check for high-Q optical modes (when solveOpt = 1) ---
            if Pe.solveOpt
                hasOpt = isfield(ds, 'ofem') && isfield(ds.ofem, 'freqs') && ...
                         ~isempty(ds.ofem.freqs);
                if ~hasOpt
                    warning('sweep_oblong_maxdef:noOptModes', ...
                        'eval %d (oblong=%.4f, maxdef=%.4f): no high-Q optical modes found', ...
                        evalCount, oblong_i, maxdef_i);
                    if strcmp(status_i, 'ok')
                        status_i = 'no_opt_modes';
                    else
                        status_i = [status_i, '+no_opt_modes'];
                    end
                end
            end

            % --- extract gOM ---
            if isfield(ds, 'cpl') && isfield(ds.cpl, 'gMax') && ...
                    ~isempty(ds.cpl.gMax)
                gOM_i = max(abs(ds.cpl.gMax));
            end

            % --- extract LSiV ---
            if isfield(ds, 'cpl') && isfield(ds.cpl, 'SiV') && ...
                    isfield(ds.cpl.SiV, 'LSiVMaxMode') && ...
                    ~isempty(ds.cpl.SiV.LSiVMaxMode)
                LSiV_i = max(real(ds.cpl.SiV.LSiVMaxMode));
            end

            % --- extract mechanical radiative Q (requires solveMechPML) ---
            [Qmech_i, ~] = local_getQmech(ds);

            % --- extract optical Q and wavelength (requires solveOpt) ---
            [Qopt_i, lambda_i] = local_getQopt(ds, OPT.targetLambda);

            % Fitness is NaN when modes are absent; sweep continues.
            if hasMech
                fit = OPT.fitnessFcn(ds);
                if isempty(fit) || ~isfinite(fit)
                    fit = NaN;
                end
            end

        catch ME
            warning('sweep_oblong_maxdef:evalFailed', ...
                'eval %d (oblong=%.4f, maxdef=%.4f) failed: %s', ...
                evalCount, oblong_i, maxdef_i, ME.message);
            status_i = 'failed';
        end

        fitnessMap(io, im) = fit;
        wM_map(io, im)     = wM_i;
        gOM_map(io, im)    = gOM_i;
        LSiV_map(io, im)   = LSiV_i;
        Qmech_map(io, im)  = Qmech_i;
        Qopt_map(io, im)   = Qopt_i;
        lambda_map(io, im) = lambda_i;

        fprintf(['         wM=%.3f GHz   gOM=%.2f kHz   LSiV=%.4f MHz   ' ...
                 'Qm=%.3g   Qo=%.3g   lam=%.1f nm   fitness=%.4g   [%s]\n'], ...
            wM_i*1e-9, gOM_i*1e-3, LSiV_i*1e-6, Qmech_i, Qopt_i, lambda_i, fit, status_i);

        % Append to CSV immediately so the log survives a mid-run crash.
        fid = fopen(logFile, 'a');
        fprintf(fid, '%d,%d,%d,%.6f,%.6f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%s\n', ...
            evalCount, io, im, oblong_i, maxdef_i, ...
            wM_i*1e-9, gOM_i*1e-3, LSiV_i*1e-6, Qmech_i, Qopt_i, lambda_i, fit, status_i);
        fclose(fid);
    end
end

%% ===================== FIND BEST POINT =====================
validMask = isfinite(fitnessMap);
if ~any(validMask(:))
    error('sweep_oblong_maxdef:noUsableResult', ...
        'No usable result found at any sweep point.');
end

[fitnessBest, linIdx] = max(fitnessMap(:));
[ioBest, imBest]      = ind2sub([OPT.nOblong, OPT.nMaxdef], linIdx);
best_oblong = oVec(ioBest);
best_maxdef = mVec(imBest);

fprintf('\n==================== SWEEP COMPLETE ====================\n');
fprintf('Best oblong   = %.4f\n', best_oblong);
fprintf('Best maxdef   = %.4f\n', best_maxdef);
fprintf('Best fitness  = %.4g\n', fitnessBest);
fprintf('Grid size     = %d x %d  (%d total evals)\n', OPT.nOblong, OPT.nMaxdef, nTotal);
fprintf('Log written   -> %s\n', logFile);

%% ===================== SWEEP PLOTS =====================

% --- 2D fitness heatmap ---
fig1 = figure('Name', 'Fitness heatmap (oblong x maxdef)', 'Color', 'w');
imagesc(mVec, oVec, fitnessMap);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('maxdef');
ylabel('oblong');
title(sprintf('Fitness map (target %.2f GHz)', OPT.targetFreq * 1e-9));
hold on;
plot(best_maxdef, best_oblong, 'r*', 'MarkerSize', 14, 'LineWidth', 1.5);
saveas(fig1, [OPT.rootLoc, 'sweep_fitness_map.png']);
savefig(fig1, [OPT.rootLoc, 'sweep_fitness_map.fig']);

% --- 2D wM heatmap ---
fig2 = figure('Name', 'wM heatmap (oblong x maxdef)', 'Color', 'w');
imagesc(mVec, oVec, wM_map * 1e-9);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('maxdef');
ylabel('oblong');
title('Mechanical frequency \omega_M (GHz)');
hold on;
plot(best_maxdef, best_oblong, 'r*', 'MarkerSize', 14, 'LineWidth', 1.5);
saveas(fig2, [OPT.rootLoc, 'sweep_wM_map.png']);
savefig(fig2, [OPT.rootLoc, 'sweep_wM_map.fig']);

% --- gOM heatmap (only if data available) ---
if any(isfinite(gOM_map(:)))
    fig3 = figure('Name', 'gOM heatmap (oblong x maxdef)', 'Color', 'w');
    imagesc(mVec, oVec, gOM_map * 1e-3);
    set(gca, 'YDir', 'normal');
    colorbar;
    xlabel('maxdef');
    ylabel('oblong');
    title('|g_{OM}| (kHz)');
    hold on;
    plot(best_maxdef, best_oblong, 'r*', 'MarkerSize', 14, 'LineWidth', 1.5);
    saveas(fig3, [OPT.rootLoc, 'sweep_gOM_map.png']);
    savefig(fig3, [OPT.rootLoc, 'sweep_gOM_map.fig']);
end

% --- LSiV heatmap (only if data available) ---
if any(isfinite(LSiV_map(:)))
    fig4 = figure('Name', 'LSiV heatmap (oblong x maxdef)', 'Color', 'w');
    imagesc(mVec, oVec, LSiV_map * 1e-6);
    set(gca, 'YDir', 'normal');
    colorbar;
    xlabel('maxdef');
    ylabel('oblong');
    title('\lambda_{SiV} (MHz)');
    hold on;
    plot(best_maxdef, best_oblong, 'r*', 'MarkerSize', 14, 'LineWidth', 1.5);
    saveas(fig4, [OPT.rootLoc, 'sweep_LSiV_map.png']);
    savefig(fig4, [OPT.rootLoc, 'sweep_LSiV_map.fig']);
end

% --- Qmech heatmap (only if PML data available) ---
if any(isfinite(Qmech_map(:)))
    fig5 = figure('Name', 'Qmech heatmap (oblong x maxdef)', 'Color', 'w');
    imagesc(mVec, oVec, Qmech_map);
    set(gca, 'YDir', 'normal');
    colorbar;
    xlabel('maxdef');
    ylabel('oblong');
    title('Mechanical radiative Q_{mech}');
    hold on;
    plot(best_maxdef, best_oblong, 'r*', 'MarkerSize', 14, 'LineWidth', 1.5);
    saveas(fig5, [OPT.rootLoc, 'sweep_Qmech_map.png']);
    savefig(fig5, [OPT.rootLoc, 'sweep_Qmech_map.fig']);
end

% --- Qopt heatmap (only if optical data available) ---
if any(isfinite(Qopt_map(:)))
    fig6 = figure('Name', 'Qopt heatmap (oblong x maxdef)', 'Color', 'w');
    imagesc(mVec, oVec, Qopt_map);
    set(gca, 'YDir', 'normal');
    colorbar;
    xlabel('maxdef');
    ylabel('oblong');
    title('Optical Q_{opt}');
    hold on;
    plot(best_maxdef, best_oblong, 'r*', 'MarkerSize', 14, 'LineWidth', 1.5);
    saveas(fig6, [OPT.rootLoc, 'sweep_Qopt_map.png']);
    savefig(fig6, [OPT.rootLoc, 'sweep_Qopt_map.fig']);
end

% --- optical wavelength heatmap (only if optical data available) ---
if any(isfinite(lambda_map(:)))
    fig7 = figure('Name', 'lambda_opt heatmap (oblong x maxdef)', 'Color', 'w');
    imagesc(mVec, oVec, lambda_map);
    set(gca, 'YDir', 'normal');
    colorbar;
    xlabel('maxdef');
    ylabel('oblong');
    title('Optical wavelength \lambda_{opt} (nm)');
    hold on;
    plot(best_maxdef, best_oblong, 'r*', 'MarkerSize', 14, 'LineWidth', 1.5);
    saveas(fig7, [OPT.rootLoc, 'sweep_lambda_map.png']);
    savefig(fig7, [OPT.rootLoc, 'sweep_lambda_map.fig']);
end

%% ===================== FINAL RE-RUN WITH PLOTS =====================
Pfinal = P;
Pfinal.oblong    = best_oblong;
Pfinal.maxdef    = best_maxdef;
Pfinal.plotgeom  = 1;
Pfinal.plotMech  = 1 * P.solveMech;
Pfinal.plotOpt   = 1 * P.solveOpt;
Pfinal.plotStrCpl = 1 * P.calcS;
Pfinal.saveplots = 1;
if isfield(Pfinal, 'fileBase'); Pfinal = rmfield(Pfinal, 'fileBase'); end

bestLoc = [OPT.rootLoc, 'BEST', filesep];
fprintf('\nRe-running best design with plots -> %s\n', bestLoc);
[dsBest, ~] = RunNanobeamFEM(Pfinal, bestLoc);

end  % function sweep_oblong_maxdef

%% ===================== LOCAL FITNESS HELPERS =====================

function fit = fitness_QxQxLambda(ds, targetLambdaNm, tolNm)
% fitness = Q_opt * Q_mech * exp(-|lambda_opt - targetLambda| / tol)
% Graceful degradation:
%   no Q_mech (PML off / absent)    -> NaN (unusable point)
%   optics off / no optical mode    -> Q_mech alone
fit = NaN;
[Qmech, ~]        = local_getQmech(ds);
[Qopt,  lambdaNm] = local_getQopt(ds, targetLambdaNm);
if isnan(Qmech); return; end
if isnan(Qopt) || isnan(lambdaNm)
    fit = Qmech;
else
    fit = Qopt * Qmech * exp(-abs(lambdaNm - targetLambdaNm) / tolNm);
end
end

function [Qmech, fMech] = local_getQmech(ds)
% Mechanical radiative Q of the best localised cavity mode (solveMechPML=1).
% Selects the localised mode with the highest localization ratio; falls back
% to the highest-Q solved mode if localisation filter was not applied.
Qmech = NaN;  fMech = NaN;
if ~isfield(ds, 'mfem'); return; end
mfem = ds.mfem;
if ~isfield(mfem, 'QAll') || isempty(mfem.QAll); return; end
if isfield(mfem, 'locInd') && ~isempty(mfem.locInd) && ...
        isfield(mfem, 'locRatio') && ~isempty(mfem.locRatio)
    [~, k] = max(mfem.locRatio(mfem.locInd));
    idx    = mfem.locInd(k);
else
    [~, idx] = max(mfem.QAll);
end
if idx <= numel(mfem.QAll)
    Qmech = mfem.QAll(idx);
end
if isfield(mfem, 'freqs') && idx <= numel(mfem.freqs)
    fMech = mfem.freqs(idx);
end
end

function [Qopt, lambdaNm] = local_getQopt(ds, targetLambdaNm)
% Optical Q and wavelength (nm) of the high-Q mode nearest targetLambdaNm.
% Uses the high-Q subset ds.ofem.Q / ds.ofem.lambda (Q>1e4); falls back to
% the full ds.ofem.QAll / ds.ofem.lambdaAll set if the subset is empty.
% ds.ofem.lambda is in meters; returned lambdaNm is in nm.
Qopt = NaN;  lambdaNm = NaN;
if ~isfield(ds, 'ofem'); return; end
ofem = ds.ofem;
if isfield(ofem, 'lambda') && ~isempty(ofem.lambda) && ...
        isfield(ofem, 'Q')   && ~isempty(ofem.Q)
    lam = ofem.lambda(:) * 1e9;  Q = ofem.Q(:);
elseif isfield(ofem, 'lambdaAll') && ~isempty(ofem.lambdaAll) && ...
        isfield(ofem, 'QAll')      && ~isempty(ofem.QAll)
    lam = ofem.lambdaAll(:) * 1e9;  Q = ofem.QAll(:);
else
    return;
end
[~, k]   = min(abs(lam - targetLambdaNm));
Qopt     = Q(k);
lambdaNm = lam(k);
end

function fit = fitness_maxGOM(ds)
% Returns max |g_OM| in Hz. NaN if gOM not computed (calcG = 0).
fit = NaN;
if isfield(ds, 'cpl') && isfield(ds.cpl, 'gMax') && ~isempty(ds.cpl.gMax)
    fit = max(abs(ds.cpl.gMax));
end
end

function fit = fitness_maxLSiV(ds)
% Returns max lambda_SiV in Hz. NaN if strain coupling not computed.
fit = NaN;
if isfield(ds, 'cpl') && isfield(ds.cpl, 'SiV') && ...
        isfield(ds.cpl.SiV, 'LSiVMaxMode') && ~isempty(ds.cpl.SiV.LSiVMaxMode)
    fit = max(real(ds.cpl.SiV.LSiVMaxMode));
end
end
