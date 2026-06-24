function [best_oblong, best_maxdef, fitnessBest, dsBest] = optimize_oblong_maxdef()
%OPTIMIZE_OBLONG_MAXDEF  Nelder-Mead optimization of nanobeam cavity taper.
%
%   Optimizes the nanobeam cavity over two parameters using fminsearch
%   (Nelder-Mead simplex):
%       x(1) = defectAspectRatio  (mapped to P.oblong, see below)
%       x(2) = maxdef             (maximum defect ratio, dimensionless)
%
%   This is the optimizer counterpart of sweep_oblong_maxdef.m: instead of
%   evaluating a fixed (ndefectAspectRatio x nMaxdef) grid, it lets
%   fminsearch search the same 2-D parameter space, maximizing the fitness
%       fitness = Q_opt * Q_mech * exp(-|lambda_opt - targetLambda| / lambdaTol)
%   via the same fitness_QxQxLambda logic.  Because fminsearch minimizes, the
%   objective returns the negated fitness.
%
%   Parameter mapping (identical to the sweep):
%       mirrorAspectRatio = P.hy / P.hx;
%       r = defectAspectRatio / mirrorAspectRatio;
%       Pe.oblong = log(r) / log(maxdef);
%
%   Each objective evaluation runs RunNanobeamFEM in a unique per-eval
%   subfolder under OPT.rootLoc so no cached result is ever reused.  Box
%   constraints are enforced with a penalty: if either parameter falls
%   outside its [min, max] bound, the objective returns 1e12 without calling
%   RunNanobeamFEM.  Any error during an evaluation is caught and returns
%   1e12 so the optimizer never crashes.
%
%   Returns:
%     best_oblong  – oblong value at the optimum (mapped from defectAspectRatio)
%     best_maxdef  – maxdef value at the optimum
%     fitnessBest  – fitness (positive, un-negated) at the optimum
%     dsBest       – ds struct from the final re-run of the best point
%
%   Based on sweep_oblong_maxdef.m; adapted to use fminsearch.

close all;

%% ===================== TUNABLE KNOBS =====================

% --- starting point (same as sweep centre) ---
% hy ~ (1-maxdef)^(1+oblong), hx ~ (1-maxdef)^(1-oblong).
% oblong=0 -> equal scaling; oblong=1 -> hx constant, only hy shrinks.
oblong_0 = 1.96;
maxdef_0 = 0.22;
defectAspectRatio_0 = (651e-9/196e-9) * (1-maxdef_0)^(2*oblong_0);  % hy/hx * (1-maxdef)^(2*oblong)

% --- defectAspectRatio box constraints (dimensionless) ---
OPT.defectAspectRatio_min = defectAspectRatio_0*0.5;
OPT.defectAspectRatio_max = defectAspectRatio_0*1.5;

% --- maxdef box constraints (dimensionless, 0–<1) ---
% Fraction by which the mirror lattice constant is reduced at cavity centre.
OPT.maxdef_min = maxdef_0*0.5;
OPT.maxdef_max = maxdef_0*1.5;

% --- target mechanical frequency (used for solver target) ---
OPT.targetFreq = 7e9;  % Hz

% --- optical-wavelength fitness target and tolerance ---
OPT.targetLambda = 1550;   % nm  — target optical resonance wavelength
OPT.lambdaTol    = 100;    % nm  — Gaussian decay scale for wavelength mismatch

% --- root folder for per-evaluation subfolders ---
currentDate = datestr(now, 'mmddyyyy');
OPT.rootLoc = ['.\test\1D_OMC_hole\optimize_oblongMaxdef_trial1_', currentDate, '\'];

% --- FITNESS FUNCTION (user-configurable) ---------------------------------
% fitness = g_OM^2 * Q_opt * Q_mech * exp(-|lambda_opt - targetLambda| / lambdaTol)
% Cooperativity-like figure of merit. Graceful degradation: if g_OM is absent
% (P.calcG=0) the g^2 factor is dropped; if optics are off, falls back to
% g_OM^2 * Q_mech. Returns NaN if Q_mech is absent.
OPT.fitnessFcn = @(ds) fitness_coopProduct(ds, OPT.targetLambda, OPT.lambdaTol);

%% ===================== P STRUCT DEFAULTS =====================
% Geometry matches test_nanobeamRectFEM_withPML.m (fabricated OMC device).
P.xsect     = 'isoFit';
P.beamMat   = 'diamond';
P.celltype  = 'hole';
P.anisoMat  = 1;
P.rxtal     = 0;
P.rxtalInFilename = 1;

P.a   = 529e-9;     % nominal lattice constant (m)
P.w   = 750e-9;     % beam width (m)
P.theta = 45;       % etch angle (degrees; no effect for isoFit/rect)
P.th  = 500e-9;     % beam thickness (m)
P.hx  = 196e-9;     % nominal mirror hole height (m)
P.hy  = 578e-9;     % nominal mirror hole width  (m)

P.nholes    = 18;   % number of holes in half-beam
P.ndef      = 8;    % number of holes in half-defect region
P.maxdef    = maxdef_0;  % initial value — will be optimized (~7.4 %)
P.oblong    = oblong_0;         % initial value — will be optimized
P.taperFunc = 'cubic';
P.holeatctr = 1;    % hole at cavity centre (matches fabricated device)

P.useManualDefect = 0;

P.wvgmir          = 0;
P.wgmTaper.func   = 'cubic';
P.wgmTaper.endtype = 'custom';
P.wgmTaper.a_end   = 650e-9;
P.wgmTaper.hx_end  = 343e-9;
P.wgmTaper.hy_end  = 300e-9;

P.asymCav = 0;      % disable asymmetric cavity for speed during optimization

P.lambda = 1550e-9;
P.nbeam  = 2.386;

P.stdDev    = [0, 0];
P.stdDevPos = 0;
P.asym      = 0;

% --- physics toggles (mechanical + optics + PML) ---
P.solveMech    = 1;
P.solveOpt     = 1;
P.calcG        = 1 * (P.solveMech && P.solveOpt);
P.calcS        = 0;
P.solveMechPML = 1;     % mechanical PML for radiative Q estimation

% Filter out non-localised modes; RunNanobeamFEM records warning when none pass.
P.extractLocMechModes = 1;

% --- plot / save: mode profiles saved to each per-eval subfolder ---
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

% PML settings match test_nanobeamRectFEM_withPML.m
P.PMLmeshDiv = 20;
P.PMLLen     = 10e-6;
P.PMLstr     = 0.008;

%% ===================== OPEN CSV LOG =====================
if ~exist(OPT.rootLoc, 'dir'); mkdir(OPT.rootLoc); end
logFile = [OPT.rootLoc, 'optimize_log.csv'];
fid = fopen(logFile, 'w');
fprintf(fid, ['eval,oblong,maxdef,wM_GHz,gOM_kHz,LSiV_MHz,' ...
              'Q_mech,Q_opt,lambda_opt_nm,fitness,status\n']);
fclose(fid);

% Eval counter (shared with the objective via nested scope).
evalCount = 0;

% Track the best (positive) fitness seen so far for the final re-run.
bestSoFar.fitness = -Inf;
bestSoFar.oblong  = NaN;
bestSoFar.maxdef  = NaN;
bestSoFar.defectAspectRatio = NaN;

%% ===================== FMINSEARCH SETUP =====================
x0 = [defectAspectRatio_0, maxdef_0];

options = optimset('Display', 'iter', ...
                   'TolX',        1e-4, ...
                   'TolFun',      1e-6, ...
                   'MaxIter',     100, ...
                   'MaxFunEvals', 200);

fprintf('\n==================== OPTIMIZATION START ====================\n');
fprintf('Start defectAspectRatio = %.6f  (oblong_0 = %.4f)\n', defectAspectRatio_0, oblong_0);
fprintf('Start maxdef            = %.6f\n', maxdef_0);
fprintf('defectAspectRatio bounds = [%.6f, %.6f]\n', ...
    OPT.defectAspectRatio_min, OPT.defectAspectRatio_max);
fprintf('maxdef bounds            = [%.6f, %.6f]\n', ...
    OPT.maxdef_min, OPT.maxdef_max);
fprintf('Log written -> %s\n\n', logFile);

% Wrap the objective in a try/catch so the optimizer never crashes.
objFcn = @(x) safeObjective(x);

[xBest, negFitBest] = fminsearch(objFcn, x0, options);

%% ===================== RESOLVE BEST POINT =====================
best_defectAspectRatio = xBest(1);
best_maxdef            = xBest(2);

% Recompute the oblong mapping at the optimum.
mirrorAspectRatio = P.hy / P.hx;
r_best            = best_defectAspectRatio / mirrorAspectRatio;
best_oblong       = log(r_best) / (2*log(1 - best_maxdef));

fitnessBest = -negFitBest;   % un-negate to report the maximized fitness

% If fminsearch returned a penalized/failed point but a valid point was seen
% earlier, prefer the best valid point that was actually evaluated.
if (~isfinite(fitnessBest) || fitnessBest <= 0) && isfinite(bestSoFar.fitness)
    fprintf(['fminsearch endpoint was penalized/invalid; ' ...
             'falling back to best valid eval seen.\n']);
    best_defectAspectRatio = bestSoFar.defectAspectRatio;
    best_maxdef            = bestSoFar.maxdef;
    best_oblong            = bestSoFar.oblong;
    fitnessBest            = bestSoFar.fitness;
end

fprintf('\n==================== OPTIMIZATION COMPLETE ====================\n');
fprintf('Best defectAspectRatio = %.6f\n', best_defectAspectRatio);
fprintf('Best oblong            = %.4f\n', best_oblong);
fprintf('Best maxdef            = %.4f\n', best_maxdef);
fprintf('Best fitness           = %.4g\n', fitnessBest);
fprintf('Total evaluations      = %d\n', evalCount);
fprintf('Log written            -> %s\n', logFile);

%% ===================== FINAL RE-RUN WITH PLOTS =====================
Pfinal = P;
Pfinal.oblong     = best_oblong;
Pfinal.maxdef     = best_maxdef;
Pfinal.plotgeom   = 1;
Pfinal.plotMech   = 1 * P.solveMech;
Pfinal.plotOpt    = 1 * P.solveOpt;
Pfinal.plotStrCpl = 1 * P.calcS;
Pfinal.saveplots  = 1;
if isfield(Pfinal, 'fileBase'); Pfinal = rmfield(Pfinal, 'fileBase'); end

bestLoc = [OPT.rootLoc, 'BEST', filesep];
fprintf('\nRe-running best design with plots -> %s\n', bestLoc);
[dsBest, ~] = RunNanobeamFEM(Pfinal, bestLoc);

%% ===================== NESTED OBJECTIVE =====================
    function J = safeObjective(x)
    % Negated-fitness objective for fminsearch with box-constraint penalty.
    % Returns 1e12 on out-of-bounds parameters or any evaluation error so
    % the optimizer never crashes.
        try
            J = objective(x);
        catch ME
            warning('optimize_oblong_maxdef:objError', ...
                'objective failed at x=[%.4f, %.4f]: %s', x(1), x(2), ME.message);
            J = 1e12;
        end
    end

    function J = objective(x)
    % Core objective: maps x -> (oblong, maxdef), runs RunNanobeamFEM in a
    % unique folder, evaluates fitness, logs to CSV, and returns -fitness.
        defectAspectRatio_i = x(1);
        maxdef_i            = x(2);

        % --- box-constraint penalty (no FEM call when out of bounds) ---
        if defectAspectRatio_i < OPT.defectAspectRatio_min || ...
           defectAspectRatio_i > OPT.defectAspectRatio_max || ...
           maxdef_i < OPT.maxdef_min || maxdef_i > OPT.maxdef_max
            J = 1e12;
            return;
        end

        evalCount = evalCount + 1;

        Pe = P;

        % --- map defect aspect ratio to oblong (same as the sweep) ---
        mirrorAspectRatio = P.hy / P.hx;
        r = defectAspectRatio_i / mirrorAspectRatio;

        Pe.maxdef = maxdef_i;
        Pe.oblong = log(r) / (2*log(1 - maxdef_i));

        % Unique folder per evaluation — RunNanobeamFEM never reuses a cache.
        evLoc = sprintf('%seval_%04d_ob%.4f_md%.4f%s', ...
            OPT.rootLoc, evalCount, Pe.oblong, maxdef_i, filesep);

        % Strip any stale fileBase so CreateFileBase regenerates it.
        if isfield(Pe, 'fileBase'); Pe = rmfield(Pe, 'fileBase'); end

        fprintf('  eval %3d: oblong=%.4f  maxdef=%.4f  -> %s\n', ...
            evalCount, Pe.oblong, maxdef_i, evLoc);

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
                error('optimize_oblong_maxdef:emptyDs', ...
                    'RunNanobeamFEM returned empty ds');
            end

            % --- check for localised mechanical modes ---
            hasMech = false;
            if isfield(ds, 'mfem') && isfield(ds.mfem, 'freqs') && ...
                    ~isempty(ds.mfem.freqs)
                hasMech = true;
            end
            if ~hasMech
                warning('optimize_oblong_maxdef:noMechModes', ...
                    'eval %d (oblong=%.4f, maxdef=%.4f): no localised mechanical modes found', ...
                    evalCount, Pe.oblong, maxdef_i);
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
                    warning('optimize_oblong_maxdef:noOptModes', ...
                        'eval %d (oblong=%.4f, maxdef=%.4f): no high-Q optical modes found', ...
                        evalCount, Pe.oblong, maxdef_i);
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

            % Fitness is NaN when modes are absent.
            if hasMech
                fit = OPT.fitnessFcn(ds);
                if isempty(fit) || ~isfinite(fit)
                    fit = NaN;
                end
            end

        catch ME
            warning('optimize_oblong_maxdef:evalFailed', ...
                'eval %d (oblong=%.4f, maxdef=%.4f) failed: %s', ...
                evalCount, Pe.oblong, maxdef_i, ME.message);
            status_i = 'failed';
        end

        fprintf(['         wM=%.3f GHz   gOM=%.2f kHz   LSiV=%.4f MHz   ' ...
                 'Qm=%.3g   Qo=%.3g   lam=%.1f nm   fitness=%.4g   [%s]\n'], ...
            wM_i*1e-9, gOM_i*1e-3, LSiV_i*1e-6, Qmech_i, Qopt_i, lambda_i, fit, status_i);

        % Append to CSV immediately so the log survives a mid-run crash.
        fidLog = fopen(logFile, 'a');
        fprintf(fidLog, '%d,%.6f,%.6f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%s\n', ...
            evalCount, Pe.oblong, maxdef_i, ...
            wM_i*1e-9, gOM_i*1e-3, LSiV_i*1e-6, Qmech_i, Qopt_i, lambda_i, fit, status_i);
        fclose(fidLog);

        % --- form the negated objective ---
        % NaN fitness (no usable mode) is penalized so fminsearch avoids it.
        if isfinite(fit)
            J = -fit;
            if fit > bestSoFar.fitness
                bestSoFar.fitness = fit;
                bestSoFar.oblong  = Pe.oblong;
                bestSoFar.maxdef  = maxdef_i;
                bestSoFar.defectAspectRatio = defectAspectRatio_i;
            end
        else
            J = 1e12;
        end
    end

end  % function optimize_oblong_maxdef

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

function fit = fitness_coopProduct(ds, targetLambdaNm, tolNm)
% fitness = g_OM^2 * Q_opt * Q_mech * exp(-|lambda_opt - targetLambda| / tol)
% Cooperativity-like figure of merit. Graceful degradation:
%   no Q_mech (PML off / absent)    -> NaN (unusable point)
%   no g_OM   (P.calcG = 0)         -> g^2 factor dropped; reverts to Q_opt*Q_mech*gate
%   optics off / no optical mode    -> g_OM^2 * Q_mech
fit = NaN;
[Qmech, ~]        = local_getQmech(ds);
[Qopt,  lambdaNm] = local_getQopt(ds, targetLambdaNm);
if isnan(Qmech); return; end
gOM     = local_getGOM(ds);
gFactor = 1;
if ~isnan(gOM); gFactor = gOM^2; end
if isnan(Qopt) || isnan(lambdaNm)
    fit = gFactor * Qmech;
else
    fit = gFactor * Qopt * exp(-abs(lambdaNm - targetLambdaNm) / tolNm);
end
end

function fit = fitness_maxGOM(ds)
% Returns max |g_OM| in Hz. NaN if gOM not computed (calcG = 0).
fit = NaN;
if isfield(ds, 'cpl') && isfield(ds.cpl, 'gMax') && ~isempty(ds.cpl.gMax)
    fit = max(abs(ds.cpl.gMax));
end
end

function gOM = local_getGOM(ds)
% Max |g_OM| in Hz of the cavity mode, or NaN if not computed (P.calcG = 0).
gOM = NaN;
if isfield(ds, 'cpl') && isfield(ds.cpl, 'gMax') && ~isempty(ds.cpl.gMax)
    gOM = max(abs(ds.cpl.gMax));
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
