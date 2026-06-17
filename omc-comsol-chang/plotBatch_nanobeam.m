% plotBatch_nanobeam.m
% -------------------------------------------------------------------------
% Batch loader + plotter for nanobeam FEM results produced by RunNanobeamFEM.
%
% Recursively finds every *.mat file under a user-specified folder, loads the
% single saved variable `ds` from each, extracts key mechanical /
% optomechanical / spin-phonon results into a summary table, prints the table
% to the console, and generates a few overview figures saved as .png.
%
% Each .mat saves one struct `ds`. Fields used here (all guarded with isfield,
% since not every run computes coupling):
%   ds.P              full parameter struct (SI units, metres)
%   ds.mfem.freqs     mechanical eigenfrequencies [Hz]
%   ds.mfem.xzpf      zero-point fluctuation [m]
%   ds.cpl.mechFreq   best mechanical frequency [Hz]      (if P.calcG)
%   ds.cpl.gMax       best |g_OM| [Hz]                    (if P.calcG)
%   ds.cpl.mSol       index/indices of best mode          (if P.calcG)
%   ds.cpl.SiV.LSiVMaxMode  spin-phonon coupling per mode [Hz]  (if P.calcS)
%
% No external toolboxes required (uses base MATLAB only).
%
% Usage: set `datLoc` below to the root results folder, then run this script.
% -------------------------------------------------------------------------

clear; close all;

%% ----------------------------- user settings ----------------------------
% Root folder to search (recursively) for *.mat result files.
datLoc = 'D:\Files\OMC-SiV\RectOMC\';

% If true, a corrupt / incompatible .mat just prints a warning and is skipped
% instead of aborting the whole batch.
skipErrors = true;

% Append a trailing filesep so dir([datLoc '**/*.mat']) is well-formed.
if ~isempty(datLoc) && ~ismember(datLoc(end), {'\','/'})
    datLoc = [datLoc, filesep];
end
if ~exist(datLoc, 'dir')
    error('plotBatch_nanobeam:badPath', 'datLoc does not exist: %s', datLoc);
end

%% --------------------------- find result files --------------------------
% '**' makes dir recurse into subfolders (base MATLAB, R2016b+).
files = dir([datLoc, '**', filesep, '*.mat']);

% Drop directory entries that happen to match (safety).
files = files(~[files.isdir]);

if isempty(files)
    fprintf('No .mat files found under %s\n', datLoc);
    return;
end
fprintf('Found %d .mat file(s) under %s\n', numel(files), datLoc);

%% ------------------------- per-file extraction --------------------------
nF = numel(files);

% Preallocate result columns (NaN = "not available for this file").
filename  = strings(nF, 1);
a_nm      = nan(nF, 1);
hx_nm     = nan(nF, 1);
hy_nm     = nan(nF, 1);
w_nm      = nan(nF, 1);
th_nm     = nan(nF, 1);
wM_GHz    = nan(nF, 1);
gOM_kHz   = nan(nF, 1);
LSiV_MHz  = nan(nF, 1);
xzpf_fm   = nan(nF, 1);

keep = false(nF, 1);   % mark files that loaded successfully

for fi = 1:nF
    fPath = fullfile(files(fi).folder, files(fi).name);

    try
        S = load(fPath, 'ds');
        if ~isfield(S, 'ds')
            error('plotBatch_nanobeam:noDs', ...
                  'no variable ''ds'' in file');
        end
        ds = S.ds;

        % Some error_*.mat runs save an empty ds — treat as unusable.
        if isempty(ds) || ~isstruct(ds)
            error('plotBatch_nanobeam:emptyDs', 'ds is empty / not a struct');
        end

        filename(fi) = string(files(fi).name);

        % --- geometry (m -> nm), guarded ---
        if isfield(ds, 'P')
            P = ds.P;
            a_nm(fi)  = getNm(P, 'a');
            hx_nm(fi) = getNm(P, 'hx');
            hy_nm(fi) = getNm(P, 'hy');
            w_nm(fi)  = getNm(P, 'w');
            th_nm(fi) = getNm(P, 'th');
        end

        % --- mechanical frequency [GHz]: prefer best-coupled mode ---
        if isfield(ds, 'cpl') && isfield(ds.cpl, 'mechFreq') && ...
                ~isempty(ds.cpl.mechFreq)
            wM_GHz(fi) = max(real(ds.cpl.mechFreq)) * 1e-9;
        elseif isfield(ds, 'mfem') && isfield(ds.mfem, 'freqs') && ...
                ~isempty(ds.mfem.freqs)
            wM_GHz(fi) = max(real(ds.mfem.freqs)) * 1e-9;
        end

        % --- best |g_OM| [kHz] ---
        if isfield(ds, 'cpl') && isfield(ds.cpl, 'gMax') && ...
                ~isempty(ds.cpl.gMax)
            gOM_kHz(fi) = max(abs(ds.cpl.gMax)) * 1e-3;
        end

        % --- best spin-phonon coupling [MHz] ---
        if isfield(ds, 'cpl') && isfield(ds.cpl, 'SiV') && ...
                isfield(ds.cpl.SiV, 'LSiVMaxMode') && ...
                ~isempty(ds.cpl.SiV.LSiVMaxMode)
            LSiV_MHz(fi) = max(real(ds.cpl.SiV.LSiVMaxMode)) * 1e-6;
        end

        % --- xzpf of the best mode [fm] ---
        % ds.cpl.mSol may be a scalar or a vector of mode indices; guard both.
        if isfield(ds, 'cpl') && isfield(ds.cpl, 'mSol') && ...
                ~isempty(ds.cpl.mSol) && ...
                isfield(ds, 'mfem') && isfield(ds.mfem, 'xzpf') && ...
                ~isempty(ds.mfem.xzpf)
            mSol = ds.cpl.mSol;
            valid = mSol >= 1 & mSol <= numel(ds.mfem.xzpf);
            if any(valid)
                xz = ds.mfem.xzpf(mSol(valid));
                xzpf_fm(fi) = max(abs(xz)) * 1e15;
            end
        end

        keep(fi) = true;

    catch ME
        if skipErrors
            warning('plotBatch_nanobeam:skip', ...
                    'Skipping "%s": %s', files(fi).name, ME.message);
        else
            rethrow(ME);
        end
    end
end

% Keep only files that loaded.
idx = find(keep);
if isempty(idx)
    fprintf('No usable result files loaded.\n');
    return;
end

% Assemble summary table (struct2table is base MATLAB).
T = table(filename(idx), a_nm(idx), hx_nm(idx), hy_nm(idx), w_nm(idx), ...
          th_nm(idx), wM_GHz(idx), gOM_kHz(idx), LSiV_MHz(idx), xzpf_fm(idx), ...
    'VariableNames', {'filename','a_nm','hx_nm','hy_nm','w_nm','th_nm', ...
                      'wM_GHz','gOM_kHz','LSiV_MHz','xzpf_fm'});

% Sort by mechanical frequency for a readable table / index axis.
[~, sortOrder] = sort(T.wM_GHz);
T = T(sortOrder, :);

%% ----------------------------- print table ------------------------------
fprintf('\n');
fprintf('%-4s %8s %7s %7s %7s %7s %9s %9s %10s %9s  %s\n', ...
        'idx','a/nm','hx/nm','hy/nm','w/nm','th/nm', ...
        'wM/GHz','gOM/kHz','LSiV/MHz','xzpf/fm','file');
fprintf('%s\n', repmat('-', 1, 110));
for k = 1:height(T)
    fprintf('%-4d %8.1f %7.1f %7.1f %7.1f %7.1f %9.3f %9.1f %10.3f %9.2f  %s\n', ...
            k, T.a_nm(k), T.hx_nm(k), T.hy_nm(k), T.w_nm(k), T.th_nm(k), ...
            T.wM_GHz(k), T.gOM_kHz(k), T.LSiV_MHz(k), T.xzpf_fm(k), ...
            shortName(char(T.filename(k))));
end
fprintf('\n');

% Flags for what coupling data is actually present in the batch.
haveGOM  = any(~isnan(T.gOM_kHz));
haveLSiV = any(~isnan(T.LSiV_MHz));

%% ----------------------- Figure 1: wM vs index --------------------------
f1 = figure('Name', 'Mechanical frequencies', 'Color', 'w');
ax1 = axes(f1); hold(ax1, 'on'); box(ax1, 'on'); grid(ax1, 'on');
xIdx = (1:height(T))';
scatter(ax1, xIdx, T.wM_GHz, 60, 'filled');
xlabel(ax1, 'file index');
ylabel(ax1, '\omega_M / GHz');
title(ax1, 'Mechanical frequencies');
% Label x-ticks with the shortened wM_* part of each filename.
shortLabels = arrayfun(@(k) wmTag(char(T.filename(k))), xIdx, ...
                       'UniformOutput', false);
set(ax1, 'XTick', xIdx, 'XTickLabel', shortLabels, ...
         'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
xlim(ax1, [0.5, height(T) + 0.5]);
saveFig(f1, [datLoc, 'batch_fig1_wM_vs_index.png']);

%% --------------------- Figure 2: gOM vs wM (colour LSiV) ----------------
if haveGOM
    f2 = figure('Name', 'g_OM vs wM', 'Color', 'w');
    ax2 = axes(f2); hold(ax2, 'on'); box(ax2, 'on'); grid(ax2, 'on');
    if haveLSiV
        scatter(ax2, T.wM_GHz, T.gOM_kHz, 60, T.LSiV_MHz, 'filled');
        cb = colorbar(ax2);
        cb.Label.String = '\lambda_{SiV} / MHz';
    else
        scatter(ax2, T.wM_GHz, T.gOM_kHz, 60, 'filled');
    end
    xlabel(ax2, '\omega_M / GHz');
    ylabel(ax2, '|g_{OM}| / kHz');
    title(ax2, 'g_{OM} vs \omega_M');
    saveFig(f2, [datLoc, 'batch_fig2_gOM_vs_wM.png']);
else
    fprintf('Figure 2 skipped: no g_OM data in batch.\n');
end

%% --------------------- Figure 3: LSiV vs wM -----------------------------
if haveLSiV
    f3 = figure('Name', 'Spin-phonon coupling vs wM', 'Color', 'w');
    ax3 = axes(f3); hold(ax3, 'on'); box(ax3, 'on'); grid(ax3, 'on');
    scatter(ax3, T.wM_GHz, T.LSiV_MHz, 60, 'filled');
    xlabel(ax3, '\omega_M / GHz');
    ylabel(ax3, '\lambda_{SiV} / MHz');
    title(ax3, 'Spin-phonon coupling vs \omega_M');
    saveFig(f3, [datLoc, 'batch_fig3_LSiV_vs_wM.png']);
else
    fprintf('Figure 3 skipped: no spin-phonon (LSiV) data in batch.\n');
end

fprintf('Done. Figures saved to %s\n', datLoc);

%% ============================ local functions ===========================
function v = getNm(P, fld)
    % Return P.(fld)*1e9 (m -> nm) if the field exists & is non-empty, else NaN.
    if isfield(P, fld) && ~isempty(P.(fld))
        v = P.(fld)(1) * 1e9;
    else
        v = NaN;
    end
end

function s = shortName(name)
    % Trim the long parameter prefix for console readability: keep from the
    % 'wM_' tag onward if present, otherwise return the full name.
    tok = regexp(name, 'wM_.*', 'match', 'once');
    if isempty(tok)
        s = name;
    else
        s = tok;
    end
end

function s = wmTag(name)
    % Extract just the 'wM_<freq>GHz' chunk for axis tick labels.
    tok = regexp(name, 'wM_[^_]*GHz', 'match', 'once');
    if isempty(tok)
        s = shortName(name);
    else
        s = tok;
    end
end

function saveFig(figH, outPath)
    % Save figure as .png, warning (not erroring) on failure.
    try
        saveas(figH, outPath);
        fprintf('  saved %s\n', outPath);
    catch ME
        warning('plotBatch_nanobeam:save', ...
                'Could not save %s: %s', outPath, ME.message);
    end
end
