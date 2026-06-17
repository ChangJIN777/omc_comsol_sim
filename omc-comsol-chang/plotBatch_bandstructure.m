%% plotBatch_bandstructure.m
% Batch re-plot of band structure data.
%
% Finds every *_bds.mat file under datLoc (recursively), loads the single
% `ds` struct each one contains, and reproduces the band structure plot that
% solveBands.m generates (lines 254-444 of solveBands.m, both 1D and 2D
% cases, TwoSymPlanes 0 or 1). Each figure is saved next to its source .mat
% file as <original_name>_bandplot.png.
%
% The `ds` struct layout depends on ds.P.TwoSymPlanes / ds.P.solveasym:
%   TwoSymPlanes=1: ds.symy_symz, ds.symy_asymz, ds.asymy_symz,
%                   ds.asymy_asymz, ds.full
%   TwoSymPlanes=0: ds.sym, [ds.asym], ds.full
% Each band struct has .F [kpts x nbands], .k_norm [kpts x 1], .midGap,
% .gapSize. ds.P carries TwoSymPlanes, solveasym, bandStruct_2D, celltype
% and the geometry fields used to build the title.

clear; close all;

%% --- user settings ---------------------------------------------------------
datLoc     = 'C:\Users\alber\Documents\Github\omc_comsol_sim\omc-comsol-chang\test\';   % root folder to scan
skipErrors = true;    % true -> log and continue past files that fail
% ---------------------------------------------------------------------------

% make sure the root path ends with a filesep so the dir pattern is well formed
if ~isempty(datLoc) && ~endsWith(datLoc,filesep)
    datLoc = [datLoc,filesep];
end

% recursive search for all band structure files
files = dir([datLoc,'**',filesep,'*_bds.mat']);

if isempty(files)
    fprintf('No *_bds.mat files found under %s\n',datLoc);
    return
end

fprintf('Found %d *_bds.mat file(s) under %s\n',numel(files),datLoc);

for ii = 1:numel(files)
    matPath = fullfile(files(ii).folder,files(ii).name);
    [~,nameNoExt,~] = fileparts(files(ii).name);
    pngPath = fullfile(files(ii).folder,[nameNoExt,'_bandplot.png']);

    fprintf('[%d/%d] %s\n',ii,numel(files),matPath);

    if skipErrors
        try
            S  = load(matPath,'ds');
            fh = plotBandStruct(S.ds);
            saveas(fh,pngPath);
            close(fh);
            fprintf('        saved %s\n',pngPath);
        catch ME
            fprintf(2,'        SKIPPED (%s): %s\n',ME.identifier,ME.message);
        end
    else
        S  = load(matPath,'ds');
        fh = plotBandStruct(S.ds);
        saveas(fh,pngPath);
        close(fh);
        fprintf('        saved %s\n',pngPath);
    end
end

fprintf('Done.\n');


%% ==========================================================================
function fh = plotBandStruct(ds)
% Reproduce the band structure plot from solveBands.m for a loaded `ds`.
% Returns the figure handle.

P = ds.P;

fh = figure; hold on

% ------------------------------------------------------------------ bands ---
if P.TwoSymPlanes
    symy_symz   = ds.symy_symz;
    symy_asymz  = ds.symy_asymz;
    asymy_symz  = ds.asymy_symz;
    asymy_asymz = ds.asymy_asymz;

    p1 = plot(symy_symz.k_norm,  symy_symz.F*1e-9,  '-k', 'linewidth',2,'DisplayName','symy_symz');
    p2 = plot(symy_asymz.k_norm, symy_asymz.F*1e-9, '--b','linewidth',2,'DisplayName','symy_asymz');
    p3 = plot(asymy_symz.k_norm, asymy_symz.F*1e-9, '--r','linewidth',2,'DisplayName','asymy_symz');
    p4 = plot(asymy_asymz.k_norm,asymy_asymz.F*1e-9,'--m','linewidth',2,'DisplayName','asymy_asymz');
    % one line handle per group is enough for the legend
    legHandles = [p1(1) p2(1) p3(1) p4(1)];
else
    sym = ds.sym;
    p1  = plot(sym.k_norm,sym.F*1e-9,'-k','linewidth',2,'DisplayName','sym');
    legHandles = p1(1);
    if P.solveasym && isfield(ds,'asym') && ~isempty(ds.asym)
        asym = ds.asym;
        p2 = plot(asym.k_norm,asym.F*1e-9,'--b','linewidth',2,'DisplayName','asym');
        legHandles = [legHandles p2(1)];
    end
end

full = ds.full;

% ------------------------------------------------------------ gaps + axes ---
if P.bandStruct_2D
    % ============================ 2D case ===============================
    % symmetric bandgaps
    if isfield(sym,'gapSize')
        for k = 1:length(sym.gapSize)
            bgp = patch([0 3 3 0],1e-9.*(sym.midGap(k) + 0.5*[sym.gapSize(k) ...
                sym.gapSize(k) -sym.gapSize(k) -sym.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
            alpha(bgp,0.5);
        end
    end

    % asymmetric bandgaps
    if exist('asym','var') && isfield(asym,'gapSize')
        for k = 1:length(asym.gapSize)
            bgp = patch([0 3 3 0],1e-9.*(asym.midGap(k) + 0.5*[asym.gapSize(k) ...
                asym.gapSize(k) -asym.gapSize(k) -asym.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
            alpha(bgp,0.2);
        end
    end

    % full bandgaps
    if isfield(full,'midGap') && ~isempty(full.midGap)
        for k = 1:length(full.gapSize)
            bgp = patch([0 3 3 0],1e-9.*(full.midGap(k) + 0.5*[full.gapSize(k) ...
                full.gapSize(k) -full.gapSize(k) -full.gapSize(k)]),180/255*[1 0 0],'EdgeColor','none');
            alpha(bgp,0.2);
        end
    end

    % symmetric midgap frequencies
    if isfield(sym,'midGap')
        for k = 1:length(sym.midGap)
            midfreqs = sym.midGap(k)*ones(length(sym.k_norm),1);
            plot(sym.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
        end
    end

    % complete midgap frequencies
    if isfield(full,'midGap')
        for k = 1:length(full.midGap)
            midfreqs = full.midGap(k)*ones(length(sym.k_norm),1);
            plot(sym.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
        end
    end

    xlabel('k','FontSize',12);
    ylabel('Frequency (GHz)','FontSize',12);
    if exist('asym','var') && isfield(asym,'F')
        amax = max([sym.F(:);asym.F(:)])*1e-9;
    else
        amax = max(sym.F(:))*1e-9;
    end
    axis([0 3 0 amax]);
    set(gca,'XTick',[0; 1; 2; 3]);
    set(gca,'XTickLabel',{'\Gamma','X','M','\Gamma'},'fontsize',12)

else
    % ============================ 1D case ===============================
    if P.TwoSymPlanes
        % --- two symmetry planes (four structs) ---
        % symmetric bandgaps
        for k = 1:length(symy_symz.gapSize)
            bgp = patch([0 1 1 0],1e-9.*(symy_symz.midGap(k) + 0.5*[symy_symz.gapSize(k) ...
                symy_symz.gapSize(k) -symy_symz.gapSize(k) -symy_symz.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
            alpha(bgp,0.5);
        end

        % asymmetric bandgaps
        for k = 1:length(asymy_asymz.gapSize)
            bgp = patch([0 1 1 0],1e-9.*(asymy_asymz.midGap(k) + 0.5*[asymy_asymz.gapSize(k) ...
                asymy_asymz.gapSize(k) -asymy_asymz.gapSize(k) -asymy_asymz.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
            alpha(bgp,0.2);
        end

        % full bandgaps
        if isfield(full,'midGap') && ~isempty(full.midGap)
            for k = 1:length(full.gapSize)
                bgp = patch([0 1 1 0],1e-9.*(full.midGap(k) + 0.5*[full.gapSize(k) ...
                    full.gapSize(k) -full.gapSize(k) -full.gapSize(k)]),180/255*[1 0 0],'EdgeColor','none');
                alpha(bgp,0.2);
            end
        end

        % symmetric midgap frequencies
        for k = 1:length(symy_symz.midGap)
            midfreqs = symy_symz.midGap(k)*ones(length(symy_symz.k_norm),1);
            plot(symy_symz.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
        end

        % asymmetric midgap frequencies
        for k = 1:length(asymy_asymz.midGap)
            midfreqs = asymy_asymz.midGap(k)*ones(length(asymy_asymz.k_norm),1);
            plot(asymy_asymz.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
        end

        % complete midgap frequencies
        if isfield(full,'midGap')
            for k = 1:length(full.midGap)
                midfreqs = full.midGap(k)*ones(length(symy_symz.k_norm),1);
                plot(symy_symz.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
            end
        end

        amax = max([symy_symz.F(:);asymy_asymz.F(:);symy_asymz.F(:);asymy_symz.F(:)])*1e-9;
    else
        % --- single symmetry plane (sym / asym) ---
        % symmetric bandgaps
        for k = 1:length(sym.gapSize)
            bgp = patch([0 1 1 0],1e-9.*(sym.midGap(k) + 0.5*[sym.gapSize(k) ...
                sym.gapSize(k) -sym.gapSize(k) -sym.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
            alpha(bgp,0.5);
        end

        % asymmetric bandgaps
        if P.solveasym && exist('asym','var')
            for k = 1:length(asym.gapSize)
                bgp = patch([0 1 1 0],1e-9.*(asym.midGap(k) + 0.5*[asym.gapSize(k) ...
                    asym.gapSize(k) -asym.gapSize(k) -asym.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
                alpha(bgp,0.2);
            end
        end

        % full bandgaps
        if isfield(full,'midGap') && ~isempty(full.midGap)
            for k = 1:length(full.gapSize)
                bgp = patch([0 1 1 0],1e-9.*(full.midGap(k) + 0.5*[full.gapSize(k) ...
                    full.gapSize(k) -full.gapSize(k) -full.gapSize(k)]),180/255*[1 0 0],'EdgeColor','none');
                alpha(bgp,0.2);
            end
        end

        % symmetric midgap frequencies
        for k = 1:length(sym.midGap)
            midfreqs = sym.midGap(k)*ones(length(sym.k_norm),1);
            plot(sym.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
        end

        % asymmetric midgap frequencies
        if P.solveasym && exist('asym','var')
            for k = 1:length(asym.midGap)
                midfreqs = asym.midGap(k)*ones(length(asym.k_norm),1);
                plot(asym.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
            end
        end

        % complete midgap frequencies
        if isfield(full,'midGap')
            for k = 1:length(full.midGap)
                midfreqs = full.midGap(k)*ones(length(sym.k_norm),1);
                plot(sym.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
            end
        end

        if P.solveasym && exist('asym','var')
            amax = max([sym.F(:);asym.F(:)])*1e-9;
        else
            amax = max(sym.F(:))*1e-9;
        end
    end

    xlabel('k','FontSize',12);
    ylabel('Frequency (GHz)','FontSize',12);
    axis([0 1 0 amax]);
    set(gca,'XTick',[0; 1]);
    set(gca,'XTickLabel',{'\Gamma','X'},'fontsize',12)

    % --------------------------------------------------------- title -------
    bandtitle = buildBandTitle(P);
    if ~isempty(bandtitle)
        title(bandtitle);
    end
end

% legend over the band lines only (not patches / midgap lines)
if ~isempty(legHandles)
    legend(legHandles,'Location','best','Interpreter','none');
end

box on
hold off
end


%% ==========================================================================
function bandtitle = buildBandTitle(P)
% Geometry string for the plot title, keyed off P.celltype.
% Mirrors the title block in solveBands.m (lines 423-501).

bandtitle = {};

switch P.celltype
    case '2D_ribs'
        bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
            'w = ',num2str(P.w*1e9,'%.0f'),'nm',...
            'hi = ',num2str(P.hi*1e9,'%.0f'),'nm, ',...
            'wi = ',num2str(P.wi*1e9,'%.0f'),'nm',...
            'ho = ',num2str(P.ho*1e9,'%.0f'),'nm, ',...
            'wo = ',num2str(P.wo*1e9,'%.0f'),'nm',...
            'ai = ',num2str(P.ai*1e9,'%.0f'),'nm'];
            ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
    case 'hole'
        bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
            'hx = ',num2str(P.hx*1e9,'%.0f'),'nm, ',...
            'hy = ',num2str(P.hy*1e9,'%.0f'),'nm',...
            'beam_width = ',num2str(P.beam_width*1e9,'%.0f'),'nm'];
            ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
    case 'cross'
        bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
            'h = ',num2str(P.h*1e9,'%.0f'),'nm, ',...
            'w = ',num2str(P.wc*1e9,'%.0f'),'nm',...
            'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
            'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm'];
            ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
    case 'boomerang'
        bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
            'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
            'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
            'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
            'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm'];
            ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
    case 'boomerang_lower'
        bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
            'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
            'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
            'd = ',num2str(P.d*1e9,'%.0f'),'nm',...
            'h = ',num2str(P.h*1e9,'%.0f'),'nm',...
            'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
            'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm'];
            ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
    case 'snowflake'
        bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
            'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
            'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
            'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
            'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm'];
            ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
    case 'Snowflake_strip'
        bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
            'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
            'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
            'wo = ',num2str(P.wo*1e9,'%.0f'),'nm',...
            'wi = ',num2str(P.wi*1e9,'%.0f'),'nm',...
            'ho = ',num2str(P.ho*1e9,'%.0f'),'nm',...
            'hi = ',num2str(P.hi*1e9,'%.0f'),'nm',...
            'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
            'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm'];
            ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
    case {'boomerang_strip','boomerang_strip_v2'}
        bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
            'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
            'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
            'wo = ',num2str(P.wo*1e9,'%.0f'),'nm',...
            'wi = ',num2str(P.wi*1e9,'%.0f'),'nm',...
            'ho = ',num2str(P.ho*1e9,'%.0f'),'nm',...
            'hi = ',num2str(P.hi*1e9,'%.0f'),'nm',...
            'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
            'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm'];
            ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
end
end
