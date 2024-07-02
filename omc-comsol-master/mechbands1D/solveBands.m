% This function computes the mechanical band structure for a nanobeam for
% both symmetric and antisymmetric bands. Plots and displacement/strain
% profiles are also generated. The input structure P is assumed to have the
% following fields, plus additional fields defined in drawUnitCell.m:
%
% P.E = Young's modulus (in Pa);
% P.rho = density (in kg/m^3);
% P.nu = Poisson's ratio;
% P.nbands = # of bands to solve for;
% P.kpts = # of k-points;
% P.bfreq = frequency at which to start looking for bands;
% P.maxdof = max # of degrees of freedom in mesh;
% P.evenz = 1 if the mode is even in z;
% P.eveny = 1 if the mode is even in y;
% P.solveasym = 1 if we want to solve asymmetric bands
% P.savedat = 1 if we want to save the data structure
% P.savebndplot = 1 if we want to save band diagrams
%
% Michael Burek, 2/22/16
%
% Modified by Cleaven Chia, 09/06/16:
% - removed TeX interpreter code
%
% TODO (20160929 1227H):
% - plot full unit cell geometry using sector 3D dataset

function ds = solveBands(P,datLoc)
%append \ to end of datLoc if not present
if ~strcmp(datLoc(end),'\')
    datLoc = [datLoc,'\'];
end
% create directory to save files
if ~exist(datLoc,'dir')
    mkdir(datLoc)
end

% utilFilesPath = 'G:\My Drive\OMC\COMSOL\FEMscripts20181224';
% addpath(utilFilesPath)

% create base folder name
fBase = ['a',num2str(P.a*1e9,'%.0f'),'nm_',...
    'w',num2str(P.w*1e9,'%.0f'),'nm_', ...
    'hx',num2str(P.hx*1e9,'%.0f'),'nm_', ...
    'hy',num2str(P.hy*1e9,'%.0f'),'nm_',...
    'rxt',num2str(P.rxtal),'o'];
if strcmp(P.xsect,'tri')
    fBase = [num2str(P.theta,'%.0f'),'o_',fBase];
elseif strcmp(P.xsect,'rect')
    fBase = ['th',num2str(P.th*1e9,'%.0f'),'nm_',fBase];
end
if isfield(P,'prefname')
    fBase = [P.prefname,'_',fBase];
end
if isfield(P,'filRad') && P.filRad > 0
    fBase = [fBase,'_Rfil_',num2str(P.filRad*1e9,'%.1f'),'nm'];
end

disp(fBase)
% if ~isfield(P,'fileBase')
%     P = CreateFileBase(P);
%     disp(P.fileBase)
% %     ds.P = P;
% end
% fBase = P.fileBase;
P.fileBase = fBase;
if isempty(dir([datLoc,fBase,'*_bds.mat']))
    tStart = tic;
    %% solve bands
    display('Solving with symmetric boundary condition');
    P.mbeveny = 1;
    if strcmp(P.xsect,'rect')
        display('in Y and Z (rect x-sect)');
        P.mbevenz = 1;
    end
    sym = RunNanobeamBands(P,datLoc);
    [sym.midGap,sym.gapSize] = findGaps(sym);   %find gaps
    
    if P.solveasym == 1
        display('Solving with anti-symmetric boundary condition');
        P.mbeveny = -1;
        if strcmp(P.xsect,'rect')
            display('in Y and Z (rect x-sect)');
            P.mbevenz = -1;
        end
        asym = RunNanobeamBands(P,datLoc);
        [asym.midGap,asym.gapSize] = findGaps(asym);    %find gaps
        
        if strcmp(P.xsect,'rect')
            display('Solving with Y anti-symmetric and Z symmetric boundary condition');
            P.mbeveny = -1;
            P.mbevenz = 1;
            asymYsymZ = RunNanobeamBands(P,datLoc);
            [asymYsymZ.midGap,asymYsymZ.gapSize] = findGaps(asymYsymZ);    %find gaps
            
            display('Solving with Y symmetric and Z antisymmetric boundary condition');
            P.mbeveny = 1;
            P.mbevenz = -1;
            symYasymZ = RunNanobeamBands(P,datLoc);
            [symYasymZ.midGap,symYasymZ.gapSize] = findGaps(symYasymZ);    %find gaps
        end
    end
    
%     close all;
    
    
    %%
    % write to data structure
    ds.sym = sym;
    if P.solveasym == 1
        ds.asym = asym;
        if strcmp(P.xsect,'rect')
            ds.asymYsymZ = asymYsymZ;
            ds.symYasymZ = symYasymZ;
        end
    end
    tEnd = toc(tStart);
    disp(['Simulation time = ',num2str(tEnd/60,'%.2f'),'mins'])
    
else
    display('Data folder exists in working directory, loading...')
%     ds.sym = [];
%     ds.asym = [];
    dsFiles = dir([datLoc,fBase,'*_bds.mat']);
    load([datLoc,dsFiles.name]);
    sym = ds.sym;
    if P.solveasym == 1
        asym = ds.asym;
        if strcmp(P.xsect,'rect')
            asymYsymZ = ds.asymYsymZ;
            symYasymZ = ds.symYasymZ;
        end
    end
end

%% find complete bandgaps
if P.completeBandGaps
    if P.solveasym == 1
        complete.F = [sym.F,asym.F];
        if strcmp(P.xsect,'rect')
            complete.F = [complete.F,asymYsymZ.F,symYasymZ.F];
        end
    end
    [complete.midGap,complete.gapSize] = findGaps(complete);
    ds.complete = complete;
end
    

%% plot bandstructure
if P.savebndplot == 1
    figure; 
    set(gcf,'position',[50 50 640 800],'units','pixels')
    hold on

    maxFreqs = [0 0 0 0];
    pHdl = [];
    % plot bands
    % call plot twice so legend appears properly
    p1 = plot(sym.kx_norm,sym.F(:,1)*1e-9,'-k','linewidth',3,'DisplayName','sym');
    plot(sym.kx_norm,sym.F*1e-9,'-k','linewidth',3);
    maxFreqs(1) = max(sym.F(:)*1e-9);
    pHdl(end+1) = p1;
    if P.solveasym == 1
        p2 = plot(asym.kx_norm,asym.F(:,1)*1e-9,'b:','linewidth',2,'DisplayName','asym');
        plot(asym.kx_norm,asym.F*1e-9,'b:','linewidth',2);
        maxFreqs(2) = max(asym.F(:)*1e-9);
        pHdl(end+1) = p2;

        if strcmp(P.xsect,'rect')
            p3 = plot(asymYsymZ.kx_norm,asymYsymZ.F(:,1)*1e-9,'-.','linewidth',0.5,'DisplayName','asym Y, sym Z','Color',80/256*[1 1 1]);
            plot(asymYsymZ.kx_norm,asymYsymZ.F*1e-9,'-.','linewidth',0.5,'Color',80/256*[1 1 1]);
            p4 = plot(symYasymZ.kx_norm,symYasymZ.F(:,1)*1e-9,'m--','linewidth',0.5,'DisplayName','sym Y, asym Z');
            plot(symYasymZ.kx_norm,symYasymZ.F*1e-9,'m--','linewidth',0.5);
            maxFreqs(3) = max(asymYsymZ.F(:)*1e-9);
            maxFreqs(4) = max(symYasymZ.F(:)*1e-9);
            pHdl(end+1) = p3;
            pHdl(end+1) = p4;
        end
    end

    % plot bandgaps
    
    for k = 1:length(sym.gapSize)
        bgp = patch([0 1 1 0],1e-9.*(sym.midGap(k) + 0.5*[sym.gapSize(k) ...
            sym.gapSize(k) -sym.gapSize(k) -sym.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
        alpha(bgp,0.5);
    end
    
    if P.completeBandGaps && ~isempty(complete.midGap)
        for k = 1:length(complete.gapSize)
            bgp = patch([0 1 1 0],1e-9.*(complete.midGap(k) + 0.5*[complete.gapSize(k) ...
                complete.gapSize(k) -complete.gapSize(k) -complete.gapSize(k)]),180/255*[1 0 0],'EdgeColor','none');
            alpha(bgp,0.2);
        end
    else

        if P.solveasym == 1 && strcmp(P.xsect,'tri')
            for k = 1:length(asym.gapSize)
                bgp = patch([0 1 1 0],1e-9.*(asym.midGap(k) + 0.5*[asym.gapSize(k) ...
                    asym.gapSize(k) -asym.gapSize(k) -asym.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
                alpha(bgp,0.4);
            end
        end
        
        % plot symmetric midgap frequencies
        for k = 1:length(sym.midGap)
            midfreqs = sym.midGap(k)*ones(length(sym.kx_norm),1);
            plot(sym.kx_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
        end
    end

    

    legend(pHdl,'location','southoutside','orientation','horizontal');
    xlabel('k_x / (\pi/a)','FontSize',12);
    ylabel('Frequency (GHz)','FontSize',12);
    maxY = max(maxFreqs);
    axis([0 1 0 maxY]);
    set(gca,'XTick',[0; 1]);
%         set(gca,'XTickLabel',{'G','C'},'fontname','symbol','fontsize',16)
    set(gca,'XTickLabel',{'\Gamma','X'},'fontsize',12)

    if strcmp(P.xsect,'tri')
        thStr = ['\theta = ',num2str(P.theta),'^o, '];
    elseif strcmp(P.xsect,'rect')
        thStr = ['th = ',num2str(P.th*1e9,'%.0f'),'nm, '];
    end
    if isfield(P,'filRad') && P.filRad > 0
        filRadTtl = [', Rf = ',num2str(P.filRad*1e9,'%.1f'),'nm'];
    else
        filRadTtl = '';
    end
    bandtitle = {[thStr,...
                  'a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                  'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
                  'rtxal = ',num2str(P.rxtal,'%.0f'),'^o, ',...
                  filRadTtl];
                 ['hx = ',num2str(P.hx*1e9,'%.0f'),'nm, ', ...
                  'hy = ',num2str(P.hy*1e9,'%.0f'),'nm, ',...
                  'hx/a = ',num2str(P.hx/P.a,'%.4f'),', ' ...
                  'hy/w = ',num2str(P.hy/P.w,'%.4f')]};
    title(bandtitle);
    box on
    hold off

    % save band diagram as .png and .fig
    set(gcf,'PaperPositionMode','auto')
    fname_fig = [fBase,'_fullBands.png'];
    pathFig = [datLoc,fname_fig];
    saveas(gcf,pathFig);
    fname_fig = [fBase,'_fullBands.fig'];
    pathFig = [datLoc,fname_fig];
    saveas(gcf,pathFig);

    close
    % plot w mech freqs
%         plot(sym.kx_norm,5.28*ones(size(sym.kx_norm)),'--g','linewidth',2);
%         plot(sym.kx_norm,8.03*ones(size(sym.kx_norm)),'--g','linewidth',2);
%         fname_fig = [fBase,'_fullBandsWMechFreq.png'];
%         pathFig = [datLoc,fname_fig];
%         saveas(gcf,pathFig);
%         fname_fig = [fBase,'_fullBandsWMechFreq.fig'];
%         pathFig = [datLoc,fname_fig];
%         saveas(gcf,pathFig);
end

% save data structures
if P.savedat
    fname_mat = [fBase,'_bds.mat'];
    pathMat = [datLoc,fname_mat];
    save(pathMat,'ds');
end

% rmpath(utilFilesPath)

