% Optical bandstructure simulations in COMSOL

function ds = solveBands(P,datLoc)
%append \ to end of datLoc if not present
if ~strcmp(datLoc(end),'\')
    datLoc = [datLoc,'\'];
end
% create directory to save files
if ~exist(datLoc,'dir')
    mkdir(datLoc)
end

% create base folder name
% create base folder name
if ~isfield(P,'fileBase')
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

    P.fileBase = fBase;
end
fBase = P.fileBase;

if isempty(dir([datLoc,fBase,'*_bds.mat']))
    %% solve bands
    display('Solving for TE modes');
    P.obeveny = -1;
    if strcmp(P.xsect,'rect') && P.obevenz ~= 0
        P.obevenz = 1;
    end
    TE = RunNanobeamBands(P,datLoc);
%     [sym.midGap,sym.gapSize] = findGaps(sym);   %find gaps
    
    if P.solveTM == 1
        display('Solving for TM modes');
        P.obeveny = 1;
        if strcmp(P.xsect,'rect') && P.obevenz ~= 0
            P.obevenz = -1;
        end
        TM = RunNanobeamBands(P,datLoc);
%         [asym.midGap,asym.gapSize] = findGaps(asym);    %find gaps
        
%         if strcmp(P.xsect,'rect')
%             display('Solving with Y anti-symmetric and Z symmetric boundary condition');
%             P.obeveny = -1;
%             P.obevenz = 1;
%             asymYsymZ = RunNanobeamBands(P,datLoc);
% %             [asymYsymZ.midGap,asymYsymZ.gapSize] = findGaps(asymYsymZ);    %find gaps
%             
%             display('Solving with Y symmetric and Z antisymmetric boundary condition');
%             P.obeveny = 1;
%             P.obevenz = -1;
%             symYasymZ = RunNanobeamBands(P,datLoc);
% %             [symYasymZ.midGap,symYasymZ.gapSize] = findGaps(symYasymZ);    %find gaps
%         end
    end
    
%     close all;
    
    
    %%
    % write to data structure
    if P.savedat
        ds.P = P;
        ds.TE = TE;
        if P.solveTM == 1
            ds.TM = TM;
        end
        fname_mat = [fBase,'_bds.mat'];
        pathMat = [datLoc,fname_mat];
        save(pathMat,'ds');
    end
%     ds.P = P;
%     ds.TE = TE;
%     if P.solveTM == 1
%         ds.TM = TM;
%     end
    
    
    
else
    display('Data folder exists in working directory, loading...')
%     ds.sym = [];
%     ds.asym = [];
    dsFiles = dir([datLoc,fBase,'*_bds.mat']);
    load([datLoc,dsFiles.name]);
    TE = ds.TE;
    P = ds.P;
    if P.solveTM == 1
        TM = ds.TM;
    end
end

%% filter data below light line
% light line - here kx_norm runs from 0 to 1 and k from 0 to pi/a
c = 299792458;

lightline = c*TE.kx_norm/2/P.a; % factor of 2 such that kx_norm runs from 0 to 0.5 * (2*pi/a)
TE.F0 = TE.F;
TEbelow = TE.F < lightline;    % check which bands are below lightline
TE.F = TE.F.*TEbelow;        % filter out data below lightline
TE.F(TE.F==0) = NaN;      % replace zeros with NaN so they don't get plotted
[TE.midGap,TE.gapSize] = findGaps(TE);   %find gaps

if P.solveTM == 1
    TM.F0 = TM.F;
    TMbelow = TM.F < lightline;
    TM.F = TM.F.*TMbelow;        % filter out data below lightline
    TM.F(TM.F==0) = NaN;      % replace zeros with NaN so they don't get plotted
    [TM.midGap,TM.gapSize] = findGaps(TM);    %find gaps
end

%% plot bandstructure
if P.savebndplot == 1
    figure; 
    set(gcf,'position',[50 50 640 640],'units','pixels')
    hold on

    maxFreqs = [0 0 0 0];
    pHdl = [];
    % plot bands
    % call plot twice so legend appears properly
    p1 = plot(TE.kx_norm,TE.F(:,1)*1e-12,'-k','linewidth',3,'DisplayName','TE');
    plot(TE.kx_norm,TE.F*1e-12,'-k','linewidth',3);
    maxFreqs(1) = max(TE.F(:)*1e-12);
    pHdl(end+1) = p1;
    if P.solveTM == 1
        p2 = plot(TM.kx_norm,TM.F(:,1)*1e-12,'b:','linewidth',2,'DisplayName','TM');
        plot(TM.kx_norm,TM.F*1e-12,'b:','linewidth',2);
        maxFreqs(2) = max(TM.F(:)*1e-12);
        pHdl(end+1) = p2;
    end

    % plot bandgaps
    for j = 1:length(TE.gapSize)
        bgp = patch([0 1 1 0],1e-12.*(TE.midGap(j) + 0.5*[TE.gapSize(j) ...
            TE.gapSize(j) -TE.gapSize(j) -TE.gapSize(j)]),180/255*[1 1 1],'EdgeColor','none');
        alpha(bgp,0.5);
    end

    if P.solveTM == 1
        for j = 1:length(TM.gapSize)
            bgp = patch([0 1 1 0],1e-12.*(TM.midGap(j) + 0.5*[TM.gapSize(j) ...
                TM.gapSize(j) -TM.gapSize(j) -TM.gapSize(j)]),180/255*[0 0 1],'EdgeColor','none');
            alpha(bgp,0.1);
        end
    end
        
    % plot symmetric midgap frequencies
    for j = 1:length(TE.midGap)
        midfreqs = TE.midGap(j)*ones(length(TE.kx_norm),1);
%         plot(TE.kx_norm,midfreqs*1e-12,'.--r','linewidth',0.5);
    end
    
    %plot light line and cone
    kxNormAll = (0:P.kpts)/P.kpts;
    lightline = c*kxNormAll/2/P.a; % factor of 2 such that kx_norm runs from 0 to 0.5 * (2*pi/a)
    plot(kxNormAll,lightline*1e-12,'k-','Linewidth',2)
    patch([0 1 0],[0,lightline(end)*1e-12,lightline(end)*1e-12],'k','FaceAlpha',0.15,'LineStyle','none');
    
    

    

    legend(pHdl,'location','southoutside','orientation','horizontal');
    xlabel('k_x / (\pi/a)','FontSize',12);
    ylabel('Frequency (THz)','FontSize',12);
    maxY = 260;%max(maxFreqs);
    axis([0 1 100 maxY]);
    set(gca,'XTick',[0; 1]);
%         set(gca,'XTickLabel',{'G','C'},'fontname','symbol','fontsize',16)
    set(gca,'XTickLabel',{'\Gamma','X'},'fontsize',12)

    if strcmp(P.xsect,'tri')
        thStr = ['\theta = ',num2str(P.theta),'^o, '];
    elseif strcmp(P.xsect,'rect')
        thStr = ['th = ',num2str(P.th*1e9,'%.0f'),'nm, '];
    end
    bandtitle = {[thStr,...
                  'a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                  'w = ',num2str(P.w*1e9,'%.0f'),'nm, '];
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
end

% save data structures
if P.savedat
    ds.P = P;
    ds.TE = TE;
    if P.solveTM == 1
        ds.TM = TM;
    end
    fname_mat = [fBase,'_bds.mat'];
    pathMat = [datLoc,fname_mat];
    save(pathMat,'ds');
end

% rmpath(utilFilesPath)

