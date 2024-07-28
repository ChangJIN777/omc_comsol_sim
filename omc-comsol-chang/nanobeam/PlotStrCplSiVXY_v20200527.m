% function to generate SiV strain coupling plots
% inputs:
% ds: data structure from simulation
% model: COMSOL model with mechanical simulations
% mModes: solution no. of mechanical eigenmodes to plot for
% xSlc, ySlc,zSlc: x-, y-, z-coordinate to plot slice at
% zSiV = [1 1 1], [-1 1 1], [-1 -1 1], [1 -1 1]:
%   z-orientation of SiV


function PlotStrCplSiVXY_v20200527(ds,model,mModes,pSlcAll,zSiVAll,datLoc)
%% front end matters
P = ds.P;
mfem = ds.mfem;
cpl = ds.cpl;

% Geometry parameters
hx = P.geom(:,1);
hy = P.geom(:,2);
xpos = P.geom(:,3);
ypos = P.geom(:,4);
w = P.geom(:,5);
a = P.geom(:,6);

% crystal orientation and strain susceptibilities
rxtal = P.rxtal; 
dg = P.dg;
fg = P.fg;

%% Filenames
if ~isfield(P,'fileBase')
    P = CreateFileBase(P);
    ds.P = P;
    fileBase = P.fileBase;
else
    fileBase = P.fileBase;
end

%% Plot over different modes
% mi/mIdx = absolute/relative indices of mode
% e.g. if mModes = (1,3,5), then mModes(mIdx) = mi
mIdx = 1;
for mi = mModes
    
    
    %% XY: Compile strain coupling profiles for plotting
    % get strain coupling
    
    for oi = 1:length(zSiVAll)
        zSiV = zSiVAll{oi};
        
        % for non-standard orientations, convert to row vector and invert signs
        [nr,nc] = size(zSiV);
        if nr==3 && nc==1
            zSiV = zSiV';
        end
        if (prod(zSiV==[-1 -1 -1]) || prod(zSiV==[ 1 -1 -1]) || ...
            prod(zSiV==[ 1  1 -1]) || prod(zSiV==[-1  1 -1]))
            zSiV = -1*zSiV;
        end

        % get string representation of orientation
        % for reference to subfield in data struct
        zSiVstr = num2str(zSiV);
        zSiVstr = strrep(zSiVstr,'-','m');
        zSiVstr = strrep(zSiVstr,' ','');
        
        %% get coordinates for strain interpolation
        xSlc = pSlcAll.(['z',zSiVstr]).pos(1,mIdx);
        ySlc = pSlcAll.(['z',zSiVstr]).pos(2,mIdx);
        zSlc = pSlcAll.(['z',zSiVstr]).pos(3,mIdx);
                
        dl = 10e-9;  % resolution interpolation

        % determine full extents of beam - especially important for triangular
        % cross-section, where width of beam is function of depth under surface
        if strcmp(P.xsect,'tri')
            wPlot = (zSlc+P.th/2)/P.th*P.w; 
        elseif strcmp(P.xsect,'rect')
            wPlot = max(w);
        end

        if isfield(P,'asymCav') && P.asymCav
            beamLen = P.beamLen;
        else
            beamLen = P.beamLenHalf;
        end

        if ~abs(P.mevenz)
            ZMin = -P.th/2;
        else
            ZMin = 0;
        end

        % generate coordinate arrays
        XDat_XY = [0:5e-9:a/2,a/2:dl:beamLen,beamLen];
        YDat_XY = [0:dl:wPlot/2,wPlot/2];
        YDat_YZ = [0:dl:max(w)/2,max(w)/2];
        ZDat_YZ = [ZMin:dl:P.th/2,P.th/2];
        
        % insert custom x/y/z coordinates
        if ~isempty(xSlc) && isempty(find(abs(XDat_XY-xSlc)<1e-9,1))
            XDat_XY = sort([XDat_XY,xSlc]);
        end
        if ~isempty(ySlc) && isempty(find(abs(YDat_XY-ySlc)<1e-9,1))
            YDat_XY = sort([YDat_XY,ySlc]);
            YDat_YZ = sort([YDat_YZ,ySlc]);
        end
        if ~isempty(zSlc) && isempty(find(abs(ZDat_YZ-zSlc)<1e-9,1))
            ZDat_YZ = sort([ZDat_YZ,zSlc]);
        end
        
        
        [~,~,dsSiV] = CalcStrCplSiV_v20200527(ds,model,mi,XDat_XY,YDat_XY,zSlc,...
                                {zSiV},[],0);
        
        %% assemble coordinates and strain coupling arrays
        % find sector where strain coupling is maximum for this orientation
        % of SiV
        secsgn = cpl.(['SiV_',zSiVstr]).full.LSiVMaxSec(:,mi);
        xSlc = xSlc*secsgn(1);
        
        % identify sectors to plot
        if zSlc < 0 && abs(P.mevenz)
            zSec = 'm';
        else
            zSec = 'p';
        end

        allSec = {'pp'};
        if abs(P.mevenx) && P.meveny == 0
            allSec = [allSec,'mp'];
        elseif abs(P.meveny) && P.mevenx == 0
            allSec = [allSec,'pm'];
        elseif abs(P.meveny) && abs(P.mevenx)
            allSec = [allSec,'mp','pm','mm'];
        end

        % assemble coordinates and reshape arrays in relevant sectors
        for j = 1:length(allSec)
            secLbl = [allSec{j},zSec];
            
            for k = 1:6
                LDat0 = dsSiV.(['z',zSiVstr]).(secLbl).LSiV;
                LDat = reshape(LDat0(mi,:),[length(YDat_XY),length(XDat_XY)]);

                % flip arrays for plotting
                if strcmp(secLbl(1),'m')    % -x axis
                    LDat = fliplr(LDat);
                    XPlt = -1*fliplr(XDat_XY);
                else
                    XPlt = XDat_XY;
                end
                if strcmp(secLbl(2),'m')    % -y axis
                    LDat = flipud(LDat);
                    YPlt = -1*fliplr(YDat_XY);
                else
                    YPlt = YDat_XY;
                end
                ePlotXY.(allSec{j}).XDat = XPlt;
                ePlotXY.(allSec{j}).YDat = YPlt;
                ePlotXY.(allSec{j}).LDat = LDat;
            end
        end

        % if plot axis contain zero coordinate, eliminate row/col corresponding
        % to that coordinate
        xid = 0;
        if ~isempty(find(XDat_XY==0,1))
            xid = 1;
        end
        yid = 0;
        if ~isempty(find(YDat_XY==0,1))
            yid = 1;
        end

        x_pp = ePlotXY.pp.XDat;
        y_pp = ePlotXY.pp.YDat;
        L_pp = ePlotXY.pp.LDat;
        ePlotXY.full = ePlotXY.pp;  % if no x & y symmetries


        % combine strain coupling data from all sectors
        if abs(P.mevenx)
            L_mp = ePlotXY.mp.LDat;
            x_mp = ePlotXY.mp.XDat;
            ePlotXY.full.LDat = [L_mp(:,1:end-xid),L_pp];
            ePlotXY.full.XDat = [x_mp(1:end-xid),x_pp];
        end
        if abs(P.meveny)
            L_pm = ePlotXY.pm.LDat;
            y_pm = ePlotXY.pm.YDat;
            ePlotXY.full.LDat = [L_pm(1:end-yid,:);L_pp];
            ePlotXY.full.YDat = [y_pm(1:end-yid),y_pp];
        end
        if abs(P.meveny) && abs(P.mevenx)
            L_mm = ePlotXY.mm.LDat;
            ePlotXY.full.LDat = [L_mm(1:end-yid,1:end-xid),L_pm(1:end-yid,:);...
                                 L_mp(:,1:end-xid),      L_pp];
        end
        
        % get max strain coupling and position
        XPlotXY = ePlotXY.full.XDat;
        YPlotXY = ePlotXY.full.YDat;
        LPlotXY = ePlotXY.full.LDat;

        % get max strain coupling
        [LmaxXY, ImaxXY] = max(LPlotXY(:));
        [ImaxY, ImaxX] = ind2sub(size(LPlotXY),ImaxXY);
        maxX = XPlotXY(ImaxX);
        maxY = YPlotXY(ImaxY);
        ePlotXY.full.LmaxXY = LmaxXY;
        ePlotXY.full.maxX = maxX;
        ePlotXY.full.maxY = maxY;
        
        % get slice along X at custom Y for plotting
        YSlice = ySlc*secsgn(2);
        [~,IslcY] = min(abs(YPlotXY-YSlice));
        LXSlice = LPlotXY(IslcY,:);
        
        %% XY surface plot
        figure;
        set(gcf,'position',[100 100 913 500])
        axpos = [0.07 0.55 0.82 0.250];
        ax = axes('position',axpos); hold(ax,'on')
        box on
        axis off

        % strain coupling
        surf(XPlotXY,YPlotXY,LPlotXY);
        grid('off'); shading('interp'); 
        daspect([1 1 1]);
        colormap(colortable('Rainbow'))
        cb = colorbar('eastoutside');
        cbp = get(cb,'Position');
        set(cb,'Position',[0.9 axpos(2) 0.03 0.25]);

        % geometry overlay
        gXL = [XPlotXY(1) XPlotXY(end)];
        gXU = gXL;
        gYL = 1/2*[-P.w -P.w];
        gYU = 1/2*[P.w P.w]; 

        if isfield(P,'celltype') && strcmp(P.celltype,'blockTet')
            % for block-tether geometry
            xdat = [];
            ydatu = [];
            for j=1:length(hy)
                xdat = [xdat,   xpos(j)-a(j)/2,...
                                xpos(j)-hx(j)/2,...
                                xpos(j)-hx(j)/2,...
                                xpos(j)+hx(j)/2,...
                                xpos(j)+hx(j)/2,...
                                xpos(j)+a(j)/2];
                ydatu = [ydatu, hy(j)/2,...
                                hy(j)/2,...
                                w(j)/2,...
                                w(j)/2,...
                                hy(j)/2,...
                                hy(j)/2];
            end
            ydatd = -ydatu;
            plot(ax,xdat,ydatu,'k:','linewidth',0.5)
            plot(ax,xdat,ydatd,'k:','linewidth',0.5)
        else
            % for beam-hole geometry
            plot(ax,gXU,gYU,'k:','linewidth',0.5)
            plot(ax,gXL,gYL,'k:','linewidth',0.5)
            plot(ax,gXL(1)*[1 1],[gYL(1) gYU(1)],'k:','linewidth',0.5)
            plot(ax,gXL(2)*[1 1],[gYL(2) gYU(2)],'k:','linewidth',0.5)
    %         plot(ax,gXL,[0 0],'r:','linewidth',0.5)
    %             xlim(ax,[min(gXU) max(gXU)])

            %airhole overlay
            for j=1:length(hy)
                xdat = linspace(-hx(j)/2,hx(j)/2,100)';
                ydatu(:,j) = sqrt(1-linspace(-hx(j)/2,hx(j)/2,100).^2/(hx(j)/2)^2)*hy(j)/2 + ypos(j);
                ydatd(:,j) = -sqrt(1-linspace(-hx(j)/2,hx(j)/2,100).^2/(hx(j)/2)^2)*hy(j)/2 + ypos(j);
                plot(ax,[xdat + xpos(j),xdat + xpos(j)], ...
                    [ydatu(:,j),ydatd(:,j)],'k:','linewidth',0.5);
            end
        end
        plot(ax,gXL,YSlice*[1 1],'r:','linewidth',0.5)
        xlim(ax,[min(gXU) max(gXU)])

        % settings for overlay to appear on top
        view(0,-90); 
        set(gca,'YDir','Reverse')

        % slice plot
        bx = axes('position',axpos-[0 0.3 0 0]); hold(bx,'on')
        box on
        axis on
        plot(bx,XPlotXY,LXSlice,'k','linewidth',1.5)
        xlim(bx,[min(gXU) max(gXU)])
        ylim(bx,[0 max(LXSlice)])
        for j=1:length(hy)
            plot(bx,[xpos(j),xpos(j)],[0 max(LXSlice)],'--','Color',[0.5 0.5 0.5],'linewidth',0.5)
        end
        title(['Cut along X-axis at y = ',num2str(YSlice*1e9,'%.1f'),...
               'nm (red dotted horizontal line)'],'fontname','arial','fontsize',10)

        % clear x ticks and labels
        xticks('')
        xticklabels({});

        % set y ticks and convert to MHz
        yTickInt = round(max(LXSlice)/5,1,'significant');
        yticks((0:5)*yTickInt)
        bYTick = get(bx,'ytick');
        bYTickLbl = {};
        for lbi = 1:length(bYTick)
            bYTickLbl{end+1} = num2str(bYTick(lbi)/1e6,'%.1f');
        end
        set(bx,'yticklabel',bYTickLbl)
        ylabel('\lambda_{SiV} / MHz')
        xlabel('x / um')
        linkaxes([bx,ax],'x')

        %% Filename and title settings
        mfilestr = ['wM_',num2str(mfem.freqs(mi)*1e-9,'%.2f'),'GHz_',...
                    'lSiV_',num2str(cpl.(['SiV_',zSiVstr]).full.LSiVMax(mi)*1e-6,'%.2f'),'MHz'];

        mtitlestr = ['\omega_M = ',num2str(mfem.freqs(mi)*1e-9,'%.2f'),'GHz, ',...
                     '\lambda_{SiV} = ',num2str(cpl.(['SiV_',zSiVstr]).full.LSiVMax(mi)*1e-6,'%.2f'),'MHz'];

        ctitlestr = '';
        if isfield(ds,'cpl')
            if isfield(cpl,'gOM')
                ctitlestr = ['g_{OM} = ',num2str(abs(real(cpl.gOM(mi)))*1e-3,'%.0f'),'kHz'];
            end
            if isfield(cpl,'F')
                ctitlestr = [ctitlestr,', ',...
                             'F = ',num2str(cpl.F(mi),'%.2e')];
            end
        else
            ctitlestr = '';
        end

        %%
        % plot title
    %     mtitlestr = [mtitlestr0,', \lambda_{SiV} = ',...
    %                  num2str(cpl.(['SiV_',zSiVstr]).full.LSiVMax*1e-6,'%.2f'),'MHz'];
    %     mfilestr = [mfilestr0,'_lSiV_',...
    %                  num2str(cpl.(['SiV_',zSiVstr]).full.xposSiVMax*1e-6,'%.2f'),'MHz'];
        plotTitle = {['Strain coupling, SiV axis along [',num2str(zSiV),'], ',...
                    mtitlestr];...
                    ctitlestr;...
                    ['Cut along XY-plane at depth = ',num2str((P.th/2-zSlc)*1e9,'%.1f'),'nm: ',...
                     'max \lambda_{SiV} in plane = ',num2str(LmaxXY*1e-6,'%.2f'),'MHz ',...
                     'at (x,y) = (',num2str((maxX-P.xc)*1e9,'%.1f'),...
                     ',',num2str(maxY*1e9,'%.1f'),...
                     ') nm'];''};
        title(ax,plotTitle,'fontname','arial','fontsize',10)

        % save file
        pathFig = [datLoc,fileBase,'_',mfilestr,'_',...
                   'd_',num2str((P.th/2-zSlc)*1e9,'%.1f'),'nm_',...
                   zSiVstr,'_XY'];
        saveas(gcf,[pathFig,'.png']);
        saveas(gcf,[pathFig,'.fig']);
        close;
        
        
        
        
        
        
        
        
        
        
        
        
    end
    
    
    
    
    
    mIdx = mIdx + 1;
end























end