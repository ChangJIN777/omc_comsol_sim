% function to generate SiV strain coupling plots
% inputs:
% ds: data structure from simulation
% model: COMSOL model with mechanical simulations
% mModes: solution no. of mechanical eigenmodes to plot for
% xSlc, ySlc,zSlc: x-, y-, z-coordinate to plot slice at
% zSiV = [1 1 1], [-1 1 1], [-1 -1 1], [1 -1 1]:
%   z-orientation of SiV


function varargout = PlotStrCplSiVYZ(ds,model,mModes,pSlcAll,zSiVAll,datLoc)
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

% generate geometry signs
[XG,YG,ZG] = meshgrid([1 -1],[1 -1],[1 -1]);
secSgnAll = transpose([XG(:),YG(:),ZG(:)]);

% generate symmetries for all sectors
[XV,YV,ZV] = meshgrid([1 P.mevenx],[1 P.meveny],[1 P.mevenz]);
secSymAll = transpose([XV(:),YV(:),ZV(:)]); 
secIdxs = find(prod(secSymAll)~=0); % find sectors where symmetries apply

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

for mIdx = 1:length(mModes) % mi = mModes
    mi = mModes(mIdx);
    
    %% YZ: Compile strain coupling profiles for plotting
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
        % pSlc is now array of size 3 x length(mModes) x length(zSiVAll)
        xSlc = pSlcAll(1,mIdx,oi);
        ySlc = pSlcAll(2,mIdx,oi);
        zSlc = pSlcAll(3,mIdx,oi);
        
        % make positive if there is symmetry plane, since strain is only
        % interpolated in symmetry sector
        xInt = xSlc*(sign(xSlc))^(abs(P.mevenx));
        yInt = ySlc*(sign(ySlc))^(abs(P.meveny));
        zInt = zSlc*(sign(zSlc))^(abs(P.mevenz));
        
        % determine full extents of beam - especially important for triangular
        % cross-section, where width of beam is function of depth under surface
        if strcmp(P.xsect,'tri')
            wPlot = (zSlc+P.th/2)/P.th*P.w; 
        elseif strcmp(P.xsect,'rect')
            wPlot = max(w);
        elseif strcmp(P.xsect,'isoFit')
            wPlot = min(max(P.w),2*(zSlc + P.th/2)*sqrt(3));
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
        dl = 10e-9;  % resolution interpolation
        XDat_XY = [0:5e-9:a/2,a/2:dl:beamLen,beamLen];
        YDat_XY = [0:dl:wPlot/2,wPlot/2];
        YDat_YZ = [0:dl:max(w)/2,max(w)/2];
        ZDat_YZ = [ZMin:dl:P.th/2,P.th/2];

        % insert custom x/y/z coordinates if no overlap with regular grid
        if ~isempty(xInt) && isempty(find(abs(XDat_XY-xInt)<1e-9,1))
            XDat_XY = sort([XDat_XY,xInt]);
        end
        if ~isempty(yInt) && isempty(find(abs(YDat_XY-yInt)<1e-9,1))
            YDat_XY = sort([YDat_XY,yInt]);
            YDat_YZ = sort([YDat_YZ,yInt]);
        end
        if ~isempty(zInt) && isempty(find(abs(ZDat_YZ-zInt)<1e-9,1))
            ZDat_YZ = sort([ZDat_YZ,zInt]);
        end
        
        [~,~,dsSiV] = CalcStrCplSiV(ds,model,mi,xInt,YDat_YZ,ZDat_YZ,...
                                {zSiV},[],0,[xInt,yInt,zInt]);
        
        %% assemble coordinates and strain coupling arrays
        % find sector where strain coupling is maximum for this orientation
        % of SiV
        maxSec = dsSiV.LSiVMaxSec(1,1,1,1);
        maxSecSgn = secSgnAll(:,maxSec);
        
        % identify sectors to plot depending on x-coordinate
        secsPlt = find((secSgnAll(1,secIdxs))==maxSecSgn(1));
        
        % extract strain coupling and coordinate arrays
        % and apply transforms depending on sector
        LDatAll = zeros(length(ZDat_YZ),length(YDat_YZ),length(secsPlt));
        ZDatAll = zeros(length(ZDat_YZ),length(secsPlt));
        YDatAll = zeros(length(YDat_YZ),length(secsPlt));
        
        for sIdx = 1:length(secsPlt)
            si = secsPlt(sIdx);
            
            % extract and reshape arrays
            LDat0 = dsSiV.LSiV(1,:,1,1,si);
            LDatS = reshape(LDat0,[length(YDat_YZ),length(ZDat_YZ)]);
            LDatS = LDatS'; % transpose array so geom Y/Z-axis --> plot X/Y-axis
            
            % flip and/or change signs of arrays
            zgt0 = (secSgnAll(3,si) == 1);
            ygt0 = (secSgnAll(2,si) == 1);
            
            ZDatS = zgt0*ZDat_YZ + ~zgt0*-1*fliplr(ZDat_YZ);
            YDatS = ygt0*YDat_YZ + ~ygt0*-1*fliplr(YDat_YZ);
            
            if ~zgt0
                LDatS = flipud(LDatS);
            end
            if ~ygt0
                LDatS = fliplr(LDatS);
            end
            
            % assemble arrays
            ZDatAll(1:length(ZDat_YZ),sIdx) = ZDatS(1:length(ZDat_YZ));
            YDatAll(1:length(YDat_YZ),sIdx) = YDatS(1:length(YDat_YZ));
            LDatAll(1:length(ZDat_YZ),1:length(YDat_YZ),sIdx) = ...
                LDatS(1:length(ZDat_YZ),1:length(YDat_YZ));
        end
        
        % get max strain coupling and position
        [LmaxYZ, ImaxYZ] = max(LDatAll(:));
        [ImaxZ, ImaxY, ImaxS] = ind2sub(size(LDatAll),ImaxYZ);
        maxZ = ZDatAll(ImaxZ,ImaxS);
        maxY = YDatAll(ImaxY,ImaxS);
        
        % get slice along Z at custom Y for plotting
        ySlcSecSgn = (sign(ySlc))^(abs(P.meveny)) + (sign(ySlc)==0);
        secsYSlc = find((secSgnAll(2,:))==ySlcSecSgn);
        [secsYSlc,secsYSlcIdx] = intersect(secsPlt,secsYSlc);
        [~,IslcY] = min(abs(YDatAll(:,secsYSlcIdx(1))-ySlc));
        
        % get slice along Y at custom Z for plotting
        zSlcSecSgn = (sign(zSlc))^(abs(P.mevenz)) + (sign(zSlc)==0);
        secsZSlc = find((secSgnAll(3,:))==zSlcSecSgn);
        [secsZSlc,secsZSlcIdx] = intersect(secsPlt,secsZSlc);
        IslcZ = min(abs(ZDatAll(:,secsZSlcIdx(1))-zSlc));
        
        % get max coupling at specified coordinates in plane
        IslcSec = intersect(secsYSlc,secsZSlc);
        IslcSecIdx = find(secsPlt == IslcSec);
        LPmax = LDatAll(IslcZ,IslcY,IslcSecIdx);
        
        %% find strain tensor at specified coordinates
%         eSiVVoigt = dsSiV.eSiVpos(1:6,1,1,1,maxSec);
        eSiVVoigt = dsSiV.eSiVpos(1:6,1,1,1,IslcSec);
        eSiVVoigt = squeeze(eSiVVoigt);
        eSiVcoord = [eSiVVoigt(1), eSiVVoigt(6), eSiVVoigt(5);
                     eSiVVoigt(6), eSiVVoigt(2), eSiVVoigt(4);
                     eSiVVoigt(5), eSiVVoigt(4), eSiVVoigt(3)];
        
        %% YZ surface plot
        figure;
        set(gcf,'position',[300 200 800 600])

        % YZ geometry overlay
        gxpos = [100 400 600 200];
        gx = axes('position',gxpos); hold(gx,'on')
        set(gx,'Units','pixels','Position',gxpos,'ActivePositionProperty','position');
        box on
        axis off
        daspect([1 1 1]);
        if isfield(P,'asymCav') && P.asymCav
            beamLen = P.beamLen;
        else
            beamLen = P.beamLenHalf;
        end
        XPDat = [0,beamLen];
        if abs(P.mevenx)
            XPlot = [-1*fliplr(XPDat(2:end)) XPDat];
        else
            XPlot = XPDat;
        end
        gXL = [XPlot(1) XPlot(end)]*1e9;
        gXU = gXL;
        gYL = -1/2*max(w)*[1 1]*1e9;
        gYU = 1/2*max(w)*[1 1]*1e9;

        if isfield(P,'celltype') && strcmp(P.celltype,'blockTet')
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
            plot(gx,xdat*1e9,ydatu*1e9,'k-','linewidth',1)
            plot(gx,xdat*1e9,ydatd*1e9,'k-','linewidth',1)
        else
            plot(gx,gXU,gYU,'k-','linewidth',1)
            plot(gx,gXL,gYL,'k-','linewidth',1)
            plot(gx,gXL(1)*[1 1],[gYL(1) gYU(1)],'k-','linewidth',1)
            plot(gx,gXL(2)*[1 1],[gYL(2) gYU(2)],'k-','linewidth',1)

            %airhole overlay
            for j=1:length(hy)
                xdat = linspace(-hx(j)/2,hx(j)/2,100)';
                ydatu(:,j) = sqrt(1-linspace(-hx(j)/2,hx(j)/2,100).^2/(hx(j)/2)^2)*hy(j)/2 + ypos(j);
                ydatd(:,j) = -sqrt(1-linspace(-hx(j)/2,hx(j)/2,100).^2/(hx(j)/2)^2)*hy(j)/2 + ypos(j);
                plot(gx,[xdat + xpos(j),xdat + xpos(j)]*1e9, ...
                    [ydatu(:,j),ydatd(:,j)]*1e9,'k-','linewidth',1);
            end

        end
        plot(gx,xSlc*1e9*[1 1],[gYL(1) gYU(1)],'r-','linewidth',1)
        xlim(gx,[min(gXU) max(gXU)])

        % settings for overlay to appear on top
        view(0,-90); 
        set(gca,'YDir','Reverse')
    %     xlabel('x/nm'); ylabel('y/nm')

        % YZ strain coupling
        axpos0 = [175 80 325 200];
        ax = axes('Units','pixels','position',axpos0); hold(ax,'on')
        box on
        axis off
        hold on
        for sIdx = 1:length(secsPlt)%si = secsPlt
            surf(YDatAll(:,sIdx),ZDatAll(:,sIdx),LDatAll(:,:,sIdx))
        end
%         surf(YPlotYZ,ZPlotYZ,LPlotYZ);
        grid('off'); shading('interp'); 
        daspect([1 1 1]);
        cmax = cpl.SiV.LSiVMaxMode(mi);
        caxis([0 0.5*cmax])
        cmap = colormap(hot);
        colormap(cmap(1:56,:));
        cb = colorbar('southoutside');
        cbp = get(cb,'Position');
        set(ax,'Units','pixels','Position',axpos0,'ActivePositionProperty','position');

        % YZ cross-section overlay
        plot(ax,[-P.w/2 P.w/2],[P.th/2 P.th/2],'k','linewidth',1)
        if strcmp(P.xsect,'tri')
            plot(ax,[P.w/2 0],[P.th/2 -P.th/2],'k','linewidth',1)
            plot(ax,[0 -P.w/2],[-P.th/2 P.th/2],'k','linewidth',1)
        elseif strcmp(P.xsect,'rect')
            plot(ax,[-max(w)/2 max(w)/2],[-P.th/2 -P.th/2],'k','linewidth',1)
            plot(ax,[-max(w)/2 -max(w)/2],[-P.th/2 P.th/2],'k','linewidth',1)
            plot(ax,[max(w)/2 max(w)/2],[-P.th/2 P.th/2],'k','linewidth',1)
        end
        xlim(ax,[-max(w)/2 max(w)/2])
        ylim(ax,[-P.th/2 P.th/2])

        % plot specified location
        plot([-max(w)/2 max(w)/2],zSlc*[1 1],'r--','linewidth',1)
        plot(ySlc*[1 1],[-P.th/2 P.th/2],'r--','linewidth',1)

        % settings for overlay to appear on top
        view(0,-90); 
        set(gca,'YDir','Reverse')

        % get axes position to size slice plots
        axpos = get(ax,'Position');

        axAR = get(ax,'PlotBoxAspectRatio');
        axARxy = axAR(1)/axAR(2);   %final aspect ratio
        axARinit = axpos(3)/axpos(4); %initial aspect ratio
        if axARxy < axARinit
            % height fills up plot axes - rescale width
            YSPlotW = axARxy/axARinit*axpos(3);
            YSPlotDx = 0.5*abs(axpos(3)-YSPlotW);
            ZSPlotH = axpos(4);
            ZSPlotDy = 0;
        else
            % width fills up plot axes - rescale height
            YSPlotW = axpos(3);
            YSPlotDx = 0;
            ZSPlotH = axARinit/axARxy*axpos(4);            
            ZSPlotDy = 0.5*abs(axpos(4)-ZSPlotH);
        end

        set(cb,'Units','pixels','Position',[axpos(1)+YSPlotDx axpos(2)-35 YSPlotW 25]);%,'Units','pixels');


        % slice plot along y
    %         ZMaxIdx = find(abs(ZPlot-LZMaxCoord(mi,3))<interpRes/2.5);
    %         LYSlice = LPlot(ZMaxIdx,:);
        YSPlotH = 75;
        bx = axes('Units','pixels','position',[axpos(1)+YSPlotDx axpos(2)+axpos(4)+50-ZSPlotDy YSPlotW YSPlotH]); hold(bx,'on')
        box on
        hold on
        for sIdx = secsZSlcIdx'%1:length(secsZSlc)
            plot(bx,YDatAll(:,sIdx),LDatAll(IslcZ,:,sIdx),'k','linewidth',1.5)
        end
%         plot(bx,YPlotYZ,LYSlice,'k','linewidth',1.5)
        xlim(bx,[-max(w)/2 max(w)/2])
        LYSlice = LDatAll(IslcZ,:,secsZSlcIdx);
        ylim(bx,[0 max(LYSlice(:))])

        YSYTick = get(bx,'ytick');
        YSYTickLbl = {};
        for lbi = 1:length(YSYTick)
            YSYTickLbl{end+1} = num2str(YSYTick(lbi)/1e6,'%.1f');
        end
        set(bx,'yticklabel',YSYTickLbl)
        XSYTick = get(bx,'xtick');
        XSYTickLbl = {};
        for lbi = 1:length(XSYTick)
            XSYTickLbl{end+1} = num2str(XSYTick(lbi)*1e9,'%.0f');
        end
        set(bx,'xticklabel',XSYTickLbl)
        ylabel('\lambda_{SiV} / MHz')
        xlabel('y / nm')
    %         title({['Slice along Y-axis at depth = ',num2str((P.th/2-LZMaxCoord(mi,3))*1e9,'%.1f'),'nm: '];...
    %                      ['max \lambda_{SiV} = ',num2str(LZMax(mi)*1e-6,'%.2f'),'MHz ',...
    %                      'at y = ',num2str(LZMaxCoord(mi,2)*1e9,'%.1f'),...
    %                      ' nm']},'fontname','arial','fontsize',10)
        linkaxes([bx,ax],'x')

        % slice plot along z
        ZSPlotW = 75;
        cx = axes('Units','pixels','position',[axpos(1)-ZSPlotW-20+YSPlotDx axpos(2)+ZSPlotDy ZSPlotW ZSPlotH]); hold(cx,'on')
        box on
        hold on
        for sIdx = secsYSlcIdx'%1:length(secsYSlc)
            plot(cx,LDatAll(:,IslcY,sIdx),(ZDatAll(:,sIdx)),'k','linewidth',1.5)
        end
        LZSlice = LDatAll(:,IslcY,secsYSlcIdx);
        ylim(cx,[-P.th/2 P.th/2])
        if ~isnan(max(LZSlice(:)))
            xlim(cx,[0 max(LZSlice(:))])
        end

        % set x-axis labels (strain coupling)
        ZSXTick = get(cx,'xtick');
        ZSXTickLbl = {};
        for lbi = 1:length(ZSXTick)
            ZSXTickLbl{end+1} = num2str(ZSXTick(lbi)/1e6,'%.1f');
        end
        set(cx,'xticklabel',ZSXTickLbl)

        % set y-axis labels (depth)
        YSXTickPos = P.th/2 - fliplr((0:50:P.th*1e9)*1e-9);
        YSXTickLbs = fliplr((0:50:P.th*1e9));
        set(cx,'ytick',YSXTickPos)
        set(cx,'yticklabel',YSXTickLbs)
        
    %     YSXTick = get(cx,'ytick');
    %     YSXTickLbl = {};
    %     for lbi = 1:length(YSXTick)
    %         YSXTickLbl{end+1} = num2str((P.th/2-YSXTick(lbi))*1e9,'%.0f');
    %     end
    %     set(cx,'yticklabel',YSXTickLbl)
        xlabel('\lambda_{SiV} / MHz')
        ylabel('depth / nm')

        %% plot strain tensor components
        % assemble full strain tensor at specified coordinate
        %% check assignment of eij's
        
        T = eSiVcoord;
%         T = dsSiV.(['z',zSiVstr]).full.eSiVMax(:,:,mi);
        % plot settings
        TPlotW = 175;
        TPlotH = axpos0(4)+YSPlotH-25;
        dx = axes('Units','pixels','position',[axpos(1)+axpos(3)+75 axpos(2)+50 TPlotW TPlotH]); hold(dx,'on')
        box on; grid on
        b = bar3(dx,real(T));
    %         set(dx,'ZAxisLocation','origin')
        view([-25 25])
        for k = 1:length(b)
            zdata = get(b(k),'ZData');
            set(b(k),'CData', zdata);
            set(b(k),'FaceColor', 'interp');
        end

        cmax = max(abs(real(T(:))));
        % colormap(bwr(101)); 
        try
            caxis([-cmax cmax])
        catch
            caxis('auto')
        end
        
        colormap(dx,bwr)
        cb = colorbar(dx,'location','southoutside');
        cbp = get(cb,'Position');
        set(cb,'Units','pixels','Position',[axpos(1)+axpos(3)+75 axpos(2)-35 TPlotW 25]);
    %         set(cb,'Units','pixels','Position',[axpos(1)+YSPlotDx axpos(2)-35 YSPlotW 25]);%,'Units','pixels');
    %         [axpos(1) axpos(2)-25 axpos(3) 25]

        set(gca,'XTickLabel',{'SiV_x','SiV_y','SiV_z'});
        set(gca,'YTickLabel',{'SiV_x','SiV_y','SiV_z'});
        title({'SiV strain tensor components';'at max coupling in cut plane';''},'fontname','arial','fontsize',10)


        %% Filename and title settings
        mfilestr = ['wM_',num2str(mfem.freqs(mi)*1e-9,'%.2f'),'GHz_',...
                    'lSiV_',num2str(cpl.SiV.LSiVMaxMode(mi)*1e-6,'%.2f'),'MHz'];

        mtitlestr = ['\omega_M = ',num2str(mfem.freqs(mi)*1e-9,'%.2f'),'GHz, ',...
                     '\lambda_{SiV} = ',num2str(cpl.SiV.LSiVMaxMode(mi)*1e-6,'%.2f'),'MHz'];

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
        % full plot title
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.5, 1,{['Strain coupling, SiV axis along [',num2str(zSiV),'], ',...
                    mtitlestr];...
                    ctitlestr;...
                    ['Cut along YZ-plane at x = ',num2str((xSlc-P.xc)*1e9,'%.1f'),'nm: ',...
                     'at (y,d) = (',num2str(ySlc*1e9,'%.1f'),...
                     ',',num2str((P.th/2-zSlc)*1e9,'%.1f'),')nm, ',...
                    '\lambda_{SiV} = ',num2str(LPmax*1e-6,'%.2f'),'MHz']},...
                'HorizontalAlignment' ,'center','VerticalAlignment', 'top','fontname','arial','fontsize',10,'FontWeight','bold')


        % save file
        pathFig = [datLoc,fileBase,'_',mfilestr,'_',...
                   'x_',num2str((xSlc-P.xc)*1e9,'%.1f'),'nm_',...
                   zSiVstr,'_YZ'];
        saveas(gcf,[pathFig,'.png']);
        saveas(gcf,[pathFig,'.fig']);
        close;
        clear dsSiV
        
        dsSiVplt.LSiVSlc(mIdx,oi) = LPmax;
        
        
        
        
        
        
        
        
        
        
    end
    dsSiVplt.LSiVSlcModes(mIdx) = mi;
    
    
    
    
%     mIdx = mIdx + 1;
end

if nargout >= 1
    varargout{1} = dsSiVplt;
end





















end