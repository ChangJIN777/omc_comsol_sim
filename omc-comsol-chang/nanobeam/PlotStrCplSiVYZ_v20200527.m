% function to generate SiV strain coupling plots
% inputs:
% ds: data structure from simulation
% model: COMSOL model with mechanical simulations
% mModes: solution no. of mechanical eigenmodes to plot for
% xSlc, ySlc,zSlc: x-, y-, z-coordinate to plot slice at
% zSiV = [1 1 1], [-1 1 1], [-1 -1 1], [1 -1 1]:
%   z-orientation of SiV


function PlotStrCplSiVYZ_v20200527(ds,model,mModes,pSlcAll,zSiVAll,datLoc)
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
%         XDat_XY = linspace(0,beamLen,round(beamLen/interpRes)+1);
%         YDat_XY = linspace(0,wPlot/2,round(wPlot/2/interpRes)+1);
%         YDat_YZ = linspace(0,max(P.w_hole)/2,round(max(P.w_hole)/2/interpRes)+1);
%         ZDat_YZ = linspace(ZMin,P.th/2,round((P.th/2-ZMin)/interpRes)+1);

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
        
        
        [~,~,dsSiV] = CalcStrCplSiV_v20200527(ds,model,mi,xSlc,YDat_YZ,ZDat_YZ,...
                                {zSiV},[],0);
        
        %% assemble coordinates and strain coupling arrays
        % find sector where strain coupling is maximum for this orientation
        % of SiV
        secsgn = cpl.(['SiV_',zSiVstr]).full.LSiVMaxSec(:,mi);
        xSlc = xSlc*secsgn(1);
        
        % identify sectors to plot
        if xSlc < 0 && abs(P.mevenx)
            xSec = 'm';
        else
            xSec = 'p';
        end

        allSec = {'pp'};
        if abs(P.mevenz) && P.meveny == 0
            allSec = [allSec,'pm'];
        elseif abs(P.meveny) && P.mevenz == 0
            allSec = [allSec,'mp'];
        elseif abs(P.meveny) && abs(P.mevenz)
            allSec = [allSec,'mp','pm','mm'];
        end

        % assemble coordinates and reshape arrays in relevant sectors
        for j = 1:length(allSec)
            secLbl = [xSec,allSec{j}];
            
            for k = 1:6
                LDat0 = dsSiV.(['z',zSiVstr]).(secLbl).LSiV;
                LDat = reshape(LDat0(mi,:),[length(YDat_YZ),length(ZDat_YZ)]);
                LDat = LDat'; % transpose array so geom Y/Z-axis --> plot X/Y-axis
                % flip arrays for plotting
                
                if strcmp(secLbl(2),'m')    % -y axis
                    LDat = fliplr(LDat);
                    YPlt = -1*fliplr(YDat_YZ);
                else
                    YPlt = YDat_YZ;
                end
                if strcmp(secLbl(3),'m')    % -z axis
                    LDat = flipud(LDat);
                    ZPlt = -1*fliplr(ZDat_YZ);
                else
                    ZPlt = ZDat_YZ;
                end
                ePlotYZ.(allSec{j}).YDat = YPlt;
                ePlotYZ.(allSec{j}).ZDat = ZPlt;
                ePlotYZ.(allSec{j}).LDat = LDat;
            end
        end

        % if plot axis contain zero coordinate, eliminate row/col corresponding
        % to that coordinate
        yid = 0;
        if ~isempty(find(YDat_YZ==0,1))
            yid = 1;
        end
        zid = 0;
        if ~isempty(find(ZDat_YZ==0,1))
            zid = 1;
        end
        
        y_pp = ePlotYZ.pp.YDat;
        z_pp = ePlotYZ.pp.ZDat;
        L_pp = ePlotYZ.pp.LDat;
        ePlotYZ.full = ePlotYZ.pp;  % if no x & y symmetries


        % combine strain coupling data from all sectors
        if abs(P.meveny)
            L_mp = ePlotYZ.mp.LDat;
            y_mp = ePlotYZ.mp.YDat;
            ePlotYZ.full.LDat = [L_mp(:,1:end-yid),L_pp];
            ePlotYZ.full.YDat = [y_mp(1:end-yid),y_pp];
        end
        if abs(P.mevenz)
            L_pm = ePlotYZ.pm.LDat;
            z_pm = ePlotYZ.pm.ZDat;
            ePlotYZ.full.LDat = [L_pm(1:end-zid,:);L_pp];
            ePlotYZ.full.ZDat = [z_pm(1:end-zid),z_pp];
        end
        if abs(P.meveny) && abs(P.mevenz)
            L_mm = ePlotYZ.mm.LDat;
            ePlotYZ.full.LDat = [L_mm(1:end-zid,1:end-yid),L_pm(1:end-zid,:);...
                                 L_mp(:,1:end-yid),      L_pp];
        end
        
        YPlotYZ = ePlotYZ.full.YDat;
        ZPlotYZ = ePlotYZ.full.ZDat;
        LPlotYZ = ePlotYZ.full.LDat;

        % get max strain coupling
        [LmaxYZ, ImaxYZ] = max(LPlotYZ(:));
        [ImaxZ, ImaxY] = ind2sub(size(LPlotYZ),ImaxYZ);
        maxY = YPlotYZ(ImaxY);
        maxZ = ZPlotYZ(ImaxZ);
        ePlotYZ.full.LmaxYZ = LmaxYZ;
        ePlotYZ.full.maxY = maxY;
        ePlotYZ.full.maxZ = maxZ;
        
        % get slice along z for max strain coupling
        YSlice = ySlc*secsgn(2);
        [~,IslcY] = min(abs(YPlotYZ-YSlice));
        LZSlice = LPlotYZ(:,IslcY);

        % get slice along Y for max strain coupling
        ZSlice = zSlc*secsgn(3);
        [~,IslcZ] = min(abs(ZPlotYZ-ZSlice));
        LYSlice = LPlotYZ(IslcZ,:);
        
        % get max coupling at specified coordinates in plane
        LPmax = LPlotYZ(IslcZ,IslcY);
        
        %% find strain tensor at specified coordinates
        secStr = 'ppp';
        if secsgn(1) == -1
            secStr(1) = 'm';
        end
        if secsgn(2) == -1
            secStr(2) = 'm';
        end
        if secsgn(3) == -1
            secStr(3) = 'm';
        end
        
        [xm,ym,zm] = meshgrid(xSlc,YDat_YZ,ZDat_YZ); % array of size length(yCoords) x length(xCoords(xi)) x length(zCoords)
%         coord = transpose([xm(:),ym(:),zm(:)]);
        yIdx = find(abs(ym(:)-ySlc)<1e-9);
        zIdx = find(abs(zm(:)-zSlc)<1e-9);
        cIdx = intersect(yIdx,zIdx);
        eSiVcoord = dsSiV.(['z',zSiVstr]).(secStr).eSiV(:,:,cIdx,mi);
        eSiVcoord = squeeze(eSiVcoord);
        
        % check
%         dg = P.dg;    
%         fg = P.fg;
%         e_Egx = dg*(eSiVcoord(1,1)-eSiVcoord(2,2)) + fg*eSiVcoord(2,3);
%         e_Egy = -2*dg*eSiVcoord(1,2) + fg*eSiVcoord(1,3);
%         ds_do = sqrt(2)*mfem.freqs/46e9;
%         LSiVm = real(ds_do(mi).*sqrt(abs(e_Egx).^2+abs(e_Egy).^2)...
%                         .*cpl.xzpf(mi)./cpl.maxDisp(mi));
%         disp(LSiVm)
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
        surf(YPlotYZ,ZPlotYZ,LPlotYZ);
        grid('off'); shading('interp'); 
        daspect([1 1 1]);
        cmax = cpl.(['SiV_',zSiVstr]).full.LSiVMax(mi);
        caxis([0 0.5*cmax])
        colormap(colortable('Rainbow'))
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

        % plot location of max strain
        plot([-max(w)/2 max(w)/2],ZSlice*[1 1],'r--','linewidth',1)
        plot(YSlice*[1 1],[-P.th/2 P.th/2],'r--','linewidth',1)

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
        plot(bx,YPlotYZ,LYSlice,'k','linewidth',1.5)
        xlim(bx,[-max(w)/2 max(w)/2])
        ylim(bx,[0 max(LYSlice)])

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
        plot(cx,LZSlice,ZPlotYZ,'k','linewidth',1.5)
%         disp(LZSlice)
        ylim(cx,[-P.th/2 P.th/2])
        if ~isnan(max(LZSlice))
            xlim(cx,[0 max(LZSlice)])
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
        
%         % get strain tensor at 
%         [~,~,dsSiV] = rotateStrainIntoSiV(ds,model,mi,xSlc,ySlc,zSlc,...
%                                 {zSiV});
%         secsgn = cpl.(['SiV_',zSiVstr]).full.LSiVMaxSec(:,mi);
%         xSlc = xSlc*secsgn(1);
%         
%         % identify sectors to plot
%         if xSlc < 0 && abs(P.mevenx)
%             xSec = 'm';
%         else
%             xSec = 'p';
%         end
% 
%         allSec = {'pp'};
%         if abs(P.mevenz) && P.meveny == 0
%             allSec = [allSec,'pm'];
%         elseif abs(P.meveny) && P.mevenz == 0
%             allSec = [allSec,'mp'];
%         elseif abs(P.meveny) && abs(P.mevenz)
%             allSec = [allSec,'mp','pm','mm'];
%         end
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
        
        
        
        
        
        
        
        
        
        
        
        
    end
    
    
    
    
    
    mIdx = mIdx + 1;
end























end