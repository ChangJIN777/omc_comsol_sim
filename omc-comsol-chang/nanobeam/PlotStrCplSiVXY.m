% function to generate SiV strain coupling plots
% inputs:
% ds: data structure from simulation
% model: COMSOL model with mechanical simulations
% mModes: solution no. of mechanical eigenmodes to plot for
% xSlc, ySlc,zSlc: x-, y-, z-coordinate to plot slice at
% zSiV = [1 1 1], [-1 1 1], [-1 -1 1], [1 -1 1]:
%   z-orientation of SiV


function varargout = PlotStrCplSiVXY(ds,model,mModes,pSlcAll,zSiVAll,datLoc)
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
            XDat_XY = sort([XDat_XY,abs(xInt)]);
        end
        if ~isempty(yInt) && isempty(find(abs(YDat_XY-yInt)<1e-9,1))
            YDat_XY = sort([YDat_XY,abs(yInt)]);
            YDat_YZ = sort([YDat_YZ,abs(yInt)]);
        end
        if ~isempty(zInt) && isempty(find(abs(ZDat_YZ-zInt)<1e-9,1))
            ZDat_YZ = sort([ZDat_YZ,zInt]);
        end
        
        [~,~,dsSiV] = CalcStrCplSiV(ds,model,mi,XDat_XY,YDat_XY,zInt,...
                                {zSiV},[],0);
        
        %% assemble coordinates and strain coupling arrays
        % find sector of max coupling
        maxSec = dsSiV.LSiVMaxSec(1,1,1,1);
        maxSecSgn = secSgnAll(:,maxSec);
        
        % identify sectors to plot depending on z-coordinate
        secsPlt = find((secSgnAll(3,secIdxs))==maxSecSgn(3));
        
        % extract strain coupling and coordinate arrays
        % and apply transforms depending on sector
        LDatAll = zeros(length(YDat_XY),length(XDat_XY),length(secsPlt));
        XDatAll = zeros(length(XDat_XY),length(secsPlt));
        YDatAll = zeros(length(YDat_XY),length(secsPlt));
        
       for sIdx = 1:length(secsPlt)
            si = secsPlt(sIdx);
            
            % extract and reshape arrays
            LDat0 = dsSiV.LSiV(1,:,1,1,si);
            LDatS = reshape(LDat0,[length(YDat_XY),length(XDat_XY)]);
            
            % flip and/or change signs of arrays
            xgt0 = (secSgnAll(1,si) == 1);
            ygt0 = (secSgnAll(2,si) == 1);
            
            XDatS = xgt0*XDat_XY + ~xgt0*-1*fliplr(XDat_XY);
            YDatS = ygt0*YDat_XY + ~ygt0*-1*fliplr(YDat_XY);
            
            if ~xgt0
                LDatS = fliplr(LDatS);
            end
            if ~ygt0
                LDatS = flipud(LDatS);
            end
            
            % assemble arrays
            XDatAll(1:length(XDat_XY),sIdx) = XDatS(1:length(XDat_XY));
            YDatAll(1:length(YDat_XY),sIdx) = YDatS(1:length(YDat_XY));
            LDatAll(1:length(YDat_XY),1:length(XDat_XY),sIdx) = ...
                LDatS(1:length(YDat_XY),1:length(XDat_XY));
        end
        
        % get max strain coupling and position
        [LmaxXY, ImaxXY] = max(LDatAll(:));
        [ImaxY, ImaxX, ImaxS] = ind2sub(size(LDatAll),ImaxXY);
        maxX = XDatAll(ImaxX,ImaxS);
        maxY = YDatAll(ImaxY,ImaxS);
        
        % get slice along X at custom Y for plotting
        % find sectors with custom Y
        ySlcSecSgn = (sign(ySlc))^(abs(P.meveny)) + (sign(ySlc)==0);
        secsYSlc = find((secSgnAll(2,:))==ySlcSecSgn);
        [secsYSlc,secsYSlcIdx] = intersect(secsPlt,secsYSlc);
%         disp(secsYSlc)
%         disp(secsYSlcIdx)
%         secsYSlcIdx = find(secsPlt == secsYSlc);
        [~,IslcY] = min(abs(YDatAll(:,secsYSlcIdx(1))-ySlc));
%         disp(IslcY)
        %% XY surface plot
        figure;
        set(gcf,'position',[100 100 913 500])
        axpos = [0.07 0.55 0.82 0.250];
        ax = axes('position',axpos); hold(ax,'on')
        box on
        axis off

        % strain coupling
        hold on
        for sIdx = 1:length(secsPlt)
            surf(XDatAll(:,sIdx),YDatAll(:,sIdx),LDatAll(:,:,sIdx))
        end
%         surf(XPlotXY,YPlotXY,LPlotXY);
        grid('off'); shading('interp'); 
        daspect([1 1 1]);
        cmax = cpl.SiV.LSiVMaxMode(mi);
        caxis([0 0.5*cmax])        
        cmap = colormap(hot);
        colormap(cmap(1:56,:));
        cb = colorbar('eastoutside');
        cbp = get(cb,'Position');
        set(cb,'Position',[0.9 axpos(2) 0.03 0.25]);

        % geometry overlay
        gXL = [min(XDatAll(:)) max(XDatAll(:))];
%         gXL = [XPlotXY(1) XPlotXY(end)];
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
        plot(ax,gXL,ySlc*[1 1],'r:','linewidth',0.5)
        xlim(ax,[min(gXU) max(gXU)])

        % settings for overlay to appear on top
        view(0,-90); 
        set(gca,'YDir','Reverse')

        % slice plot
        bx = axes('position',axpos-[0 0.3 0 0]); hold(bx,'on')
        box on
        axis on
        hold on
        for sIdx = secsYSlcIdx'
            plot(bx,XDatAll(:,sIdx),LDatAll(IslcY,:,sIdx),'k','linewidth',1.5)
        end
%         plot(bx,XPlotXY,LXSlice,'k','linewidth',1.5)
        xlim(bx,[min(gXU) max(gXU)])
        LXSlice = LDatAll(IslcY,:,secsYSlcIdx);
        ylim(bx,[0 max(LXSlice(:))])
%         ylim(bx,[0 max(LXSlice)])
        for j=1:length(hy)
            plot(bx,[xpos(j),xpos(j)],[0 max(LXSlice(:))],'--','Color',[0.5 0.5 0.5],'linewidth',0.5)
        end
        title(['Cut along X-axis at y = ',num2str(ySlc*1e9,'%.1f'),...
               'nm (red dotted horizontal line)'],'fontname','arial','fontsize',10)
%         xticks(xpos)
        ylabel('\lambda_{SiV}')
        xlabel('x')
        linkaxes([bx,ax],'x')

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
        % plot title
        plotTitle = {['Strain coupling, SiV axis along [',num2str(zSiV),'], ',...
                    mtitlestr];...
                    ctitlestr;...
                    ['Cut along XY-plane at depth = ',num2str((P.th/2-zSlc)*1e9,'%.1f'),'nm: ',...
                     'at (x,y) = (',num2str((maxX-P.xc)*1e9,'%.1f'),...
                     ',',num2str(maxY*1e9,'%.1f'),') nm, ',...
                     '\lambda_{SiV} = ',num2str(LmaxXY*1e-6,'%.2f'),'MHz'];''};
        title(ax,plotTitle,'fontname','arial','fontsize',10)

        % save file
        pathFig = [datLoc,fileBase,'_',mfilestr,'_',...
                   'd_',num2str((P.th/2-zSlc)*1e9,'%.1f'),'nm_',...
                   zSiVstr,'_XY'];
        saveas(gcf,[pathFig,'.png']);
        saveas(gcf,[pathFig,'.fig']);
        close;
        clear dsSiV
        
        dsSiVplt.LSiVSlc(mIdx,oi) = LmaxXY;
        
        
        
        
        
        
        
        
        
    end
    dsSiVplt.LSiVSlcModes(mIdx) = mi;
    
    
    
    

end

if nargout >= 1
    varargout{1} = dsSiVplt;
end






















end