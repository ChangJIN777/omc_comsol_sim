% function to generate SiV strain coupling plots
% inputs:
% ds: data structure from simulation
% model: COMSOL model with mechanical simulations
% mModes: solution no. of mechanical eigenmodes to plot for
% xSlc, ySlc,zSlc: x-, y-, z-coordinate to plot slice at
% zSiV = [1 1 1], [-1 1 1], [-1 -1 1], [1 -1 1]:
%   z-orientation of SiV


function PlotStrCplSiV(ds,model,dsSiV,mModes,xSlc,ySlc,zSlc,zSiVAll,datLoc)
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
    
    %% get coordinates for strain interpolation
    interpRes = 10e-9;  % resolution interpolation
    
    % determine full extents of beam - especially important for triangular
    % cross-section, where width of beam is function of depth under surface
    if strcmp(P.xsect,'tri')
        wPlot = (zSlc(mIdx)+P.th/2)/P.th*P.w; 
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
    XDat_XY = linspace(0,beamLen,round(beamLen/interpRes)+1);
    YDat_XY = linspace(0,wPlot/2,round(wPlot/2/interpRes)+1);
    YDat_YZ = linspace(0,max(P.w_hole)/2,round(max(P.w_hole)/2/interpRes)+1);
    ZDat_YZ = linspace(ZMin,P.th/2,round((P.th/2-ZMin)/interpRes)+1);
    
    % insert custom x/y/z coordinates
    if ~isempty(xSlc(mIdx)) && isempty(find(abs(XDat_XY-xSlc(mIdx))<1e-9,1))
        XDat_XY = sort([XDat_XY,xSlc(mIdx)]);
    end
    if ~isempty(ySlc(mIdx)) && isempty(find(abs(YDat_XY-ySlc(mIdx))<1e-9,1))
        YDat_XY = sort([YDat_XY,ySlc(mIdx)]);
        YDat_YZ = sort([YDat_YZ,ySlc(mIdx)]);
    end
    if ~isempty(zSlc(mIdx)) && isempty(find(abs(ZDat_YZ-zSlc(mIdx))<1e-9,1))
        ZDat_YZ = sort([ZDat_YZ,zSlc(mIdx)]);
    end
    
    %% XY: Compile strain coupling profiles for plotting
    % get strain coupling
    [~,~,dsSiV] = CalcStrCplSiV(ds,model,mi,XDat_XY,YDat_XY,zSlc(mIdx),...
                                zSiVAll,runStats);
    
    % identify sectors to plot
    if zSlc(mIdx) < 0 && abs(P.mevenz)
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
    
    % assemble coordinates and strain coupling arrays
    for j = 1:length(allSec)
        secLbl = [allSec{j},zSec];
        
        % assemble coordinate arrays for full beam
        if strcmp(secLbl(1),'m')    % -x axis
            eDat1 = fliplr(eDat1);
            XPlt = -1*fliplr(XDat_XY);
        else
            XPlt = XDat_XY;
        end
        if strcmp(secLbl(2),'m')    % -y axis
            eDat1 = flipud(eDat1);
            YPlt = -1*fliplr(YDat_XY);
        else
            YPlt = YDat_XY;
        end
        
        ePlotXY.(allSec{j}).XDat = XPlt;
        ePlotXY.(allSec{j}).YDat = YPlt;
    end
    
    % if plot axis contain zero coordinate, eliminate row/col corresponding
    % to that coordinate
    xid = 1;
    if ~isempty(find(XDat_XY==0,1))
        xid = 2;
    end
    yid = 1;
    if ~isempty(find(YDat_XY==0,1))
        yid = 2;
    end
    
    % reshape strain coupling data for all sectors
    if abs(P.mevenx)
        eij_mp = ePlotXY.mp.(eijLbl);
        x_mp = ePlotXY.mp.XDat;
        ePlotXY.full.(eijLbl) = [eij_mp(:,xid:end),eij_pp];
        ePlotXY.full.XDat = [x_mp(xid:end),x_pp];
    end
    if abs(P.meveny)
        eij_pm = ePlotXY.pm.(eijLbl);
        y_pm = ePlotXY.pm.YDat;
        ePlotXY.full.(eijLbl) = [eij_pm(yid:end,:);eij_pp];
        ePlotXY.full.YDat = [y_pm(yid:end),y_pp];
    end
    if abs(P.meveny) && abs(P.mevenx)
        eij_mm = ePlotXY.mm.(eijLbl);
        ePlotXY.full.(eijLbl) = [eij_mm(yid:end,xid:end),eij_pm(yid:end,:);...
                                 eij_mp(:,xid:end),      eij_pp];
    end
    
    
    
    
    
    
end























end