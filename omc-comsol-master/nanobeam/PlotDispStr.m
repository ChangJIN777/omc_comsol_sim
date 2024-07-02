% Script to plot displacement and strain profile of mechanical modes
% Takes COMSOL model created from runNanobeamOptMech


function PlotDispStr(model,ds,mModes,datLoc,fullBeamPlot)
close all;

P = ds.P;
mfem = ds.mfem;
if isfield(ds,'cpl')
    cpl = ds.cpl;
end

phase = 0;  %in degrees
model.result.dataset('mdset').set('phase',phase);

freqs = mfem.freqs;
xzpf = mfem.xzpf;
meff = mfem.meff;

xySymFac = 2^(abs(P.mevenx)+abs(P.meveny));
zSym = abs(P.mevenz*strcmp(P.xsect,'rect'));
zSymFac = 2^zSym;

%% Filenames
if ~isfield(P,'fileBase')
    P = CreateFileBase(P);
    ds.P = P;
    fileBase = P.fileBase;
else
    fileBase = P.fileBase;
end

%% plot full beam
if fullBeamPlot
    dsetTags = mphmodel(model.result.dataset);
    
    if ~isfield(dsetTags,'mdset_sec')
        mdset_sec = model.result.dataset.create('mdset_sec', 'Sector3D');
    else
        mdset_sec = model.result.dataset('mdset_sec');
    end
    % in XY-plane: create full structure with both rotation and reflection
    % by creating 4 sectors
    mdset_sec.label('M Sol Full Beam XY');
    mdset_sec.set('data', 'mdset');
    mdset_sec.set('method', 'twopoint');
    mdset_sec.setIndex('genpoints', '0', 0, 0); % set axis data
    mdset_sec.setIndex('genpoints', '0', 0, 1);
    mdset_sec.setIndex('genpoints', '0', 0, 2);
    mdset_sec.setIndex('genpoints', '0', 1, 0);
    mdset_sec.setIndex('genpoints', '0', 1, 1);
    mdset_sec.setIndex('genpoints', '1', 1, 2);
    mdset_sec.set('sectors', xySymFac);
    mdset_sec.set('trans', 'rotrefl');          % transformation: rotate and reflect
    mdset_sec.set('reflaxis', {'1' '0' '0'});   % reflection axis
    if (P.mevenx==-1)
        mdset_sec.set('rotinv', 'on');  % invert phase for x antisymmetry
    end
    if (P.meveny==-1)
        mdset_sec.set('reflinv', 'on');  % invert phase for y antisymmetry
    end
    mdsetPlot = 'mdset_sec';
    
    if zSym
        % reflect about XY-plane to create 2 sectors of full XY structure
        if ~isfield(dsetTags,'mdset_secZ')
            mdset_sec = model.result.dataset.create('mdset_secZ', 'Sector3D');
        else
            mdset_sec = model.result.dataset('mdset_secZ');
        end
        mdset_sec.label('M Sol Full Beam XYZ');
        mdset_sec.set('data', 'mdset_sec');
        mdset_sec.set('method', 'twopoint');
        mdset_sec.setIndex('genpoints', '0', 0, 0); % set axis data
        mdset_sec.setIndex('genpoints', '0', 0, 1);
        mdset_sec.setIndex('genpoints', '0', 0, 2);
        mdset_sec.setIndex('genpoints', '1', 1, 0);
        mdset_sec.setIndex('genpoints', '0', 1, 1);
        mdset_sec.setIndex('genpoints', '0', 1, 2);
        mdset_sec.set('sectors', zSymFac);
        mdset_sec.set('trans', 'rotrefl');          % transformation: rotate and reflect
        mdset_sec.set('reflaxis', {'0' '1' '0'});   % reflection axis
        if (P.mevenz==-1)
            mdset_sec.set('reflinv', 'on'); % invert phase for z antisymmetry
        end
        mdsetPlot = 'mdset_secZ';
        
    end
    
else
    mdsetPlot = 'mdset';
end


%% Get max displacements
dispMaxAll = zeros(1,P.mneigs);
for mi = mModes
    dispMaxAll(mi) = mphmax(model,'abs(solid.disp)','volume','dataset','mdset','solnum',mi);
end

%% Save mechanical displacement profile

% hide air cylinder if present
geomnames = fieldnames(mphmodel(model.geom));
geomname = geomnames{1};
NDoms = model.geom(geomname).getNDomains;
if NDoms == 2
    if isfield(mfem,'air_domind')
        airdom = mfem.air_domind;
    else
        airdom = setdiff([1,2],mfem.dia_domind);
    end
    
    hideTags = mphmodel(model.view('view1').hideEntities);
    if ~isfield(hideTags,'hide1')
        model.view('view1').hideEntities.create('hide1');
        model.view('view1').hideEntities('hide1').set(airdom);
    end
end

resTags = mphmodel(model.result);

if ~isfield(resTags,'dispplot')
    % create 3D deformation plot group for displacement field
    dispplot = model.result.create('dispplot', 'PlotGroup3D');
    dispplot_vol = dispplot.create('dispplot_vol', 'Volume');
    dispplot_def = dispplot_vol.create('dispplot_def', 'Deform');
else
    dispplot = model.result('dispplot');
    dispplot_vol = dispplot.feature('dispplot_vol');
    dispplot_def = dispplot_vol.feature('dispplot_def');
end
dispplot.label('Displacement Deformation Plots');
dispplot.set('data', mdsetPlot);
dispplot.set('edges', 'off');
dispplot.set('titletype', 'none');
dispplot_vol.set('data', 'parent');
dispplot_vol.set('rangecoloractive', 'on');
dispplot_def.set('expr', {'u' 'v' 'w'});
dispplot_def.set('scaleactive','on');
% dispplot.run;

for mi = mModes
    % Filename and title settings
    mfilestr = ['wM_',num2str(mfem.freqs(mi)*1e-9,'%.2f'),'GHz'];
    
    mtitlestr = ['\omega_M = ',num2str(mfem.freqs(mi)*1e-9,'%.2f'),'GHz, ',...
                 'm_{eff} = ',num2str(mfem.meff(mi)*1e18,'%.1f'),'fg, ',...
                 'x_{zpf} = ',num2str(mfem.xzpf(mi)*1e15,'%.2f'),'fm'];
    
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
    
    %% disp plot
    dispMax = dispMaxAll(mi);
    dispExpr = ['solid.disp/',num2str(dispMax),'*',num2str(xzpf(mi))];
%     dispExpr = ['solid.disp/',num2str(dispMax)];
    
    dispplot.set('solnum',mi);
    dispplot_vol.set('expr',dispExpr);
    dispplot_vol.set('rangecolormax',num2str(xzpf(mi)));
%     dispplot_vol.set('rangecolormax',num2str(1));
    dispplot_vol.set('rangecolormin',num2str(0));
    dispplot_def.set('scale', [num2str(dispMax),'*1e-9']);
    dispplot.run;
    
    figure;
    subplot(2,1,1);
    try % Plot XY displacement profile
        mphplot(model,'dispplot');
    catch err
        display(err.message)
        for lsi=1:length(err.stack)
            err.stack(lsi)
        end
    end
    axis off
    daspect([1 1 1])
    view(2)
    zoom(1)
    
    %title
    plotTitle = {['Displacement, ',mtitlestr];ctitlestr};
    title(plotTitle,'fontname','arial','fontsize',10)

    subplot(2,1,2);
    try % Plot YZ displacement profile
        mphplot(model,'dispplot');
    catch err
        display(err.message)
        for lsi=1:length(err.stack)
            err.stack(lsi)
        end
    end
    axis off
    daspect([1 1 1])
    view(-90,0)
%     zoom(20)
    camzoom(8)  % for newer MATLAB versions
    colorbar('eastoutside');
    caxis([0 xzpf(mi)]);
%     caxis([0 1]);
    colormap(colortable('Rainbow'));
    
    pathFig = [datLoc,fileBase,'_',mfilestr,'_Disp.png'];
    saveas(gcf,pathFig);
    close;
end

% close all;

%% Save mechanical strain profile

strExpr = '(solid.eXX+solid.eYY+solid.eZZ)';
% strExpr = 'solid.eYY';

if ~isfield(P,'celltype') 
    if P.holeatctr
        ctrHoleRX = P.xc + ds.P.hx_hole(1)/2;   % central hole right-x-coordinate
        adjHoleLX = P.xc + ds.P.a_hole(2) - ds.P.hx_hole(2)/2; % adjacent hole left-x-coordinate
        xSlc = 0.5*(ctrHoleRX+adjHoleLX);
    else
        xSlc = 0;
    end
else
    xSlc = 0;
end

resTags = mphmodel(model.result);
if ~isfield(resTags,'strplot')
    % create 3D slice plot group for strain field
    strplot = model.result.create('strplot', 'PlotGroup3D');
    strplot_slcXY = strplot.create('strplot_slcXY', 'Slice');
    strplot_slcYZ = strplot.create('strplot_slcYZ', 'Slice');
else
    strplot = model.result('strplot');
    strplot_slcXY = strplot.feature('strplot_slcXY');
    strplot_slcYZ = strplot.feature('strplot_slcYZ');
end
strplot.label('Strain Profile');
strplot.set('data', mdsetPlot);
strplot.set('titletype', 'none');
strplot_slcXY.set('planetype', 'quick').set('data', 'parent');
strplot_slcXY.set('quickzmethod', 'coord').set('quickz',P.th/2);
strplot_slcXY.set('rangecoloractive', 'on');
strplot_slcXY.set('colortable', 'WaveLight');
strplot_slcYZ.set('planetype', 'quick').set('data', 'parent');
strplot_slcYZ.set('quickxmethod', 'coord').set('quickx',xSlc);%P.a*(1-P.maxdef)/2+P.xc);
strplot_slcYZ.set('rangecoloractive', 'on');
strplot_slcYZ.set('colortable', 'WaveLight');
%     strplot.run;

for mi = mModes
    % Filename and title settings
    mfilestr = ['wM_',num2str(mfem.freqs(mi)*1e-9,'%.2f'),'GHz'];
    
    mtitlestr = ['\omega_M = ',num2str(mfem.freqs(mi)*1e-9,'%.2f'),'GHz, ',...
                 'm_{eff} = ',num2str(mfem.meff(mi)*1e18,'%.1f'),'fg, ',...
                 'x_{zpf} = ',num2str(mfem.xzpf(mi)*1e15,'%.2f'),'fm'];
    
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
    
    %% str plot
    dispMax = dispMaxAll(mi);
    strMax = mphmax(model,strExpr,'volume','dataset','mdset','solnum',mi);
    strxpfExpr = [strExpr,'/',num2str(dispMax),'*',num2str(xzpf(mi))];

    strplot.set('solnum', mi);
    strplot_slcXY.set('expr',strxpfExpr);
    strplot_slcXY.set('rangecolormax',num2str(strMax/dispMax*xzpf(mi)));
    strplot_slcXY.set('rangecolormin',num2str(-strMax/dispMax*xzpf(mi)));
    strplot_slcYZ.set('expr',strxpfExpr);
    strplot_slcYZ.set('rangecolormax',num2str(strMax/dispMax*xzpf(mi)));
    strplot_slcYZ.set('rangecolormin',num2str(-strMax/dispMax*xzpf(mi)));
    strplot.run;

    % Plot XY displacement profile
    figure;
    subplot(2,1,1);
    try 
        mphplot(model,'strplot');
    catch err
        display(err.message)
        for lsi=1:length(err.stack)
            err.stack(lsi)
        end
    end
    box on
    axis tight
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');
    set(gca,'ZTickLabel','');
    daspect([1 1 1])
    view(2)
    zoom(1)

    plotTitle = {['Strain (XX+YY+ZZ), ',mtitlestr];ctitlestr};
    title(plotTitle,'fontname','arial','fontsize',10)


    % Plot YZ displacement profile
    subplot(2,1,2);
    try 
        mphplot(model,'strplot');
    catch err
        display(err.message)
        for lsi=1:length(err.stack)
            err.stack(lsi)
        end
    end
    box on
    axis tight
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');
    set(gca,'ZTickLabel','');
    daspect([1 1 1])
    view(-90,0)
%     zoom(20)
    camzoom(8)  % for newer MATLAB versions
    colorbar('eastoutside');
    caxis([-strMax/dispMax*xzpf(mi) strMax/dispMax*xzpf(mi)]);
    colormap(colortable('WaveLight'));

    pathFig = [datLoc,fileBase,'_',mfilestr,'_Strain.png'];
    saveas(gcf,pathFig);
    close;
end
    
close all; 

if P.solveOpt
    hideTags = mphmodel(model.view('view1').hideEntities);
    if isfield(hideTags,'hide1')
        model.view('view1').hideEntities.remove('hide1');
    end
end

end