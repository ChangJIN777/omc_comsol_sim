% Script to plot displacement and strain profile of mechanical modes
% Takes COMSOL model created from runNanobeamOptMech


function PlotEy(model,ds,oModes,datLoc,fullBeamPlot)
close all;

P = ds.P;
ofem = ds.ofem;
c = 299792458;

xySymFac = 2^(abs(P.oevenx)+abs(P.oeveny));
zSym = abs(P.oevenz*strcmp(P.xsect,'rect'));
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
    if ~isfield(dsetTags,'odset_sec')
        odset_sec = model.result.dataset.create('odset_sec', 'Sector3D');
    else
        odset_sec = model.result.dataset('odset_sec');
    end
    % in XY-plane: create full structure with both rotation and reflection
    % by creating 4 sectors
    odset_sec.label('O Sol Full Beam XY');
    odset_sec.set('data', 'odset');
    odset_sec.set('method', 'twopoint');
    odset_sec.setIndex('genpoints', '0', 0, 0); % set axis data
    odset_sec.setIndex('genpoints', '0', 0, 1);
    odset_sec.setIndex('genpoints', '0', 0, 2);
    odset_sec.setIndex('genpoints', '0', 1, 0);
    odset_sec.setIndex('genpoints', '0', 1, 1);
    odset_sec.setIndex('genpoints', '1', 1, 2);
    odset_sec.set('sectors', xySymFac);
    odset_sec.set('trans', 'rotrefl');          % transformation: rotate and reflect
    odset_sec.set('reflaxis', {'1' '0' '0'});   % reflection axis
%     if (P.oevenx==-1) && (P.oeveny==-1)
%         odset_sec.set('rotinv', 'on');  % odd symmetry about x-plane
%     end
    if (P.oevenx==-1)
        odset_sec.set('rotinv', 'on');  % invert phase for x antisymmetry
    else
        odset_sec.set('rotinv', 'off');
    end
    if (P.oeveny==-1)
        odset_sec.set('reflinv', 'off');  % invert phase for y antisymmetry
    else
%         odset_sec.set('reflinv', 'off');
    end
    odsetPlot = 'odset_sec';
    
    if zSym
        % reflect about XY-plane to create 2 sectors of full XY structure
        if ~isfield(dsetTags,'odset_secZ')
            odset_sec = model.result.dataset.create('odset_secZ', 'Sector3D');
        else
            odset_sec = model.result.dataset('odset_secZ');
        end
        odset_sec.label('O Sol Full Beam XYZ');
        odset_sec.set('data', 'odset_sec');
        odset_sec.set('method', 'twopoint');
        odset_sec.setIndex('genpoints', '0', 0, 0); % set axis data
        odset_sec.setIndex('genpoints', '0', 0, 1);
        odset_sec.setIndex('genpoints', '0', 0, 2);
        odset_sec.setIndex('genpoints', '1', 1, 0);
        odset_sec.setIndex('genpoints', '0', 1, 1);
        odset_sec.setIndex('genpoints', '0', 1, 2);
        odset_sec.set('sectors', zSymFac);
        odset_sec.set('trans', 'rotrefl');          % transformation: rotate and reflect
        odset_sec.set('reflaxis', {'0' '1' '0'});   % reflection axis
%         if (P.oevenz==-1)
%             odset_sec.set('rotinv', 'on');  % odd symmetry about x-plane
%             odset_sec.set('reflinv', 'on');  
%         end
        if (P.oevenz==-1)
            odset_sec.set('reflinv', 'on'); % invert phase for z antisymmetry
        end
        odsetPlot = 'odset_secZ';
        
    end
    
else
    odsetPlot = 'odset';
end

% draw nanobeam and cylinder
airdom = ofem.air_domind;
hideTags = mphmodel(model.view('view1').hideEntities);
if isfield(hideTags,'hide1')
    model.view('view1').hideEntities.remove('hide1');
end

% create 3D plot group for Ey on top plane and cross section of beam 
resTags = mphmodel(model.result);
if ~isfield(resTags,'eytopplot')
    eytopplot = model.result.create('eytopplot', 'PlotGroup3D');
    eytopplot_slc = eytopplot.create('eytopplot_slc', 'Slice');
else
    eytopplot = model.result('eytopplot');
    eytopplot_slc = eytopplot.feature('eytopplot_slc');
end
eytopplot.label('Electric Field, y-component, top view');
eytopplot.set('data', odsetPlot);
eytopplot.set('titletype', 'none');
eytopplot_slc.set('planetype', 'quick');
eytopplot_slc.set('quickzmethod', 'coord').set('quickz', P.th/2);
eytopplot_slc.set('rangecoloractive', 'on');
eytopplot_slc.set('rangecolormax',1).set('rangecolormin',-1);
eytopplot_slc.set('colortable', 'WaveLight');

if strcmp(P.celltype,'snowflake') || strcmp(P.celltype,'boomerang_v2')
    dielCtrX = P.xc;
else
    ctrHoleRX = P.xc + ds.P.hx_hole(1)/2;
    adjHoleLX = P.xc + ds.P.a_hole(2) - ds.P.hx_hole(2)/2;
    dielCtrX = (ctrHoleRX+adjHoleLX)/2;
end

if ~isfield(resTags,'eyxsecplot')
    eyxsecplot = model.result.create('eyxsecplot', 'PlotGroup3D');
    eyxsecplot_slc = eyxsecplot.create('eyxsecplot_slc', 'Slice');
else
    eyxsecplot = model.result('eyxsecplot');
    eyxsecplot_slc = eyxsecplot.feature('eyxsecplot_slc');
end
eyxsecplot.label('Electric Field, y-component, cross section');
eyxsecplot.set('data', odsetPlot);
eyxsecplot.set('titletype', 'none');
eyxsecplot_slc.set('planetype', 'quick');
eyxsecplot_slc.set('quickxmethod', 'coord').set('quickx', dielCtrX);
eyxsecplot_slc.set('rangecoloractive', 'on');
eyxsecplot_slc.set('rangecolormax',1).set('rangecolormin',-1);
eyxsecplot_slc.set('colortable', 'WaveLight');



for oi = oModes
    
    absEyMax = mphmax(model,'abs(emw.Ey)','volume','dataset','odset','solnum',oi);

    eytopplot.set('solnum', oi);
    eytopplot_slc.set('expr', ['emw.Ey/',num2str(absEyMax)]);
    eytopplot.run;
    
    eyxsecplot.set('solnum', oi);
    eyxsecplot_slc.set('expr', ['emw.Ey/',num2str(absEyMax)]);
    eyxsecplot.run;

    figure;
    subplot(2,1,1)
    try % Plot Ey in triangular cross-section
        mphplot(model,'eytopplot');
    catch err
        display(err.message)
        for lsi=1:length(err.stack)
            err.stack(lsi)
        end
    end
    
    view(2)
    zoom(1);
%     colorbar('southoutside');
%     caxis([-1 1]);
    colormap(colortable('WaveLight'));
    axis off
    ylim([-2e-6 2e-6])
    daspect([1 1 1])
    
    % plot title
    otitlestr = ['\lambda_O = ',num2str(ofem.lambdaAll(oi)*1e9,'%.0f'),'nm, ',...
                 'Q = ',num2str(ofem.QAll(oi)*1e-5,'%.2f'),'x 10^5'];
    plotTitle = ['E_y, ',otitlestr];
    title(plotTitle,'fontname','arial','fontsize',10)
    
    subplot(2,1,2)
    try % Plot Ey in triangular cross-section
        mphplot(model,'eyxsecplot');
    catch err
        display(err.message)
        for lsi=1:length(err.stack)
            err.stack(lsi)
        end
    end
    
%     zoom(1);
    colorbar('eastoutside');
    caxis([-1 1]);
    colormap(colortable('WaveLight'));
    axis off
    zlim(2*P.th*[-1 1])
    ylim(2*P.w*[-1 1])
    daspect([100 1 1])
    view(-90,0)
    
    %% Saving
    % filename
    ofilestr = ['wl_',num2str(ofem.lambdaAll(oi)*1e9,'%.0f'),'nm_',...
                'Q_',num2str(ofem.QAll(oi)*1e-5,'%.2f'),'E5'];
    pathFig = [datLoc,fileBase,'_',ofilestr,'_Ey.png'];
    saveas(gcf,pathFig);
end

end