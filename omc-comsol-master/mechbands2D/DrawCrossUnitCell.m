
function [model,P] = DrawCrossUnitCell(model,P)

% extract parameters from P
xsect = P.xsect;
celltype = P.celltype;
a = P.a;
w = P.w;
th = P.th;
cx = P.cx;
cy = P.cy;
tx = P.tx;
ty = P.ty;
filRad = P.filRad;
evenz = P.mbevenz;


% lithographic constraints
if strcmp(celltype,'hollow')
    if cx >= a || cy >= w
        error('cross dimensions larger than unit cell dimensions')
    end
    if tx >= cx || ty >= cy
        error('cross leg dimensions larger than cross dimensions')
    end
    if filRad > 0 
        % calculate maximum permissible fillet radius
        l1 = 0.5*tx;
        l2 = 0.25*(cx-tx);
        l3 = 0.25*(cy-ty);
        l4 = 0.5*ty;
        minFilRad = min([l1,l2,l3,l4]);
        if filRad > minFilRad
            error(['Fillet radius of ',...
                   num2str(filRad*1e9,'%.1f'),'nm too large, ',...
                   'maximum permissible radius is ',...
                   num2str(minFilRad*1e9,'%.3f'),'nm'])
%             filRad = minFilRad;
        end
    end
end

if strcmp(celltype,'solid')
    if cx >= a || cy >= w
        error('block dimensions larger than unit cell dimensions')
    end
    if tx >= cx || ty >= cy
        error('linkage dimensions larger than block dimensions')
    end
    if filRad > 0 
        % calculate maximum permissible fillet radius
        l1 = 0.5*tx;
        l2 = 0.25*(w-cy);
        l3 = 0.25*(cy-ty);
        l4 = 0.25*(a-cx);
        l5 = 0.25*(cx-tx);
        l6 = 0.5*ty;
        minFilRad = min([l1,l2,l3,l4,l5,l6]);
        if filRad > minFilRad
            error(['Fillet radius of ',...
                   num2str(filRad*1e9,'%.1f'),'nm too large, ',...
                   'maximum permissible radius is ',...
                   num2str(minFilRad*1e9,'%.3f'),'nm'])
%             filRad = minFilRad;
        end
    end
end

%% Create rectangular unit cell
ucellgeom = model.geom.create('ucell', 3);
% create workplane for unit cell
ucellWP = ucellgeom.feature.create('ucellWP', 'WorkPlane');
ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
ucellPlane = ucellWP.geom.feature.create('ucellPlane', 'Rectangle');
ucellPlane.set('type', 'solid').set('base', 'center');
ucellPlane.set('pos', [0 0]).set('size', [a w]);

%% Create solid cross pattern (block + thin linkages)
% generate crosses as for hollow pattern, then use block - crosses at four
% corners
if strcmp(celltype,'solid')
    % create workplane for cross, set at (-a/2,w/2)
    horzbarPlane = ucellWP.geom.create('horzbarPlane', 'Rectangle');
    horzbarPlane.set('type', 'solid').set('base', 'center');
    horzbarPlane.set('pos', [-a/2 -w/2]).set('size', [(a-tx) (w-cy)]);
    ucellgeom.runCurrent;
    vertbarPlane = ucellWP.geom.create('vertbarPlane', 'Rectangle');
    vertbarPlane.set('type', 'solid').set('base', 'center');
    vertbarPlane.set('pos', [-a/2 -w/2]).set('size', [(a-cx) (w-ty)]);
    ucellgeom.runCurrent;
    
    % compose cross
    crossComp = ucellWP.geom.feature.create('crossComp', 'Compose');
    crossComp.selection('input').set({'horzbarPlane' 'vertbarPlane'});
    crossComp.set('formula', 'horzbarPlane + vertbarPlane');
    crossComp.set('intbnd', 'off');
    ucellgeom.runCurrent;
    
    % array crosses at four corners of unit cell in xy plane
    crossArr = ucellWP.geom.feature.create('crossArr', 'Array');
    crossArr.selection('input').set('crossComp');
    crossArr.set('displ', [a w]).set('size', [2 2]);
    ucellgeom.runCurrent;
    
    % compose unit cell plane
    ucellComp = ucellWP.geom.feature.create('ucellComp', 'Compose');
    ucellComp.selection('input').set({'crossArr' 'ucellPlane'});
    ucellComp.set('formula', 'ucellPlane - crossArr(1,1) - crossArr(1,2) - crossArr(2,1) - crossArr(2,2)');
    ucellgeom.runCurrent;
    
    % fillet
    if filRad > 0
        ucellFillet = ucellWP.geom.create('ucellFillet', 'Fillet');
        ucellFillet.set('radius', filRad);
        
        % selecting vertices to fillet - vertices that do not lie on
        % periodic boundary edges
        % 1. get edges along periodic boundaries
        g = model.geom('ucell').obj('ucellWP');
        edgXL = edgeindex(g,[-a/2,   0,-th/2],[1,0,0]); % X left edge
        edgXR = edgeindex(g,[ a/2,   0,-th/2],[1,0,0]); % X right edge
        edgYT = edgeindex(g,[   0, w/2,-th/2],[0,1,0]); % Y top edge
        edgYB = edgeindex(g,[   0,-w/2,-th/2],[0,1,0]); % Y bottom edge
        edgPBCs = [edgXL,edgXR,edgYT,edgYB];
        
        % 2. find vertices connected by edges 
        vtx_adjto_edg = g.getAdj(1,0);
        vtxPBCs = [];
        for edg = edgPBCs
            vtxPBCs = [vtxPBCs;vtx_adjto_edg{edg+1}]; % add 1 to edge index as first element in vtx_adjto_edg is {}
        end
        vtxPBCs = unique(vtxPBCs);
        
        % 3. exclude these vertices from fillet
        vtxFil = setdiff(1:g.getNVertices(),vtxPBCs);
        ucellFillet.selection('point').set('ucellComp', vtxFil);
        ucellgeom.runCurrent;
    end
end

%% Create hollow cross pattern (block - cross)
if strcmp(celltype,'hollow')
    % create workplane for cross, set at (-a/2,w/2)
    horzbarPlane = ucellWP.geom.create('horzbarPlane', 'Rectangle');
    horzbarPlane.set('type', 'solid').set('base', 'center');
    horzbarPlane.set('pos', [0 0]).set('size', [cx ty]);
    ucellgeom.runCurrent;
    vertbarPlane = ucellWP.geom.create('vertbarPlane', 'Rectangle');
    vertbarPlane.set('type', 'solid').set('base', 'center');
    vertbarPlane.set('pos', [0 0]).set('size', [tx cy]);
    ucellgeom.runCurrent;
    
    % compose cross
    crossComp = ucellWP.geom.feature.create('crossComp', 'Compose');
    crossComp.selection('input').set({'horzbarPlane' 'vertbarPlane'});
    crossComp.set('formula', 'horzbarPlane + vertbarPlane');
    crossComp.set('intbnd', 'off');
    ucellgeom.runCurrent;
    
    % compose unit cell plane
    ucellComp = ucellWP.geom.feature.create('ucellComp', 'Compose');
    ucellComp.selection('input').set({'crossComp' 'ucellPlane'});
    ucellComp.set('formula', 'ucellPlane - crossComp');
    ucellgeom.runCurrent;
    
    % fillet
    if filRad > 0
        ucellFillet = ucellWP.geom.create('ucellFillet', 'Fillet');
        ucellFillet.set('radius', filRad);
        
        % selecting vertices to fillet - vertices that do not lie on
        % periodic boundary edges
        % 1. get edges along periodic boundaries
        g = model.geom('ucell').obj('ucellWP');
        edgXL = edgeindex(g,[-a/2,   0,-th/2],[1,0,0]); % X left edge
        edgXR = edgeindex(g,[ a/2,   0,-th/2],[1,0,0]); % X right edge
        edgYT = edgeindex(g,[   0, w/2,-th/2],[0,1,0]); % Y top edge
        edgYB = edgeindex(g,[   0,-w/2,-th/2],[0,1,0]); % Y bottom edge
        edgPBCs = [edgXL,edgXR,edgYT,edgYB];
        
        % 2. find vertices connected by edges 
        vtx_adjto_edg = g.getAdj(1,0);
        vtxPBCs = [];
        for edg = edgPBCs
            vtxPBCs = [vtxPBCs;vtx_adjto_edg{edg+1}]; % add 1 to edge index as first element in vtx_adjto_edg is {}
        end
        vtxPBCs = unique(vtxPBCs);
        
        % 3. exclude these vertices from fillet
        vtxFil = setdiff(1:g.getNVertices(),vtxPBCs);
        ucellFillet.selection('point').set('ucellComp', vtxFil);
        ucellgeom.runCurrent;
    end
end

%% Extrude unit cell
ucellPlaneExt = ucellgeom.feature.create('ucellPlaneExt', 'Extrude');
ucellPlaneExt.set('distance', th);
ucellgeom.runCurrent;
ucellName = 'ucellPlaneExt';

%% Symmetries in y and z
if evenz
    % create symmetry block
    symZ = abs(evenz);
    ucsymZWP = ucellgeom.feature.create('ucsymZWP', 'WorkPlane');
    ucsymZWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
    ucsymZPlane = ucsymZWP.geom.feature.create('ucsymZPlane', 'Rectangle');
    ucsymZPlane.set('type', 'solid').set('base', 'corner');
    ucsymZPlane.set('pos', [-a/2 -w/2]).set('size', [a w]);

    % extrude symmetry block
    ucellgeom.runCurrent;
    ucsymZPlaneExt = ucellgeom.feature.create('ucsymZPlaneExt', 'Extrude');
    ucsymZPlaneExt.set('distance', th/2^symZ);
    
    % compose: unit cell - symmetry block
    ucsymZComp = ucellgeom.feature.create('ucsymZComp', 'Compose');
    ucsymZComp.selection('input').set(ucellName);
    ucsymZComp.selection('input').set('ucsymZPlaneExt');
    ucsymZComp.set('formula', [ucellName,' - ucsymZPlaneExt']);
    ucellgeom.run;
    ucellName = 'ucsymZComp';
end

%% Plot
if P.plotgeom
    ucellgeom.run;
    figure;
    subplot(1,2,2)
    mphgeom(model,'ucell');
    daspect([1 1 1]);
    view(3)     %default 3D view in MATLAB
    camzoom(1.15)
    axis off
    
    subplot(1,2,1)
    mphgeom(model,'ucell');
    daspect([1 1 1]);
    view(0,90)     %XY view in MATLAB
    camzoom(0.8)
    axis off
    
    if P.filRad > 0
        filRadTtl = [', Rf = ',num2str(P.filRad*1e9,'%.1f'),'nm'];
    else
        filRadTtl = '';
    end
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,{['\bf ',P.celltype,' cross, ',...
                  'a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                  'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
                  'th = ',num2str(P.th*1e9,'%.0f'),'nm'];
                 ['cx = ',num2str(P.cx*1e9,'%.0f'),'nm, ', ...
                  'cy = ',num2str(P.cy*1e9,'%.0f'),'nm, ',...
                  'tx = ',num2str(P.tx*1e9,'%.0f'),'nm, ' ...
                  'ty = ',num2str(P.ty*1e9,'%.0f'),'nm',...
                  filRadTtl]},...
            'HorizontalAlignment' ,'center','VerticalAlignment', 'top')
end

end