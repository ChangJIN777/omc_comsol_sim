
function [model,P] = DrawUnitCell2DSqTets(model,P)

% extract parameters from P
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

% full name for use in plot titles
P.celltypeTitle = 'Square with tethers';

% lithographic constraints
if cx >= a || cy >= w
    error('block dimensions larger than unit cell dimensions')
end
if tx >= cx || ty >= cy
    error('linkage dimensions larger than block dimensions')
end

% calculate maximum permissible fillet radius
if filRad > 0
    % calculate extra lengths in tether due to protrusion of square
    px = ty/cy*cx;   % using similar triangles: cx/cy = px/ty
    py = tx/cx*cy;   % using similar triangles: cy/cx = py/tx
    
    l1 = 0.5*(a-cx+px); % lengths of tethers
    l2 = 0.5*(w-cy+py);
    l3 = 0.5*sqrt((cx-tx-px)^2 + (cy-ty-py)^2);   % half of length of sloped side
    
    % get maximum length of edge that can be filleted
    s1 = min([l1,l3]);
    s2 = min([l2,l3]);
    
    % determine max angle subtended by fillet
    filTheta1 = pi/2 + atan(cx/cy);  % angle between horizontal tether and edge of square
    filTheta2 = pi/2 + atan(cy/cx);  % angle between vertical tether and edge of square
    filAlpha1 = pi - filTheta1;
    filAlpha2 = pi - filTheta2;
    
    % get max fillet radius
    r1 = s1/tan(filAlpha1/2);
    r2 = s2/tan(filAlpha2/2);
    minFilRad = min([r1,r2]);
    
    if filRad > minFilRad
        error(['Fillet radius of ',...
               num2str(filRad*1e9,'%.1f'),'nm too large, ',...
               'maximum permissible radius is ',...
               num2str(minFilRad*1e9,'%.3f'),'nm'])
    end
end


%% Create square body
ucellgeom = model.geom.create('ucell', 3);
% create workplane for unit cell
ucellWP = ucellgeom.feature.create('ucellWP', 'WorkPlane');
ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
CSq = ucellWP.geom.feature.create('CSq', 'Polygon');
CSq.set('x', 0.5*[0 cx 0 -cx]).set('y', 0.5*[cy 0 -cy 0]);

%% Create linkages
Llink = ucellWP.geom.create('Llink', 'Rectangle');
Llink.set('type', 'solid').set('base', 'corner');
Llink.set('pos', [-a/2 -ty/2]).set('size', [a/2 ty]);
ucellgeom.runCurrent;

Rlink = ucellWP.geom.create('Rlink', 'Rectangle');
Rlink.set('type', 'solid').set('base', 'corner');
Rlink.set('pos', [0 -ty/2]).set('size', [a/2 ty]);
ucellgeom.runCurrent;

Ulink = ucellWP.geom.create('Ulink', 'Rectangle');
Ulink.set('type', 'solid').set('base', 'corner');
Ulink.set('pos', [-tx/2 0]).set('size', [tx w/2]);
ucellgeom.runCurrent;

Blink = ucellWP.geom.create('Blink', 'Rectangle');
Blink.set('type', 'solid').set('base', 'corner');
Blink.set('pos', [-tx/2 -w/2]).set('size', [tx w/2]);
ucellgeom.runCurrent;

% compose unit cell plane
ucellComp = ucellWP.geom.feature.create('ucellComp', 'Compose');
ucellComp.selection('input').set({'CSq' 'Llink' 'Rlink' 'Ulink' 'Blink'});
ucellComp.set('formula', 'CSq + Llink + Rlink + Ulink + Blink');
ucellComp.set('intbnd', 'off');
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
    text(0.5, 1,{['\bf ',P.celltypeTitle,', ',...
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