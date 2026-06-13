% Function to construct geometry (nanobeam, air/mech PML etc.) for FEM
% simulations
% Cleaven Chia, 11/21/18

function [model,P] = BuildNanobeamFEM(model,P)

% extract geometry parameters from P
wid = P.w;
thi = P.th;
% asym = P.asym;

if isfield(P,'asymCav') && P.asymCav
    geom = P.geom;
    beamLen = P.beamLen;
    P.xc = P.beamLenHalfL;
else
    geom = P.geomHalf;
    beamLen = P.beamLenHalf;
    P.xc = 0;
end
   

nholes = size(geom,1);
hx = geom(:,1)';
hy = geom(:,2)';
xpos = geom(:,3)';
ypos = geom(:,4)';
w = geom(:,5)';
a = geom(:,6)';

% % if using rectangular cross section, use half of height due to symmetry
% if strcmp(P.xsect,'rect')
%     thi = P.th/2;
% end
% 
% if using triangular cross section, recalculate height of beam
if strcmp(P.xsect,'tri')
    thi = wid/(2*tan(P.theta*pi/180));
    P.th = thi;
end

%% Create component
comp = model.modelNode.create('comp');
comp.label('Nanobeam FEM simulation');
% if P.solveOpt
%     beamlabel = 'Nanobeam with air cylinder';
%     beamname = 'beamCyl';
% else
%     beamlabel = 'Nanobeam';
%     beamname = 'beam';
% end

P.geomname = 'beam';
beamgeom = model.geom.create(P.geomname, 3);
beamgeom.label('Nanobeam geometry');

%% Create beam with rectangular cross-section
beamWP = beamgeom.feature.create('beamWP', 'WorkPlane');
beamWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -thi/2);


if isfield(P,'celltype') && strcmp(P.celltype,'blockTet')
    % tethers
    beamplane = beamWP.geom.feature.create('beamplane', 'Rectangle');
    beamplane.set('type', 'solid').set('base', 'corner');
    beamplane.set('pos', [0 0]).set('size', [beamLen P.hy/2]);
    beamgeom.runCurrent;
    
    holeList = {};
    holeFormula = [];
    for k = 1:nholes
        if k==1 && ~P.holeatctr
            hxk = hx(k)/2;
            dx = 0;
        else
            hxk = hx(k);
            dx = hxk/2;
        end
    %     hole = holes.geom.feature.create(['hole_' num2str(k)], 'Ellipse');
        holeID = ['block_' num2str(k)];
        hole = beamWP.geom.feature.create(holeID, 'Rectangle');
        hole.set('type', 'solid').set('base', 'corner');
        hole.set('pos', [xpos(k)-dx ypos(k)]).set('size', [hxk w(k)/2]);
        holeList = [holeList,holeID];
        holeFormula = [holeFormula,' + ',holeID];
    end
    beamgeom.runCurrent;
    
    
    % compose workplane
    beamComp = beamWP.geom.feature.create('beamComp', 'Compose');
    beamComp.selection('input').set(['beamplane',holeList]);
    beamComp.set('formula', ['beamplane',holeFormula]).set('intbnd',false);
    beamgeom.runCurrent;
    
else
    
    beamplane = beamWP.geom.feature.create('beamplane', 'Rectangle');
    beamplane.set('type', 'solid').set('base', 'corner');
    beamplane.set('pos', [0 0]).set('size', [beamLen wid/2]);
    beamgeom.runCurrent;

    holeList = {};
    holeFormula = [];
    for k = 1:nholes
    %     hole = holes.geom.feature.create(['hole_' num2str(k)], 'Ellipse');
        holeID = ['hole_' num2str(k)];
        hole = beamWP.geom.feature.create(holeID, 'Ellipse');
        hole.set('pos', [xpos(k) ypos(k)]);
        hole.set('semiaxes', [hx(k)/2 hy(k)/2]);
        holeList = [holeList,holeID];
        holeFormula = [holeFormula,' - ',holeID];
    end
    beamgeom.runCurrent;

    % compose workplane
    beamComp = beamWP.geom.feature.create('beamComp', 'Compose');
    beamComp.selection('input').set(['beamplane',holeList]);
    beamComp.set('formula',['beamplane',holeFormula]).set('intbnd',false);
    beamgeom.runCurrent;    
end

% extrude
finBeamTag = 'beamHoles';
beamHoles = beamgeom.feature.create(finBeamTag, 'Extrude');
beamHoles.set('distance', thi);
beamgeom.runCurrent;  
displayBeamStr = 'Nanobeam, rectangular cross-section';

% track max dimensions for selections
totLen = beamLen;
maxWid = wid;
maxThi = thi;

%% Create beam with triangular cross-section
% create workplane for triangular cross-section
if strcmp(P.xsect,'tri')
    P.mevenz = 0;
    beamTriWP = beamgeom.feature.create('beamTriWP', 'WorkPlane');
    beamTriWP.set('planetype', 'coordinates');
    beamTriWP.set('genpoints', [0 0 0; 0 1 0; 0 0 1]);
    beamTri = beamTriWP.geom.feature.create('beamTri', 'Polygon');
    beamTri.set('type', 'solid');
    beamTri.set('x', [0 0 wid/2]).set('y', [-thi/2 thi/2 thi/2]);
    
    beamgeom.runCurrent;    % run workplane geometry before extrusion
    beamTriExt = beamgeom.feature.create('beamTriExt', 'Extrude');
    beamTriExt.set('distance', beamLen);
    
    % Compose final beam with holes and triangular cross-section: (beam - holes) * triangle cross-section
    finBeamTag = 'beamHolesTri';
    BeamHolesTri = beamgeom.feature.create(finBeamTag, 'Compose');
    BeamHolesTri.selection('input').set('beamTriExt');
    BeamHolesTri.selection('input').set('beamHoles');
    BeamHolesTri.set('formula', 'beamTriExt * beamHoles');
    beamgeom.runCurrent;  
    displayBeamStr = 'Nanobeam, triangular cross-section';
end

% track max dimensions for selections
totLen = beamLen;
maxWid = max(w);
maxThi = thi;

%% Create beam with cross-section more realistic from iso etch
% create workplane for triangular cross-section
if strcmp(P.xsect,'isoFit')
    P.mevenz = 0;
    beamTriWP = beamgeom.feature.create('beamisoFitWP', 'WorkPlane');
    beamTriWP.set('planetype', 'coordinates');
    beamTriWP.set('genpoints', [0 0 0; 0 1 0; 0 0 1]);
    beamTri = beamTriWP.geom.feature.create('beamIso', 'Polygon');
    beamTri.set('type', 'solid');
    beamTri.set('x', [0 wid/2 wid/2]).set('y', [-thi/2 -thi/2 -thi/2+sqrt(3)*wid/6]);
    
    beamgeom.runCurrent;    % run workplane geometry before extrusion
    beamTriExt = beamgeom.feature.create('beamIsoExt', 'Extrude');
    beamTriExt.set('distance', beamLen);
    
    % Compose final beam with holes and triangular cross-section: (beam - holes) * triangle cross-section
    finBeamTag = 'beamHolesIso';
    BeamHolesTri = beamgeom.feature.create(finBeamTag, 'Compose');
    BeamHolesTri.selection('input').set('beamIsoExt');
    BeamHolesTri.selection('input').set('beamHoles');
    BeamHolesTri.set('formula', 'beamHoles - beamIsoExt');
    beamgeom.runCurrent;  
    displayBeamStr = 'Nanobeam, isotropic fitted cross-section';
end

% track max dimensions for selections
totLen = beamLen;
maxWid = max(w);
maxThi = thi;

%% Run geometry
beamgeom.run;

%% to implement: adding phononic mirrors
% use createNanobeamGeom
% then update total length

%% Create air cylinder around beam
displayCylStr = '';
if isfield(P,'airrad') && P.solveOpt
    cyl_cut_wp = beamgeom.feature.create('cyl_cut_wp', 'WorkPlane');
    cyl_cut_wp.set('planetype', 'quick').set('quickplane', 'xz');
    cyl_cut = cyl_cut_wp.geom.feature.create('cyl_cut','Rectangle');
    cyl_cut.set('type','solid').set('pos',[0,0]).set('size',[totLen,P.airrad]);

    air_cyl = beamgeom.feature.create('air_cyl', 'Revolve');
    air_cyl.selection('input').set('cyl_cut_wp');
    
    air_cyl.set('angle1','-180');
    air_cyl.set('axis',[1,0]).set('pos',[0,0]).set('angle2','0');

    % Compose final geometry
    beamHolesAir = beamgeom.feature.create('beamHolesAir', 'Compose');
    beamHolesAir.selection('input').set(finBeamTag);
    beamHolesAir.selection('input').set('air_cyl');
    beamHolesAir.set('formula', [finBeamTag,' + air_cyl']);
    beamgeom.run;
    displayCylStr = ', air cylinder';
    finBeamTag = 'beamHolesAir';
    
    % track max dimensions for selections
    totLen = beamLen;
    maxWid = max(maxWid,2*P.airrad);
    maxThi = max(maxThi,2*P.airrad);
end



%% Create mechanical PML  (L-frame around the +x beam end)
xL = 0;
displayPMLStr = '';
if P.solveMech && isfield(P,'solveMechPML') && P.solveMechPML
    % end-cap arm: +x beyond beam end, lower part (touches y=0 symmetry plane)
    PMLxArm = beamgeom.feature.create('PMLxArm', 'Block');
    PMLxArm.set('base', 'corner');
    PMLxArm.set('pos', [beamLen, 0, -thi/2]);
    PMLxArm.set('size', [P.PMLLen, P.PMLLen, thi]);

    beamgeom.runCurrent;
    totLen = totLen + P.PMLLen;

    PMLtag = 'PMLxArm';

    % mirrored left end-cap for asymmetric cavity
    if isfield(P,'asymCav') && P.asymCav
        PMLxArmL = beamgeom.feature.create('PMLxArmL', 'Block');
        PMLxArmL.set('base', 'corner');
        PMLxArmL.set('pos', [-P.PMLLen, 0, -thi/2]);
        PMLxArmL.set('size', [P.PMLLen, P.PMLLen, thi]);

        PMLcornerL = beamgeom.feature.create('PMLcornerL', 'Block');
        PMLcornerL.set('base', 'corner');
        PMLcornerL.set('pos', [-P.PMLLen, max(w)/2, -thi/2]);
        PMLcornerL.set('size', [P.PMLLen, P.PMLLen, thi]);
        beamgeom.runCurrent;
        totLen = totLen + P.PMLLen;
        xL = -P.PMLLen;

        PMLAll = beamgeom.feature.create('MechPML', 'Compose');
        PMLAll.selection('input').set({'MechPMLR','PMLxArmL','PMLcornerL'});
        PMLAll.set('formula', 'MechPMLR + PMLxArmL + PMLcornerL');
        beamgeom.runCurrent;
        PMLtag = 'MechPML';
    end

    % Compose final geometry: beam + PML
    beamHolesAir = beamgeom.feature.create('beamHolesPML', 'Compose');
    beamHolesAir.selection('input').set(finBeamTag);
    beamHolesAir.selection('input').set(PMLtag);
    beamHolesAir.set('formula', [finBeamTag,' + ',PMLtag]);
    beamgeom.run;
    displayPMLStr = ', mech PML';
    finBeamTag = 'beamHolesPML';

    % track max dimensions for selections (y now reaches max(w)/2 + PMLLen)
    maxWid = max(maxWid, max(w) + 2*P.PMLLen);
    maxThi = max(maxThi, thi);
end

%% Symmetry in z
symZOn = strcmp(P.xsect,'rect') && ...
        ((P.solveMech && ~P.solveOpt && abs(P.mevenz)) || ...
         (P.solveOpt && ~P.solveMech && abs(P.oevenz)) || ...
         (P.solveMech && P.solveOpt && abs(P.mevenz) && abs(P.oevenz)));
if symZOn
%     symZthList = thi/2;
%     symWList = wid/2;
%     if isfield(P,'airrad') && P.solveOpt
%         symZthList(end+1) = P.airrad;
%         symWList(end+1) = P.airrad;
%     end
%     if P.solveMech && isfield(P,'solveMechPML') && P.solveMechPML
%         symZthList(end+1) = thi/2;
%         symWList(end+1) = P.PMLLen;
%     end
    % create symmetry block
    symZWP = beamgeom.feature.create('symZWP', 'WorkPlane');
    symZWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -maxThi/2);
    symZPlane = symZWP.geom.feature.create('symZPlane', 'Rectangle');
    symZPlane.set('type', 'solid').set('base', 'corner');
    symZPlane.set('pos', [xL 0]).set('size', [totLen maxWid/2]);
    
    % extrude symmetry block
    beamgeom.runCurrent;
    symZPlaneExt = beamgeom.feature.create('symZPlaneExt', 'Extrude');
    symZPlaneExt.set('distance', maxThi/2);
    
    % compose: unit cell - symmetry block
    symZComp = beamgeom.feature.create('symZComp', 'Compose');
    symZComp.selection('input').set(finBeamTag);
    symZComp.selection('input').set('symZPlaneExt');
    symZComp.set('formula', [finBeamTag,' - symZPlaneExt']);
    beamgeom.run;
    
    maxThi = maxThi/2;
end
    


display(['Geometry created - ',displayBeamStr,displayCylStr,displayPMLStr]);


%% Create domain selections after full geometry is constructed

% Create selection with beam only
% define box slightly larger and fully containing full beam volume
delta = 10e-9; 
beamSel = beamgeom.create('beamSel', 'BoxSelection');
beamSel.set('xmin', -delta).set('xmax', beamLen + delta);
beamSel.set('ymin', -delta).set('ymax', max(w)/2 + delta);
beamSel.set('zmin', -thi/2-delta).set('zmax', thi/2 + delta);
beamSel.set('entitydim', 3).set('condition', 'inside');
beamgeom.runCurrent;
inds = model.selection([P.geomname,'_beamSel']).inputEntities();
P.domSel.beam = inds';

% Create selection with cylinder only
% - define box slightly larger and fully containing full cylinder volume
% then take difference selection with beam
if isfield(P,'airrad') && P.solveOpt
    delta = 10e-9; 
    beamCylSel = beamgeom.create('beamCylSel', 'BoxSelection');
    beamCylSel.set('xmin', -delta).set('xmax', beamLen + delta);
    beamCylSel.set('ymin', -delta).set('ymax', P.airrad + delta);
    beamCylSel.set('zmin', -P.airrad-delta).set('zmax', P.airrad + delta);
    beamCylSel.set('condition', 'inside');
    beamgeom.runCurrent;
    
    cylSel = beamgeom.create('cylSel', 'DifferenceSelection');
    cylSel.set('entitydim',3).set('add','beamCylSel').set('subtract','beamSel');
    beamgeom.runCurrent;
    inds = model.selection([P.geomname,'_cylSel']).inputEntities();
    P.domSel.cyl = inds';
end

% Create selection with PML pads only (L-frame)
% L-shape cannot be captured by one 'inside' box; subtract beam from full footprint.
if P.solveMech && isfield(P,'solveMechPML') && P.solveMechPML
    delta = 10e-9;
    beamPMLSel = beamgeom.create('beamPMLSel', 'BoxSelection');
    beamPMLSel.set('xmin', xL-delta).set('xmax', beamLen+P.PMLLen+delta);
    beamPMLSel.set('ymin', -delta).set('ymax', P.PMLLen+delta);
    beamPMLSel.set('zmin', -thi/2-delta).set('zmax', thi/2+delta);
    beamPMLSel.set('entitydim', 3).set('condition', 'inside');
    beamgeom.runCurrent;

    PMLSel = beamgeom.create('PMLSel', 'DifferenceSelection');
    PMLSel.set('entitydim', 3).set('add', 'beamPMLSel').set('subtract', 'beamSel');
    beamgeom.runCurrent;
    inds = model.selection([P.geomname,'_PMLSel']).inputEntities();
    P.domSel.PML = inds';
    disp(['PML domain indices: ', num2str(P.domSel.PML)]);
end






%% Create boundary selections after full geometry is constructed
% flat planes: use box selection with condition that all vertices are in box
% curved planes: use box selection with condition that box intersects plane

delta = 5e-9;

% beam x-symmetry plane
beamXsymSel = beamgeom.create('beamXsymSel', 'BoxSelection');
beamXsymSel.set('xmin', -delta).set('xmax', delta);
beamXsymSel.set('ymin', -delta).set('ymax', max(w)/2+delta);
beamXsymSel.set('zmin', -thi/2*(1-symZOn)-delta).set('zmax', thi/2+delta);
beamXsymSel.set('entitydim', 2).set('condition', 'allvertices');
beamgeom.runCurrent;
inds = model.selection([P.geomname,'_beamXsymSel']).inputEntities();
P.bndSel.beamXsym = inds';

% beam x-end plane
beamXendSel = beamgeom.create('beamXendSel', 'BoxSelection');
beamXendSel.set('xmin', beamLen-delta).set('xmax', beamLen+delta);
beamXendSel.set('ymin', -delta).set('ymax', max(w)/2+delta);
beamXendSel.set('zmin', -thi/2*(1-symZOn)-delta).set('zmax', thi/2+delta);
beamXendSel.set('entitydim', 2).set('condition', 'allvertices');
beamgeom.runCurrent;
inds = model.selection([P.geomname,'_beamXendSel']).inputEntities();
P.bndSel.beamXend = inds';

% beam y-symmetry plane
beamYsymSel = beamgeom.create('beamYsymSel', 'BoxSelection');
beamYsymSel.set('xmin', -delta).set('xmax', beamLen+delta);
beamYsymSel.set('ymin', -delta).set('ymax', delta);
beamYsymSel.set('zmin', -thi/2*(1-symZOn)-delta).set('zmax', thi/2+delta);
beamYsymSel.set('entitydim', 2).set('condition', 'allvertices');
beamgeom.runCurrent;
inds = model.selection([P.geomname,'_beamYsymSel']).inputEntities();
P.bndSel.beamYsym = inds';

% beam z-symmetry plane
beamZsymSel = beamgeom.create('beamZsymSel', 'BoxSelection');
beamZsymSel.set('xmin', delta).set('xmax', 2*delta);
beamZsymSel.set('ymin', min(w)/2-2*delta).set('ymax', min(w)/2-delta);
beamZsymSel.set('zmin', -delta).set('zmax', delta);
beamZsymSel.set('entitydim', 2).set('condition', 'intersects'); % only want beam surface, exclude air holes
beamgeom.runCurrent;
inds = model.selection([P.geomname,'_beamZsymSel']).inputEntities();
P.bndSel.beamZsym = inds';

if isfield(P,'airrad') && P.solveOpt
    % cyl x-symmetry plane
    cylXsymSel = beamgeom.create('cylXsymSel', 'BoxSelection');
    cylXsymSel.set('xmin', -delta).set('xmax', delta);
    cylXsymSel.set('ymin', -delta).set('ymax', P.airrad+delta);
    cylXsymSel.set('zmin', -P.airrad-delta).set('zmax', P.airrad+delta);
    cylXsymSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;
    inds = model.selection([P.geomname,'_cylXsymSel']).inputEntities();
    P.bndSel.cylXsym = inds';
        
    % cyl x-end plane
    cylXendSel = beamgeom.create('cylXendSel', 'BoxSelection');
    cylXendSel.set('xmin', beamLen-delta).set('xmax', beamLen+delta);
    cylXendSel.set('ymin', -delta).set('ymax', P.airrad+delta);
    cylXendSel.set('ymin', -delta).set('ymax', P.airrad+delta);
    cylXendSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;
    inds = model.selection([P.geomname,'_cylXendSel']).inputEntities();
    P.bndSel.cylXend = inds';
    
    % cyl y-symmetry plane
    cylYsymSel = beamgeom.create('cylYsymSel', 'BoxSelection');
    cylYsymSel.set('xmin', -delta).set('xmax', beamLen+delta);
    cylYsymSel.set('ymin', -delta).set('ymax', delta);
    cylYsymSel.set('zmin', -P.airrad*(1-symZOn)-delta).set('zmax', P.airrad+delta);
    cylYsymSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;
    inds = model.selection([P.geomname,'_cylYsymSel']).inputEntities();
    P.bndSel.cylYsym = inds';
    
    % cyl z-symmetry plane
    cylZsymSel = beamgeom.create('cylZsymSel', 'BoxSelection');
    cylZsymSel.set('xmin', -delta).set('xmax', beamLen+delta);
    cylZsymSel.set('ymin', -delta).set('ymax', P.airrad+delta);
    cylZsymSel.set('zmin', -delta).set('zmax', delta);
    cylZsymSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;
    inds = model.selection([P.geomname,'_cylZsymSel']).inputEntities();
    P.bndSel.cylZsym = inds';
    
%     % beam holes z-symmetry plane
%     bhlZsymSel = beamgeom.create('bhlZsymSel', 'BoxSelection');
%     bhlZsymSel.set('xmin', -delta).set('xmax', beamLen+delta);
%     bhlZsymSel.set('ymin', -delta).set('ymax', max(P.geom(:,2))/2+delta);
%     bhlZsymSel.set('zmin', -delta).set('zmax', delta);
%     bhlZsymSel.set('entitydim', 2).set('condition', 'allvertices');
%     beamgeom.runCurrent;
%     inds = model.selection([P.geomname,'_bhlZsymSel']).inputEntities();
%     P.bndSel.cylZsym = [P.bndSel.cylZsym,inds'];
    
    % cyl top curved plane
    cylCurvSel = beamgeom.create('cylCurvSel', 'BoxSelection');
    cylCurvSel.set('xmin', delta).set('xmax', 2*delta);
    cylCurvSel.set('ymin', P.airrad/sqrt(2)-10*delta).set('ymax', P.airrad/sqrt(2)+10*delta);
    cylCurvSel.set('zmin', P.airrad/sqrt(2)-10*delta).set('zmax', P.airrad/sqrt(2)+10*delta);
    cylCurvSel.set('entitydim', 2).set('condition', 'intersects');
    beamgeom.runCurrent;
    inds = model.selection([P.geomname,'_cylCurvSel']).inputEntities();
    P.bndSel.cylCurv = inds';
    
    if ~strcmp(P.xsect,'rect')
        cylCurv2Sel = beamgeom.create('cylCurv2Sel', 'BoxSelection');
        cylCurv2Sel.set('xmin', delta).set('xmax', 2*delta);
        cylCurv2Sel.set('ymin', P.airrad/sqrt(2)-10*delta).set('ymax', P.airrad/sqrt(2)+10*delta);
        cylCurv2Sel.set('zmin', -P.airrad/sqrt(2)-10*delta).set('zmax', -P.airrad/sqrt(2)+10*delta);
        cylCurv2Sel.set('entitydim', 2).set('condition', 'intersects');
        beamgeom.runCurrent;
        inds = model.selection([P.geomname,'_cylCurv2Sel']).inputEntities();
        P.bndSel.cylCurv = [P.bndSel.cylCurv,inds'];
    end
end

if P.solveMech && isfield(P,'solveMechPML') && P.solveMechPML
    % PML top face (z = +thi/2): whole L-frame footprint
    PMLtopSel = beamgeom.create('PMLtopSel', 'BoxSelection');
    PMLtopSel.set('xmin', xL-delta).set('xmax', beamLen+P.PMLLen+delta);
    PMLtopSel.set('ymin', -delta).set('ymax', max(w)/2+delta);
    PMLtopSel.set('zmin', thi/2-delta).set('zmax', thi/2+delta);
    PMLtopSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;
    inds = model.selection([P.geomname,'_PMLtopSel']).inputEntities();
    P.bndSel.PMLZtop = inds';

    % PML bottom face (z = -thi/2*(1-symZOn)): z-symmetry BC plane
    PMLbotSel = beamgeom.create('PMLbotSel', 'BoxSelection');
    PMLbotSel.set('xmin', xL-delta).set('xmax', beamLen+P.PMLLen+delta);
    PMLbotSel.set('ymin', -delta).set('ymax', max(w)/2+delta);
    PMLbotSel.set('zmin', -thi/2*(1-symZOn)-delta).set('zmax', -thi/2*(1-symZOn)+delta);
    PMLbotSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;
    inds = model.selection([P.geomname,'_PMLbotSel']).inputEntities();
    P.bndSel.PMLZsym = inds';

    % PML y-symmetry plane (y = 0): end-cap faces only
    PMLYsymSel = beamgeom.create('PMLYsymSel', 'BoxSelection');
    PMLYsymSel.set('xmin', xL-delta).set('xmax', beamLen+P.PMLLen+delta);
    PMLYsymSel.set('ymin', -delta).set('ymax', delta);
    PMLYsymSel.set('zmin', -thi/2*(1-symZOn)-delta).set('zmax', thi/2+delta);
    PMLYsymSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;
    inds = model.selection([P.geomname,'_PMLYsymSel']).inputEntities();
    P.bndSel.PMLYsym = inds';

    % outer +x face at x = beamLen+PMLLen (traction-free by default)
    PMLxoutSel = beamgeom.create('PMLxoutSel', 'BoxSelection');
    PMLxoutSel.set('xmin', beamLen+P.PMLLen-delta).set('xmax', beamLen+P.PMLLen+delta);
    PMLxoutSel.set('ymin', -delta).set('ymax', max(w)/2+delta);
    PMLxoutSel.set('zmin', -thi/2*(1-symZOn)-delta).set('zmax', thi/2+delta);
    PMLxoutSel.set('entitydim', 2).set('condition', 'allvertices');
    beamgeom.runCurrent;
    inds = model.selection([P.geomname,'_PMLxoutSel']).inputEntities();
    P.bndSel.PMLcurv = inds';

    % outer -x face at x = -PMLLen (asymCav left end-cap only)
    if isfield(P,'asymCav') && P.asymCav
        PMLxoutLSel = beamgeom.create('PMLxoutLSel', 'BoxSelection');
        PMLxoutLSel.set('xmin', -P.PMLLen-delta).set('xmax', -P.PMLLen+delta);
        PMLxoutLSel.set('ymin', -delta).set('ymax', max(w)/2+P.PMLLen+delta);
        PMLxoutLSel.set('zmin', -thi/2*(1-symZOn)-delta).set('zmax', thi/2+delta);
        PMLxoutLSel.set('entitydim', 2).set('condition', 'allvertices');
        beamgeom.runCurrent;
        inds = model.selection([P.geomname,'_PMLxoutLSel']).inputEntities();
        P.bndSel.PMLcurv = [P.bndSel.PMLcurv, inds'];
    end
end






















