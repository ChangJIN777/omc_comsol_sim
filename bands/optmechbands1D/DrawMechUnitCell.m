% This script creates a single unit cell of triangular nanobeam
% Modified from previous version (Sean Meenehan, 01/26/11) for
% compatibility with COMSOL 5.x.
% Input arguments: model, P
% - model: COMSOL model
% - P: data structure consisting of following geometry parameters
%      a: lattice constant;  w: beam width; th: beam height;
%      hx: hole height; hy: hole width; 
%      hole: hole type ('rect' for rectangular, or 'elps' for elliptical); 
%      optBands: if included, create unit cell for optical bandstructure
%      calculations;
%      airrad: radius of air cylinder for optical bandstructure
%      calculations
% Output: model
% - model: updated COMSOL model with built nanobeam geometry 
%
% Cleaven Chia, 09/06/16

function [model,P] = DrawUnitCell(model,P,mfem)

nperiod = P.nperiod;
holeatedge = P.holeatedge; % 1/0 if unit cell terminates in middle of hole/dielectric
nholes = nperiod + holeatedge;

w = P.w;
th = P.th;
a = P.a;
hx = P.hx;
hy = P.hy;
P.beamLen = nperiod*a;
evenz = mfem.mbevenz;

if strcmp(P.xsect,'tri')
    P.mevenz = 0;
end

if P.a-P.hx < 50e-9 || P.w-P.hy < 50e-9
    error('Hole height or width too large');
end

%% Create component
ucellcomp = model.modelNode.create('ucell');
ucellcomp.label('Unit cell FEM simulation');
ucelllabel = 'Unit Cell';
ucellname = 'ucell';

ucellgeom = model.geom.create(ucellname, 3);
ucellgeom.label(ucelllabel);

%% Create unit cell with rectangular cross-section
ucellWP = ucellgeom.feature.create('ucellWP', 'WorkPlane');
ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
ucellplane = ucellWP.geom.feature.create('ucellplane', 'Rectangle');
ucellplane.set('type', 'solid').set('base', 'corner');
ucellplane.set('pos', [0 0]).set('size', [P.beamLen w/2]);
ucellgeom.runCurrent;

% create holes
holeList = {};
holeFormula = [];
for k = 1:nholes
%     hole = holes.geom.feature.create(['hole_' num2str(k)], 'Ellipse');
    holeID = ['hole_' num2str(k)];
    hole = ucellWP.geom.feature.create(holeID, 'Ellipse');
    hole.set('pos', [(1-holeatedge)*a/2+(k-1)*a 0]);
    hole.set('semiaxes', [hx/2 hy/2]);
    holeList = [holeList,holeID];
    holeFormula = [holeFormula,' - ',holeID];
end
ucellgeom.runCurrent;

% compose workplane
ucellComp = ucellWP.geom.feature.create('ucellComp', 'Compose');
ucellComp.selection('input').set(['ucellplane',holeList]);
ucellComp.set('formula', ['ucellplane',holeFormula]);
ucellgeom.runCurrent;    

% extrude
finBeamTag = 'ucellHoles';
ucellHoles = ucellgeom.feature.create(finBeamTag, 'Extrude');
ucellHoles.set('distance', th);
ucellgeom.runCurrent;  
displayBeamStr = 'Unit cell, rectangular cross-section';

%% Create unit cell with triangular cross-section
% create workplane for triangular cross-section
if strcmp(P.xsect,'tri')
    ucellTriWP = ucellgeom.feature.create('ucellTriWP', 'WorkPlane');
    ucellTriWP.set('planetype', 'coordinates');
    ucellTriWP.set('genpoints', [0 0 0; 0 1 0; 0 0 1]);
    ucellTri = ucellTriWP.geom.feature.create('ucellTri', 'Polygon');
    ucellTri.set('type', 'solid');
    ucellTri.set('x', [0 0 w/2]).set('y', [-th/2 th/2 th/2]);
    
    ucellgeom.runCurrent;    % run workplane geometry before extrusion
    ucellTriExt = ucellgeom.feature.create('ucellTriExt', 'Extrude');
    ucellTriExt.set('distance', P.beamLen);
    
    % Compose final beam with holes and triangular cross-section: (beam - holes) * triangle cross-section
    finBeamTag = 'ucellHolesTri';
    ucellHolesTri = ucellgeom.feature.create(finBeamTag, 'Compose');
    ucellHolesTri.selection('input').set('ucellTriExt');
    ucellHolesTri.selection('input').set('ucellHoles');
    ucellHolesTri.set('formula', 'ucellTriExt * ucellHoles');
    ucellgeom.runCurrent;  
    displayBeamStr = 'Unit Cell, triangular cross-section';
end

ucellgeom.run;

%% Symmetry in z
if strcmp(P.xsect,'rect') && abs(mfem.mbevenz)
    symZth = th/2;
    symW = w/2;
    % create symmetry block
    symZWP = ucellgeom.feature.create('symZWP', 'WorkPlane');
    symZWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -symZth);
    symZPlane = symZWP.geom.feature.create('symZPlane', 'Rectangle');
    symZPlane.set('type', 'solid').set('base', 'corner');
    symZPlane.set('pos', [0 0]).set('size', [P.beamLen symW]);
    
    % extrude symmetry block
    ucellgeom.runCurrent;
    symZPlaneExt = ucellgeom.feature.create('symZPlaneExt', 'Extrude');
    symZPlaneExt.set('distance', symZth);
    
    % compose: unit cell - symmetry block
    symZComp = ucellgeom.feature.create('symZComp', 'Compose');
    symZComp.selection('input').set(finBeamTag);
    symZComp.selection('input').set('symZPlaneExt');
    symZComp.set('formula', [finBeamTag,' - symZPlaneExt']);
    ucellgeom.run;
    
    
    % beam z-symmetry plane
    delta = 10e-9;
    ZsymSel = ucellgeom.create('ZsymSel', 'BoxSelection');
    ZsymSel.set('xmin', -delta).set('xmax', nperiod*a+delta);
    ZsymSel.set('ymin', -delta).set('ymax', w/2+delta);
    ZsymSel.set('zmin', -delta).set('zmax', delta);
    ZsymSel.set('entitydim', 2).set('condition', 'allvertices'); % only want beam surface, exclude air holes
    ucellgeom.runCurrent;
    inds = model.selection([ucellname,'_ZsymSel']).inputEntities();
    P.bndSel.Zsym = inds';
    
end

%% Create selections for boundary conditions
delta = 10e-9;
% beam x-symmetry plane
XsymSel = ucellgeom.create('XsymSel', 'BoxSelection');
XsymSel.set('xmin', -delta).set('xmax', delta);
XsymSel.set('ymin', -delta).set('ymax', w/2+delta);
XsymSel.set('zmin', -th/2*(1-evenz)-delta).set('zmax', th/2+delta);
XsymSel.set('entitydim', 2).set('condition', 'allvertices');
ucellgeom.runCurrent;
inds = model.selection([ucellname,'_XsymSel']).inputEntities();
P.bndSel.Xsym = inds';

% beam x-end plane
XendSel = ucellgeom.create('XendSel', 'BoxSelection');
XendSel.set('xmin', nperiod*a-delta).set('xmax', nperiod*a+delta);
XendSel.set('ymin', -delta).set('ymax', w/2+delta);
XendSel.set('zmin', -th/2*(1-evenz)-delta).set('zmax', th/2+delta);
XendSel.set('entitydim', 2).set('condition', 'allvertices');
ucellgeom.runCurrent;
inds = model.selection([ucellname,'_XendSel']).inputEntities();
P.bndSel.Xend = inds';

% beam y-symmetry plane
YsymSel = ucellgeom.create('YsymSel', 'BoxSelection');
YsymSel.set('xmin', -delta).set('xmax', nperiod*a+delta);
YsymSel.set('ymin', -delta).set('ymax', delta);
YsymSel.set('zmin', -th/2*(1-evenz)-delta).set('zmax', th/2+delta);
YsymSel.set('entitydim', 2).set('condition', 'allvertices');
ucellgeom.runCurrent;
inds = model.selection([ucellname,'_YsymSel']).inputEntities();
P.bndSel.Ysym = inds';
