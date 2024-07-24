function [model,P] = buildLowerBoomerangUnitCell(model,P)
%
% buildBoomerangUnitCell.m
%
% Model exported on Jul 21 2024, 16:57 by COMSOL 6.1.0.357.

%% read the input parameters 
a = P.a;        % lattice constant 
w = P.w;        % the width of the hole 
th = P.th;        % the height of the hole
r = P.r;        % the height of the hole 
h = P.h;        % the height of the hole in the lower portion
d = P.d;        % the width of the hole in the lower portion
r1 = P.r1;      % the fillet radius of the edges of the hole 
r2 = P.r2;      % the fillet radius of the center of the hole 
abssym = abs(P.mbevenz);    % symmetry in the z direction 

%% Create component 
ucellcomp = model.modelNode.create('comp1');
ucellcomp.label('Unit cell FEM simulation');
ucelllabel = 'Unit Cell';
ucellname = 'geom1';
ucellgeom = model.geom.create(ucellname, 3);
ucellgeom.label(ucelllabel);

%% create unit cell with boomerang geometry in the lower cavity
ucellWP = ucellgeom.create('wp1', 'WorkPlane');
ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
ucellplane = ucellWP.geom.feature.create('sq1', 'Square');
ucellplane.set('pos', [-a/2 -a/2]);
ucellplane.set('size', a);
rec1 = ucellWP.geom.feature.create('r1', 'Rectangle');
rec1.set('pos', [a/2-a/2 r/2+a/2+sqrt(3)*w/4-a/2]);
rec1.set('base', 'center');
rec1.set('size', [w r]);
rec2 = ucellWP.geom.feature.create('r2', 'Rectangle');
rec2.set('pos', [a/2-a/2 a/2-sqrt(3)*w/2+sqrt(3)*w/4-a/2]);
rec2.set('rot', 120);
rec2.set('size', [w r]);
rec3 = ucellWP.geom.feature.create('r3', 'Rectangle');
rec3.set('pos', [a/2+w/2-a/2 a/2+sqrt(3)*w/4-a/2]);
rec3.set('rot', 240);
rec3.set('size', [w r]);
rec4 = ucellWP.geom.feature.create('r4', 'Rectangle');
rec4.set('pos', [w/2+sqrt(3)*r/2+a/2-a/2 a/2-r/2+sqrt(3)*w/4-a/2]);
rec4.set('rot', 180);
rec4.set('size', [d h]);
rec5 = ucellWP.geom.feature.create('r5', 'Rectangle');
rec5.set('pos', [-w/2-sqrt(3)*r/2+a/2+d-a/2 a/2-r/2+sqrt(3)*w/4-a/2]);
rec5.set('rot', 180);
rec5.set('size', [d h]);
pol1 = ucellWP.geom.feature.create('pol1', 'Polygon');
pol1.set('tableconstr', {'off' 'off'});
pol1.set('source', 'table');
pol1.set('table', [-w/2+a/2-a/2 a/2+sqrt(3)*w/4-a/2; w/2+a/2-a/2 a/2+sqrt(3)*w/4-a/2; a/2-a/2 a/2-sqrt(3)*w/4-a/2; -w/2+a/2-a/2 a/2+sqrt(3)*w/4-a/2]);
compose1 = ucellWP.geom.feature.create('co1', 'Compose');
compose1.set('formula', 'sq1-r1-r2-r3-r4-r5-pol1');
fillet1 = ucellWP.geom.feature.create('fil1', 'Fillet');
fillet1.set('radius', r1);
fillet1.selection('point').set('co1(1)', [3 5 8 11 12 14]);
fillet2 = ucellWP.geom.feature.create('fil2', 'Fillet');
fillet2.set('radius', r2);
fillet2.selection('point').set('fil1(1)', [4 8 9 12 14 17 21]);
extrude = model.component('comp1').geom('geom1').feature.create('ext1', 'Extrude');
extrude.setIndex('distance', th, 0);
extrude.selection('input').set({'wp1'});

%% Symmetry in z 
if abs(P.mbevenz)
    symZth = th/2;
    symW = a;
    % create symmetry block
    symZWP = ucellgeom.feature.create('symZWP', 'WorkPlane');
    symZWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -symZth);
    symZPlane = symZWP.geom.feature.create('sq2', 'Square');
    symZPlane.set('pos', [-a/2 -a/2]);
    symZPlane.set('size', a);

    % extrude symmetry block
    ucellgeom.runCurrent;
    symZPlaneExt = ucellgeom.feature.create('symZPlaneExt', 'Extrude');
    symZPlaneExt.set('distance', symZth);

    % compose: unit cell - symmetry block
    symZComp = ucellgeom.feature.create('symZComp', 'Compose');
    symZComp.selection('input').set('ext1');
    symZComp.selection('input').set('symZPlaneExt');
    symZComp.set('formula', ['ext1 - symZPlaneExt']);
    ucellgeom.runCurrent;


    % beam z-symmetry plane
    delta = 10e-9;
    ZsymSel = ucellgeom.create('ZsymSel', 'BoxSelection');
    ZsymSel.set('xmin', -a/2-delta).set('xmax', a/2+delta);
    ZsymSel.set('ymin', -a/2-delta).set('ymax', a/2+delta);
    ZsymSel.set('zmin', -delta).set('zmax', delta);
    ZsymSel.set('entitydim', 2).set('condition', 'allvertices'); % only want beam surface, exclude air holes
    ucellgeom.runCurrent;
    inds = model.selection([ucellname,'_ZsymSel']).inputEntities();
    P.bndSel.Zsym = inds';

end
%% Making selections (manual)
mphgeom(model);
P.xEnd1 =  bndindex(ucellgeom, [-a/2 0 0], [1 0 0]);
P.xEnd2 = bndindex(ucellgeom, [ a/2 0 0], [1 0 0]);

P.yEnd1 =  bndindex(ucellgeom, [0 -a/2 0], [0 1 0]);
P.yEnd2 = bndindex(ucellgeom, [0 a/2 0], [0 1 0]);
% Note that this will return no indices if there is no boundary at z=0
P.zEnd = bndindex(ucellgeom, [0 0 0], [0 0 1]);
disp(P) % debugging
out = model;
