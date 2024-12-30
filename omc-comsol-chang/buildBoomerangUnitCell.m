function [model,P] = buildBoomerangUnitCell(model,P)
%
% buildBoomerangUnitCell.m
%
% Model exported on Jul 21 2024, 16:57 by COMSOL 6.1.0.357.

%% read the input parameters 
a = P.a;        % lattice constant 
w = P.w;        % the width of the hole 
th = P.th;        % the height of the hole
r = P.r;        % the height of the hole 
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

%% create unit cell with boomerang geometry
ucellWP = ucellgeom.feature.create('wp1', 'WorkPlane');
ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
% ucellWP.set('unite', true);

ucellplane = ucellWP.geom.feature.create('pol1', 'Polygon');
ucellplane.label('Base plane');
ucellplane.set('source', 'table');
ucellplane.set('table', [0 0; a/2 (a/2)*sqrt(3); a*(3/2) (a/2)*sqrt(3); a 0; 0 0]);
rec_1 = ucellWP.geom.feature.create('r1', 'Rectangle');
rec_1.set('pos', [a*(1/2+1/4) a*sqrt(3)/4+r/2]);
rec_1.set('base', 'center');
rec_1.set('size', [w r]);
rec_2 = ucellWP.geom.feature.create('r2', 'Rectangle');
rec_2.set('pos', [a*(1/2+1/4)-sqrt(3)*r/4 a*sqrt(3)/4-r/4]);
rec_2.set('base', 'center');
rec_2.set('rot', 120);
rec_2.set('size', [w r]);
rec_3 = ucellWP.geom.feature.create('r3', 'Rectangle');
rec_3.set('pos', [a*(1/2+1/4)+sqrt(3)*r/4 a*sqrt(3)/4-r/4]);
rec_3.set('base', 'center');
rec_3.set('rot', 240);
rec_3.set('size', [w r]);
composit_geom = ucellWP.geom.feature.create('co1', 'Compose');
composit_geom.set('formula', 'pol1-r1-r2-r3');

% add fillets 
hole_pos = [a*(1/2+1/4) a*sqrt(3)/4];
selection_width = 50e-9;
addFillet(P,ucellWP.geom,hole_pos,selection_width)
% ucellWP.geom.create('fil1', 'Fillet');
% ucellWP.geom.feature('fil1').set('radius', r1);
% ucellWP.geom.feature('fil1').selection('point').set('r1(1)', [3 4]);
% ucellWP.geom.feature('fil1').selection('point').set('r3(1)', [3 4]);
% ucellWP.geom.feature('fil1').selection('point').set('r2(1)', [3 4]);
% centerTriangle = ucellWP.geom.feature.create('pol2', 'Polygon');
% centerTriangle.set('source', 'table');
% centerTriangle.set('table', [a/2 (w/2)*sqrt(3)/2+(a/2)/sqrt(3); w+a/2 (w/2)*sqrt(3)/2+(a/2)/sqrt(3); a/2+w/2 -(w/2)*sqrt(3)/2+(a/2)/sqrt(3); a/2 (w/2)*sqrt(3)/2+(a/2)/sqrt(3)]);

extrude = ucellgeom.feature.create('ext1', 'Extrude');
extrude.setIndex('distance', th, 0);
extrude.selection('input').set({'wp1'});
ucellgeom.runAll;

%% Symmetry in z 
if abs(P.mbevenz)
    symZth = th/2;
    symW = a;
    % create symmetry block
    symZWP = ucellgeom.feature.create('symZWP', 'WorkPlane');
    symZWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -symZth);
    symZPlane = symZWP.geom.feature.create('pol2', 'Polygon');
    symZPlane.set('source', 'table');
    symZPlane.set('table', [0 0; symW/2 (symW/2)*sqrt(3); symW*(3/2) (symW/2)*sqrt(3); symW 0; 0 0]);

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
P.xEnd1 =  bndindex(ucellgeom, [0 0 0], [sqrt(3)*a/2 -a/2 0]);
P.xEnd2 = bndindex(ucellgeom, [a 0 0], [sqrt(3)*a/2 -a/2 0]);

P.yEnd1 =  bndindex(ucellgeom, [0 0 0], [0 1 0]);
P.yEnd2 = bndindex(ucellgeom, [a (a/2)*sqrt(3) 0], [0 1 0]);
% Note that this will return no indices if there is no boundary at z=0
P.zEnd = bndindex(ucellgeom, [0 0 0], [0 0 1]);
% % debugging
% mphplot(model);
disp(P) % debugging
out = model;
end

function addFillet(P,ucellgeom,hole_pos,selection_width)
    w = P.w;
    r = P.r;
    r1 = P.r1;
    r2 = P.r2;
    disksel1_label = sprintf('h_disksel1');
    disksel2_label = sprintf('h_disksel2');
    fillet1_label = sprintf('h_fil1');
    fillet2_label = sprintf('h_fil2');
    % hole 
    disksel1 = ucellgeom.feature.create(disksel1_label, 'DiskSelection');
    disksel1.set('entitydim', 0);
    disksel1.set('posx', hole_pos(1));
    disksel1.set('posy', hole_pos(2));
    disksel1.set('r', w);
    disksel1.set('rin', w/(2*sqrt(2)));
    disksel1.set('condition', 'allvertices');
    fil1 = ucellgeom.feature.create(fillet1_label, 'Fillet');
    fil1.set('radius', r1);
    fil1.selection('point').named(disksel1_label);
    disksel2 = ucellgeom.feature.create(disksel2_label, 'DiskSelection');
    disksel2.set('entitydim', 0);
    disksel2.set('posy', hole_pos(2));
    disksel2.set('posx', hole_pos(1));
    disksel2.set('r', r+selection_width/2);
    disksel2.set('rin', r-selection_width/2);
    disksel2.set('condition', 'allvertices');
    fil2 = ucellgeom.feature.create(fillet2_label, 'Fillet');
    fil2.set('radius', r2);
    fil2.selection('point').named(disksel2_label);
end


