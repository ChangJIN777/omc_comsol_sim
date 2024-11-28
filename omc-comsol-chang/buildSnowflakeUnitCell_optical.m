function [model,P] = buildSnowflakeUnitCell_optical(model,P)
%
% snowFlakeBuildUnitCellCode.m
%
% Model exported on Jul 29 2024, 18:17 by COMSOL 6.1.0.357.

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

%% create unit cell with snowflake geometry
ucellWP = ucellgeom.feature.create('wp1', 'WorkPlane');
ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', 0);
% ucellWP.set('unite', true);

ucellplane = ucellWP.geom.feature.create('pol1', 'Polygon');
ucellplane.label('Base plane');
ucellplane.set('source', 'table');
ucellplane.set('table', [0 0; a/2 (a/2)*sqrt(3); a*(3/2) (a/2)*sqrt(3); a 0; 0 0]);

rec_1 = ucellWP.geom.feature.create('r1', 'Rectangle');
rec_1.set('pos', [a*(1/2+1/4) a*sqrt(3)/4]);
rec_1.set('base', 'center');
rec_1.set('size', [2*r w]);

rec_2 = ucellWP.geom.feature.create('r2', 'Rectangle');
rec_2.set('pos', [a*(1/2+1/4) a*sqrt(3)/4]);
rec_2.set('rot', 60);
rec_2.set('base', 'center');
rec_2.set('size', [2*r w]);

rec_3 = ucellWP.geom.feature.create('r3', 'Rectangle');
rec_3.set('pos', [a*(1/2+1/4) a*sqrt(3)/4]);
rec_3.set('rot', 120);
rec_3.set('base', 'center');
rec_3.set('size', [2*r w]);

composit_geom = ucellWP.geom.feature.create('co1', 'Compose');
composit_geom.set('formula', 'pol1-r1-r2-r3');

selection_width = 10e-9;
hole_pos = [a*(1/2+1/4) a*sqrt(3)/4];
addFillet(P,ucellWP,hole_pos,selection_width)

extrude = ucellgeom.feature.create('ext1', 'Extrude');
extrude.setIndex('distance', th/2, 0);
extrude.selection('input').set({'wp1'});
ucellgeom.runAll;

%% add the air block
airZth = 1e-6;
symW = a;
% create symmetry block
airWP = ucellgeom.feature.create('airWP', 'WorkPlane');
airWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', 0);
airZPlane = airWP.geom.feature.create('pol2', 'Polygon');
airZPlane.set('source', 'table');
airZPlane.set('table', [0 0; symW/2 (symW/2)*sqrt(3); symW*(3/2) (symW/2)*sqrt(3); symW 0; 0 0]);

% extrude symmetry block
ucellgeom.runCurrent;
airZPlaneExt = ucellgeom.feature.create('airZPlaneExt', 'Extrude');
airZPlaneExt.set('distance', airZth);

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

%% add the pml block
pmlZth = 0.5e-6;
% create symmetry block
pmlWP = ucellgeom.feature.create('pmlWP', 'WorkPlane');
pmlWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', airZth);
pmlZPlane = pmlWP.geom.feature.create('pol2', 'Polygon');
pmlZPlane.set('source', 'table');
pmlZPlane.set('table', [0 0; symW/2 (symW/2)*sqrt(3); symW*(3/2) (symW/2)*sqrt(3); symW 0; 0 0]);

% extrude symmetry block
ucellgeom.runCurrent;
pmlZPlaneExt = ucellgeom.feature.create('pmlZPlaneExt', 'Extrude');
pmlZPlaneExt.set('distance', pmlZth);

ucellgeom.runCurrent;

%% Making selections (with box select and index select)
mphgeom(model);

% box selection for the boundary condition 
x_boundary_disksel_r = ucellgeom.feature.create('x_boundary_boxsel_r', 'BoxSelection');
x_boundary_disksel_r.set('entitydim', 2);
x_boundary_disksel_r.set('xmin', 5*a/4-selection_width);
x_boundary_disksel_r.set('xmax', 5*a/4+selection_width);
x_boundary_disksel_r.set('ymin', sqrt(3)*a/4-selection_width);
x_boundary_disksel_r.set('ymax', sqrt(3)*a/4+selection_width);
x_boundary_disksel_r.set('zmin', th/2-selection_width);
x_boundary_disksel_r.set('zmax', airZth-selection_width);
x_boundary_disksel_r.set('inputent', 'all');
x_boundary_disksel_r.set('condition', 'intersects');

x_boundary_disksel_l = ucellgeom.feature.create('x_boundary_boxsel_l', 'BoxSelection');
x_boundary_disksel_l.set('entitydim', 2);
x_boundary_disksel_l.set('xmin', a/4-selection_width);
x_boundary_disksel_l.set('xmax', a/4+selection_width);
x_boundary_disksel_l.set('ymin', sqrt(3)*a/4-selection_width);
x_boundary_disksel_l.set('ymax', sqrt(3)*a/4+selection_width);
x_boundary_disksel_l.set('zmin', th/2-selection_width);
x_boundary_disksel_l.set('zmax', airZth-selection_width);
x_boundary_disksel_l.set('inputent', 'all');
x_boundary_disksel_l.set('condition', 'intersects');

ucellgeom.selection.create('xboundaries','CumulativeSelection');
ucellgeom.selection('xboundaries').label('Cumulative Selection x boundaries');
x_boundary_disksel_r.set('contributeto','xboundaries');
x_boundary_disksel_l.set('contributeto','xboundaries');

y_boundary_disksel_t = ucellgeom.feature.create('y_boundary_boxsel_t', 'BoxSelection');
y_boundary_disksel_t.set('entitydim', 2);
y_boundary_disksel_t.set('xmin', 3*a/4-selection_width);
y_boundary_disksel_t.set('xmax', 3*a/4+selection_width);
y_boundary_disksel_t.set('ymin', sqrt(3)*a/2-selection_width/2);
y_boundary_disksel_t.set('ymax', sqrt(3)*a/2+selection_width/2);
y_boundary_disksel_t.set('zmin', th/2-selection_width);
y_boundary_disksel_t.set('zmax', airZth-selection_width);
y_boundary_disksel_t.set('inputent', 'all');
y_boundary_disksel_t.set('condition', 'intersects');

y_boundary_disksel_b = ucellgeom.feature.create('y_boundary_boxsel_b', 'BoxSelection');
y_boundary_disksel_b.set('entitydim', 2);
y_boundary_disksel_b.set('xmin', a/4-selection_width);
y_boundary_disksel_b.set('xmax', a/4+selection_width);
y_boundary_disksel_b.set('ymin', -selection_width/2);
y_boundary_disksel_b.set('ymax', selection_width/2);
y_boundary_disksel_b.set('zmin', th/2-selection_width);
y_boundary_disksel_b.set('zmax', airZth-selection_width);
y_boundary_disksel_b.set('inputent', 'all');
y_boundary_disksel_b.set('condition', 'intersects');

ucellgeom.selection.create('yboundaries','CumulativeSelection');
ucellgeom.selection('yboundaries').label('Cumulative Selection y boundaries');
y_boundary_disksel_t.set('contributeto','yboundaries');
y_boundary_disksel_b.set('contributeto','yboundaries');

P.xEnd1 =  bndindex(ucellgeom, [0 0 0], [sqrt(3)*a/2 -a/2 0]);
P.xEnd2 = bndindex(ucellgeom, [a 0 0], [sqrt(3)*a/2 -a/2 0]);

P.yEnd1 =  bndindex(ucellgeom, [0 0 0], [0 1 0]);
P.yEnd2 = bndindex(ucellgeom, [a (a/2)*sqrt(3) 0], [0 1 0]);
% Note that this will return no indices if there is no boundary at z=0
P.zEnd = bndindex(ucellgeom, [0 0 0], [0 0 1]);
P.zEnd2 = bndindex(ucellgeom, [0 0 airZth], [0 0 1]);
% % debugging
% mphplot(model);
disp(P) % debugging
out = model;
end

function addFillet(P,ucellWP,hole_pos,selection_width)
    w = P.w;
    r = P.r;
    r1 = P.r1;
    r2 = P.r2;
    disksel1_label = sprintf('h_disksel1');
    disksel2_label = sprintf('h_disksel2');
    fillet1_label = sprintf('h_fil1');
    fillet2_label = sprintf('h_fil2');
    % hole 
    disksel1 = ucellWP.geom.feature.create(disksel1_label, 'DiskSelection');
    disksel1.set('entitydim', 0);
    disksel1.set('posx', hole_pos(1));
    disksel1.set('posy', hole_pos(2));
    disksel1.set('r', w);
    disksel1.set('rin', w/(2*sqrt(2)));
    disksel1.set('condition', 'allvertices');
    fil1 = ucellWP.geom.feature.create(fillet1_label, 'Fillet');
    fil1.set('radius', r1);
    fil1.selection('point').named(disksel1_label);
    disksel2 = ucellWP.geom.feature.create(disksel2_label, 'DiskSelection');
    disksel2.set('entitydim', 0);
    disksel2.set('posy', hole_pos(2));
    disksel2.set('posx', hole_pos(1));
    disksel2.set('r', r+selection_width/2);
    disksel2.set('rin', r-selection_width/2);
    disksel2.set('condition', 'allvertices');
    fil2 = ucellWP.geom.feature.create(fillet2_label, 'Fillet');
    fil2.set('radius', r2);
    fil2.selection('point').named(disksel2_label);
end