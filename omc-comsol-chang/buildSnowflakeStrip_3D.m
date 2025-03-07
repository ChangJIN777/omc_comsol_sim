function [model,P] = buildSnowflakeStrip_3D(model,P)
%
% buildBoomerangUnitCell.m
%
% Model exported on Jul 21 2024, 16:57 by COMSOL 6.1.0.357.

%% read the input parameters 
a = P.a;        % lattice constant 
b = P.b;        % the offset for the optical mode
b_base = P.b_base; 
th = P.th;        % the height of the hole
w = P.w;        % the width of the hole
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
ucellWP.set('unite', true);

% build the base plate for the strip
ucellplane = ucellWP.geom.feature.create('r_base', 'Rectangle');
ucellplane.label('Base plane');
ucellplane.set('pos', [-a/2 0]);
ucellplane.set('base','corner');
ucellplane.set('size',[a sqrt(3)*(2+1/2)*a+b_base+b]);

% add the holes for the strip
% hole #1 
c1_rec_1 = ucellWP.geom.feature.create('c1_r1', 'Rectangle');
c1_rec_1.set('pos', [0 a*sqrt(3)/2+b+b_base]);
c1_rec_1.set('base', 'center');
c1_rec_1.set('size', [2*r w]);

c1_rec_2 = ucellWP.geom.feature.create('c1_r2', 'Rectangle');
c1_rec_2.set('pos', [0 a*sqrt(3)/2+b+b_base]);
c1_rec_2.set('rot', 60);
c1_rec_2.set('base', 'center');
c1_rec_2.set('size', [2*r w]);

c1_rec_3 = ucellWP.geom.feature.create('c1_r3', 'Rectangle');
c1_rec_3.set('pos', [0 a*sqrt(3)/2+b+b_base]);
c1_rec_3.set('rot', 120);
c1_rec_3.set('base', 'center');
c1_rec_3.set('size', [2*r w]);

% hole #2 
c2_rec_1 = ucellWP.geom.feature.create('c2_r1', 'Rectangle');
c2_rec_1.set('pos', [0 -(a*sqrt(3)/2+b+b_base)]);
c2_rec_1.set('base', 'center');
c2_rec_1.set('size', [2*r w]);

c2_rec_2 = ucellWP.geom.feature.create('c2_r2', 'Rectangle');
c2_rec_2.set('pos', [0 -(a*sqrt(3)/2+b+b_base)]);
c2_rec_2.set('rot', 60);
c2_rec_2.set('base', 'center');
c2_rec_2.set('size', [2*r w]);

c2_rec_3 = ucellWP.geom.feature.create('c2_r3', 'Rectangle');
c2_rec_3.set('pos', [0 -(a*sqrt(3)/2+b+b_base)]);
c2_rec_3.set('rot', 120);
c2_rec_3.set('base', 'center');
c2_rec_3.set('size', [2*r w]);

% hole #3 
c3_rec_1 = ucellWP.geom.feature.create('c3_r1', 'Rectangle');
c3_rec_1.set('pos', [0 a*sqrt(3)/2+sqrt(3)*a+b_base]);
c3_rec_1.set('base', 'center');
c3_rec_1.set('size', [2*r w]);

c3_rec_2 = ucellWP.geom.feature.create('c3_r2', 'Rectangle');
c3_rec_2.set('pos', [0 a*sqrt(3)/2+sqrt(3)*a+b_base]);
c3_rec_2.set('rot', 60);
c3_rec_2.set('base', 'center');
c3_rec_2.set('size', [2*r w]);

c3_rec_3 = ucellWP.geom.feature.create('c3_r3', 'Rectangle');
c3_rec_3.set('pos', [0 a*sqrt(3)/2+sqrt(3)*a+b_base]);
c3_rec_3.set('rot', 120);
c3_rec_3.set('base', 'center');
c3_rec_3.set('size', [2*r w]);

% hole #4 
c4_rec_1 = ucellWP.geom.feature.create('c4_r1', 'Rectangle');
c4_rec_1.set('pos', [0 a*sqrt(3)/2+2*sqrt(3)*a+b_base]);
c4_rec_1.set('base', 'center');
c4_rec_1.set('size', [2*r w]);

c4_rec_2 = ucellWP.geom.feature.create('c4_r2', 'Rectangle');
c4_rec_2.set('pos', [0 a*sqrt(3)/2+2*sqrt(3)*a+b_base]);
c4_rec_2.set('rot', 60);
c4_rec_2.set('base', 'center');
c4_rec_2.set('size', [2*r w]);

c4_rec_3 = ucellWP.geom.feature.create('c4_r3', 'Rectangle');
c4_rec_3.set('pos', [0 a*sqrt(3)/2+2*sqrt(3)*a+b_base]);
c4_rec_3.set('rot', 120);
c4_rec_3.set('base', 'center');
c4_rec_3.set('size', [2*r w]);

% hole #5
c5_rec_1 = ucellWP.geom.feature.create('c5_r1', 'Rectangle');
c5_rec_1.set('pos', [0 -(a*sqrt(3)/2+sqrt(3)*a+b_base)]);
c5_rec_1.set('base', 'center');
c5_rec_1.set('size', [2*r w]);

c5_rec_2 = ucellWP.geom.feature.create('c5_r2', 'Rectangle');
c5_rec_2.set('pos', [0 -(a*sqrt(3)/2+sqrt(3)*a+b_base)]);
c5_rec_2.set('rot', 60);
c5_rec_2.set('base', 'center');
c5_rec_2.set('size', [2*r w]);

c5_rec_3 = ucellWP.geom.feature.create('c5_r3', 'Rectangle');
c5_rec_3.set('pos', [0 -(a*sqrt(3)/2+sqrt(3)*a+b_base)]);
c5_rec_3.set('rot', 120);
c5_rec_3.set('base', 'center');
c5_rec_3.set('size', [2*r w]);

% hole #6
c6_rec_1 = ucellWP.geom.feature.create('c6_r1', 'Rectangle');
c6_rec_1.set('pos', [0 -(a*sqrt(3)/2+2*sqrt(3)*a+b_base)]);
c6_rec_1.set('base', 'center');
c6_rec_1.set('size', [2*r w]);

c6_rec_2 = ucellWP.geom.feature.create('c6_r2', 'Rectangle');
c6_rec_2.set('pos', [0 -(a*sqrt(3)/2+2*sqrt(3)*a+b_base)]);
c6_rec_2.set('rot', 60);
c6_rec_2.set('base', 'center');
c6_rec_2.set('size', [2*r w]);

c6_rec_3 = ucellWP.geom.feature.create('c6_r3', 'Rectangle');
c6_rec_3.set('pos', [0 -(a*sqrt(3)/2+2*sqrt(3)*a+b_base)]);
c6_rec_3.set('rot', 120);
c6_rec_3.set('base', 'center');
c6_rec_3.set('size', [2*r w]);

% hole #7
c7_rec_1 = ucellWP.geom.feature.create('c7_r1', 'Rectangle');
c7_rec_1.set('pos', [a/2 sqrt(3)*a+b_base]);
c7_rec_1.set('base', 'center');
c7_rec_1.set('size', [2*r w]);

c7_rec_2 = ucellWP.geom.feature.create('c7_r2', 'Rectangle');
c7_rec_2.set('pos', [a/2 sqrt(3)*a+b_base]);
c7_rec_2.set('rot', 60);
c7_rec_2.set('base', 'center');
c7_rec_2.set('size', [2*r w]);

c7_rec_3 = ucellWP.geom.feature.create('c7_r3', 'Rectangle');
c7_rec_3.set('pos', [a/2 sqrt(3)*a+b_base]);
c7_rec_3.set('rot', 120);
c7_rec_3.set('base', 'center');
c7_rec_3.set('size', [2*r w]);

% hole #8
c8_rec_1 = ucellWP.geom.feature.create('c8_r1', 'Rectangle');
c8_rec_1.set('pos', [-a/2 sqrt(3)*a+b_base]);
c8_rec_1.set('base', 'center');
c8_rec_1.set('size', [2*r w]);

c8_rec_2 = ucellWP.geom.feature.create('c8_r2', 'Rectangle');
c8_rec_2.set('pos', [-a/2 sqrt(3)*a+b_base]);
c8_rec_2.set('rot', 60);
c8_rec_2.set('base', 'center');
c8_rec_2.set('size', [2*r w]);

c8_rec_3 = ucellWP.geom.feature.create('c8_r3', 'Rectangle');
c8_rec_3.set('pos', [-a/2 sqrt(3)*a+b_base]);
c8_rec_3.set('rot', 120);
c8_rec_3.set('base', 'center');
c8_rec_3.set('size', [2*r w]);

% hole #9
c9_rec_1 = ucellWP.geom.feature.create('c9_r1', 'Rectangle');
c9_rec_1.set('pos', [a/2 -sqrt(3)*a-b_base]);
c9_rec_1.set('base', 'center');
c9_rec_1.set('size', [2*r w]);

c9_rec_2 = ucellWP.geom.feature.create('c9_r2', 'Rectangle');
c9_rec_2.set('pos', [a/2 -sqrt(3)*a-b_base]);
c9_rec_2.set('rot', 60);
c9_rec_2.set('base', 'center');
c9_rec_2.set('size', [2*r w]);

c9_rec_3 = ucellWP.geom.feature.create('c9_r3', 'Rectangle');
c9_rec_3.set('pos', [a/2 -sqrt(3)*a-b_base]);
c9_rec_3.set('rot', 120);
c9_rec_3.set('base', 'center');
c9_rec_3.set('size', [2*r w]);

% hole #10
c10_rec_1 = ucellWP.geom.feature.create('c10_r1', 'Rectangle');
c10_rec_1.set('pos', [-a/2 -sqrt(3)*a-b_base]);
c10_rec_1.set('base', 'center');
c10_rec_1.set('size', [2*r w]);

c10_rec_2 = ucellWP.geom.feature.create('c10_r2', 'Rectangle');
c10_rec_2.set('pos', [-a/2 -sqrt(3)*a-b_base]);
c10_rec_2.set('rot', 60);
c10_rec_2.set('base', 'center');
c10_rec_2.set('size', [2*r w]);

c10_rec_3 = ucellWP.geom.feature.create('c10_r3', 'Rectangle');
c10_rec_3.set('pos', [-a/2 -sqrt(3)*a-b_base]);
c10_rec_3.set('rot', 120);
c10_rec_3.set('base', 'center');
c10_rec_3.set('size', [2*r w]);

% hole #11
c11_rec_1 = ucellWP.geom.feature.create('c11_r1', 'Rectangle');
c11_rec_1.set('pos', [a/2 sqrt(3)*a+a*sqrt(3)+b_base]);
c11_rec_1.set('base', 'center');
c11_rec_1.set('size', [2*r w]);

c11_rec_2 = ucellWP.geom.feature.create('c11_r2', 'Rectangle');
c11_rec_2.set('pos', [a/2 sqrt(3)*a+a*sqrt(3)+b_base]);
c11_rec_2.set('rot', 60);
c11_rec_2.set('base', 'center');
c11_rec_2.set('size', [2*r w]);

c11_rec_3 = ucellWP.geom.feature.create('c11_r3', 'Rectangle');
c11_rec_3.set('pos', [a/2 sqrt(3)*a+a*sqrt(3)+b_base]);
c11_rec_3.set('rot', 120);
c11_rec_3.set('base', 'center');
c11_rec_3.set('size', [2*r w]);

% hole #12
c12_rec_1 = ucellWP.geom.feature.create('c12_r1', 'Rectangle');
c12_rec_1.set('pos', [-a/2 sqrt(3)*a+a*sqrt(3)+b_base]);
c12_rec_1.set('base', 'center');
c12_rec_1.set('size', [2*r w]);

c12_rec_2 = ucellWP.geom.feature.create('c12_r2', 'Rectangle');
c12_rec_2.set('pos', [-a/2 sqrt(3)*a+a*sqrt(3)+b_base]);
c12_rec_2.set('rot', 60);
c12_rec_2.set('base', 'center');
c12_rec_2.set('size', [2*r w]);

c12_rec_3 = ucellWP.geom.feature.create('c12_r3', 'Rectangle');
c12_rec_3.set('pos', [-a/2 sqrt(3)*a+a*sqrt(3)+b_base]);
c12_rec_3.set('rot', 120);
c12_rec_3.set('base', 'center');
c12_rec_3.set('size', [2*r w]);

% hole #13
c13_rec_1 = ucellWP.geom.feature.create('c13_r1', 'Rectangle');
c13_rec_1.set('pos', [a/2 -sqrt(3)*a-a*sqrt(3)-b_base]);
c13_rec_1.set('base', 'center');
c13_rec_1.set('size', [2*r w]);

c13_rec_2 = ucellWP.geom.feature.create('c13_r2', 'Rectangle');
c13_rec_2.set('pos', [a/2 -sqrt(3)*a-a*sqrt(3)-b_base]);
c13_rec_2.set('rot', 60);
c13_rec_2.set('base', 'center');
c13_rec_2.set('size', [2*r w]);

c13_rec_3 = ucellWP.geom.feature.create('c13_r3', 'Rectangle');
c13_rec_3.set('pos', [a/2 -sqrt(3)*a-a*sqrt(3)-b_base]);
c13_rec_3.set('rot', 120);
c13_rec_3.set('base', 'center');
c13_rec_3.set('size', [2*r w]);

% hole #14
c14_rec_1 = ucellWP.geom.feature.create('c14_r1', 'Rectangle');
c14_rec_1.set('pos', [-a/2 -sqrt(3)*a-a*sqrt(3)-b_base]);
c14_rec_1.set('base', 'center');
c14_rec_1.set('size', [2*r w]);

c14_rec_2 = ucellWP.geom.feature.create('c14_r2', 'Rectangle');
c14_rec_2.set('pos', [-a/2 -sqrt(3)*a-a*sqrt(3)-b_base]);
c14_rec_2.set('rot', 60);
c14_rec_2.set('base', 'center');
c14_rec_2.set('size', [2*r w]);

c14_rec_3 = ucellWP.geom.feature.create('c14_r3', 'Rectangle');
c14_rec_3.set('pos', [-a/2 -sqrt(3)*a-a*sqrt(3)-b_base]);
c14_rec_3.set('rot', 120);
c14_rec_3.set('base', 'center');
c14_rec_3.set('size', [2*r w]);

composit_geom = ucellWP.geom.feature.create('co1', 'Compose');
composit_geom.set('formula', 'r_base-(c1_r1 + c1_r2 + c1_r3 + c2_r1 + c2_r2 + c2_r3 + c3_r1 + c3_r2 + c3_r3 + c4_r1 + c4_r2 + c4_r3 + c5_r1 + c5_r2 + c5_r3 + c6_r1 + c6_r2 + c6_r3 + c7_r1 + c7_r2 + c7_r3 + c8_r1 + c8_r2 + c8_r3 + c9_r1 + c9_r2 + c9_r3 + c10_r1 + c10_r2 + c10_r3 + c11_r1 + c11_r2 + c11_r3 + c12_r1 + c12_r2 + c12_r3 + c13_r1 + c13_r2 + c13_r3 + c14_r1 + c14_r2 + c14_r3)');

% the position of all the holes 
hole_pos = [0 a*sqrt(3)/2+b+b_base; 0 -(a*sqrt(3)/2+b+b_base); 0 a*sqrt(3)/2+sqrt(3)*a+b_base;
    0 a*sqrt(3)/2+2*sqrt(3)*a+b_base; 0 -(a*sqrt(3)/2+sqrt(3)*a+b_base);0 -(a*sqrt(3)/2+2*sqrt(3)*a+b_base);
    a/2 sqrt(3)*a+b_base; -a/2 sqrt(3)*a+b_base; a/2 -sqrt(3)*a-b_base; -a/2 -sqrt(3)*a-b_base;
    a/2 sqrt(3)*a+a*sqrt(3)+b_base; -a/2 sqrt(3)*a+a*sqrt(3)+b_base; a/2 -sqrt(3)*a-a*sqrt(3)-b_base;
    -a/2 -sqrt(3)*a-a*sqrt(3)-b_base];

selection_width = 5e-9;
for i=1:14
    addFillet(P,ucellWP.geom,i,hole_pos(i,:),selection_width);
end

extrude = ucellgeom.feature.create('ext1', 'Extrude');
% extrude.set('extrudefrom', 'faces');
extrude.setIndex('distance', th, 0);
extrude.selection('input').set({'wp1'});

% holeplane = ucellgeom.feature.create('r_air', 'Rectangle');
% holeplane.label('Air plane');
% holeplane.set('pos', [0 0]);
% holeplane.set('base','center');
% holeplane.set('size',[a sqrt(3)*5*a+b*2+b_base*2]);

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
    symZPlane.set('table', [-symW/2 0; -symW/2 sqrt(3)*(2+1/2)*a+b_base+b; symW/2 sqrt(3)*(2+1/2)*a+b_base+b; symW/2 0; -symW/2 0]);

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
    ZsymSel.set('ymin', -delta).set('ymax', sqrt(3)*(2+1/2)*a+b_base+b+delta);
    ZsymSel.set('zmin', -delta).set('zmax', delta);
    ZsymSel.set('entitydim', 2).set('condition', 'allvertices'); % only want beam surface, exclude air holes
    ucellgeom.runCurrent;
    inds = model.selection([ucellname,'_ZsymSel']).inputEntities();
    P.bndSel.Zsym = inds';

end

%% Making selections (with box select)
mphgeom(model);

% box selection for the boundary condition 
x_boundary_disksel_r = ucellgeom.feature.create('x_boundary_boxsel_r', 'BoxSelection');
x_boundary_disksel_r.set('entitydim', 2);
x_boundary_disksel_r.set('xmin', a/2-selection_width/2);
x_boundary_disksel_r.set('xmax', a/2+selection_width/2);
x_boundary_disksel_r.set('inputent', 'all');
x_boundary_disksel_r.set('condition', 'inside');

x_boundary_disksel_l = ucellgeom.feature.create('x_boundary_boxsel_l', 'BoxSelection');
x_boundary_disksel_l.set('entitydim', 2);
x_boundary_disksel_l.set('xmin', -a/2-selection_width/2);
x_boundary_disksel_l.set('xmax', -a/2+selection_width/2);
x_boundary_disksel_l.set('inputent', 'all');
x_boundary_disksel_l.set('condition', 'inside');

ucellgeom.selection.create('xboundaries','CumulativeSelection');
ucellgeom.selection('xboundaries').label('Cumulative Selection x boundaries');
x_boundary_disksel_r.set('contributeto','xboundaries');
x_boundary_disksel_l.set('contributeto','xboundaries');

% y_boundary_disksel_t = ucellgeom.feature.create('y_boundary_boxsel_t', 'BoxSelection');
% y_boundary_disksel_t.set('entitydim', 2);
% y_boundary_disksel_t.set('ymin', sqrt(3)*5*a/2+b+b_base-selection_width/2);
% y_boundary_disksel_t.set('ymax', sqrt(3)*5*a/2+b+b_base+selection_width/2);
% y_boundary_disksel_t.set('inputent', 'all');
% y_boundary_disksel_t.set('condition', 'inside');

y_boundary_disksel_b = ucellgeom.feature.create('y_boundary_boxsel_b', 'BoxSelection');
y_boundary_disksel_b.set('entitydim', 2);
y_boundary_disksel_b.set('ymin', -selection_width/2);
y_boundary_disksel_b.set('ymax', selection_width/2);
y_boundary_disksel_b.set('inputent', 'all');
y_boundary_disksel_b.set('condition', 'inside');

ucellgeom.selection.create('yboundaries','CumulativeSelection');
ucellgeom.selection('yboundaries').label('Cumulative Selection y boundaries');
% y_boundary_disksel_t.set('contributeto','yboundaries');
y_boundary_disksel_b.set('contributeto','yboundaries');

P.xEnd1 = [1 3 4 5 6 7 8 9 10];
P.xEnd2 = [16 17 18 19 20 21 22 23 24];

P.yEnd1 = [2 12 14];
P.yEnd2 = [11 13 15];
% Note that this will return no indices if there is no boundary at z=0

% % debugging
% mphplot(model);
disp(P) % debugging
out = model;
end

function addFillet(P,ucellgeom,hole_indx,hole_pos,selection_width)
    w = P.w;
    r = P.r;
    r1 = P.r1;
    r2 = P.r2;
    disksel1_label = sprintf('h%d_disksel1',hole_indx);
    disksel2_label = sprintf('h%d_disksel2',hole_indx);
    fillet1_label = sprintf('h%d_fil1',hole_indx);
    fillet2_label = sprintf('h%d_fil2',hole_indx);
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