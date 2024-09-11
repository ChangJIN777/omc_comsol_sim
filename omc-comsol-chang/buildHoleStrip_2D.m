function [model,P] = buildHoleStrip_2D(model,P)
%
% buildBoomerangUnitCell.m
%
% Model exported on Jul 21 2024, 16:57 by COMSOL 6.1.0.357.

%% read the input parameters 
a = P.a;        % lattice constant 
b = P.b;        % the width of the hole 
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
ucellgeom = model.geom.create(ucellname, 2);
ucellgeom.label(ucelllabel);

%% create unit cell with boomerang geometry
% ucellWP = ucellgeom.feature.create('wp1', 'WorkPlane');
% ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
% ucellWP.set('unite', true);

ucellplane = ucellgeom.feature.create('r1', 'Rectangle');
ucellplane.label('Base plane');
ucellplane.set('pos', [0 0]);
ucellplane.set('base','center');
ucellplane.set('size',[a sqrt(3)*5*a]);

airHole1 = ucellgeom.feature.create('c1','Circle');
airHole1.set('r',r);
airHole1.set('pos', [0 a*sqrt(3)/2+b]);
airHole1.set('base', 'center');

airHole2 = ucellgeom.feature.create('c2','Circle');
airHole2.set('r',r);
airHole2.set('pos', [0 -(a*sqrt(3)/2+b)]);
airHole2.set('base', 'center');

airHole3 = ucellgeom.feature.create('c3','Circle');
airHole3.set('r',r);
airHole3.set('pos', [0 a*sqrt(3)/2+sqrt(3)*a]);
airHole3.set('base', 'center');

airHole4 = ucellgeom.feature.create('c4','Circle');
airHole4.set('r',r);
airHole4.set('pos', [0 a*sqrt(3)/2+2*sqrt(3)*a]);
airHole4.set('base', 'center');

airHole5 = ucellgeom.feature.create('c5','Circle');
airHole5.set('r',r);
airHole5.set('pos', [0 -(a*sqrt(3)/2+sqrt(3)*a)]);
airHole5.set('base', 'center');

airHole6 = ucellgeom.feature.create('c6','Circle');
airHole6.set('r',r);
airHole6.set('pos', [0 -(a*sqrt(3)/2+2*sqrt(3)*a)]);
airHole6.set('base', 'center');

airHole7 = ucellgeom.feature.create('c7','Circle');
airHole7.set('r',r);
airHole7.set('pos', [a/2 sqrt(3)*a]);
airHole7.set('base', 'center');

airHole8 = ucellgeom.feature.create('c8','Circle');
airHole8.set('r',r);
airHole8.set('pos', [-a/2 sqrt(3)*a]);
airHole8.set('base', 'center');

airHole9 = ucellgeom.feature.create('c9','Circle');
airHole9.set('r',r);
airHole9.set('pos', [a/2 -sqrt(3)*a]);
airHole9.set('base', 'center');

airHole10 = ucellgeom.feature.create('c10','Circle');
airHole10.set('r',r);
airHole10.set('pos', [-a/2 -sqrt(3)*a]);
airHole10.set('base', 'center');

airHole11 = ucellgeom.feature.create('c11','Circle');
airHole11.set('r',r);
airHole11.set('pos', [a/2 sqrt(3)*a+a*sqrt(3)]);
airHole11.set('base', 'center');

airHole12 = ucellgeom.feature.create('c12','Circle');
airHole12.set('r',r);
airHole12.set('pos', [-a/2 sqrt(3)*a+a*sqrt(3)]);
airHole12.set('base', 'center');

airHole13 = ucellgeom.feature.create('c13','Circle');
airHole13.set('r',r);
airHole13.set('pos', [a/2 -sqrt(3)*a-a*sqrt(3)]);
airHole13.set('base', 'center');

airHole14 = ucellgeom.feature.create('c14','Circle');
airHole14.set('r',r);
airHole14.set('pos', [-a/2 -sqrt(3)*a-a*sqrt(3)]);
airHole14.set('base', 'center');

compose1 = ucellgeom.feature.create('co1', 'Compose');
compose1.set('formula', 'r1-c1-c2-c3-c4-c5-c6-c7-c8-c9-c10-c11-c12-c13-c14');

holeplane = ucellgeom.feature.create('r2', 'Rectangle');
holeplane.label('Air plane');
holeplane.set('pos', [0 0]);
holeplane.set('base','center');
holeplane.set('size',[a sqrt(3)*5*a]);

ucellgeom.runAll;

%% Making selections (with box select)
mphgeom(model);
selection_width = 10e-9;

% box selection for the boundary condition 
x_boundary_disksel_r = ucellgeom.feature.create('x_boundary_boxsel_r', 'BoxSelection');
x_boundary_disksel_r.set('entitydim', 1);
x_boundary_disksel_r.set('xmin', a/2-selection_width/2);
x_boundary_disksel_r.set('xmax', a/2+selection_width/2);
x_boundary_disksel_r.set('inputent', 'all');
x_boundary_disksel_r.set('condition', 'inside');

x_boundary_disksel_l = ucellgeom.feature.create('x_boundary_boxsel_l', 'BoxSelection');
x_boundary_disksel_l.set('entitydim', 1);
x_boundary_disksel_l.set('xmin', -a/2-selection_width/2);
x_boundary_disksel_l.set('xmax', -a/2+selection_width/2);
x_boundary_disksel_l.set('inputent', 'all');
x_boundary_disksel_l.set('condition', 'inside');

ucellgeom.selection.create('xboundaries','CumulativeSelection');
ucellgeom.selection('xboundaries').label('Cumulative Selection x boundaries');
x_boundary_disksel_r.set('contributeto','xboundaries');
x_boundary_disksel_l.set('contributeto','xboundaries');

y_boundary_disksel_t = ucellgeom.feature.create('y_boundary_boxsel_t', 'BoxSelection');
y_boundary_disksel_t.set('entitydim', 1);
y_boundary_disksel_t.set('ymin', sqrt(3)*5*a/2-selection_width/2);
y_boundary_disksel_t.set('ymax', sqrt(3)*5*a/2+selection_width/2);
y_boundary_disksel_t.set('inputent', 'all');
y_boundary_disksel_t.set('condition', 'inside');

y_boundary_disksel_b = ucellgeom.feature.create('y_boundary_boxsel_b', 'BoxSelection');
y_boundary_disksel_b.set('entitydim', 1);
y_boundary_disksel_b.set('ymin', -sqrt(3)*5*a/2-selection_width/2);
y_boundary_disksel_b.set('ymax', -sqrt(3)*5*a/2+selection_width/2);
y_boundary_disksel_b.set('inputent', 'all');
y_boundary_disksel_b.set('condition', 'inside');

ucellgeom.selection.create('yboundaries','CumulativeSelection');
ucellgeom.selection('yboundaries').label('Cumulative Selection y boundaries');
y_boundary_disksel_t.set('contributeto','yboundaries');
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