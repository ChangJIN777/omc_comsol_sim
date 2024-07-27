function [model,P] = buildLowerBoomerangUnitCell_2D(model,P)
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
ucellgeom = model.geom.create(ucellname, 2);
ucellgeom.label(ucelllabel);

%% create unit cell with boomerang geometry in the lower cavity
% ucellWP = ucellgeom.create('wp1', 'WorkPlane');
% ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);

ucellplane = ucellgeom.feature.create('sq1', 'Square');
ucellplane.set('pos', [-a/2 -a/2]);
ucellplane.set('size', a);
rec1 = ucellgeom.feature.create('r1', 'Rectangle');
rec1.set('pos', [a/2-a/2 r/2+a/2+sqrt(3)*w/4-a/2]);
rec1.set('base', 'center');
rec1.set('size', [w r]);
rec2 = ucellgeom.feature.create('r2', 'Rectangle');
rec2.set('pos', [a/2-a/2 a/2-sqrt(3)*w/2+sqrt(3)*w/4-a/2]);
rec2.set('rot', 120);
rec2.set('size', [w r]);
rec3 = ucellgeom.feature.create('r3', 'Rectangle');
rec3.set('pos', [a/2+w/2-a/2 a/2+sqrt(3)*w/4-a/2]);
rec3.set('rot', 240);
rec3.set('size', [w r]);
rec4 = ucellgeom.feature.create('r4', 'Rectangle');
rec4.set('pos', [w/2+sqrt(3)*r/2+a/2-a/2 a/2-r/2+sqrt(3)*w/4-a/2]);
rec4.set('rot', 180);
rec4.set('size', [d h]);
rec5 = ucellgeom.feature.create('r5', 'Rectangle');
rec5.set('pos', [-w/2-sqrt(3)*r/2+a/2+d-a/2 a/2-r/2+sqrt(3)*w/4-a/2]);
rec5.set('rot', 180);
rec5.set('size', [d h]);
pol1 = ucellgeom.feature.create('pol1', 'Polygon');
pol1.set('tableconstr', {'off' 'off'});
pol1.set('source', 'table');
pol1.set('table', [-w/2+a/2-a/2 a/2+sqrt(3)*w/4-a/2; w/2+a/2-a/2 a/2+sqrt(3)*w/4-a/2; a/2-a/2 a/2-sqrt(3)*w/4-a/2; -w/2+a/2-a/2 a/2+sqrt(3)*w/4-a/2]);
compose1 = ucellgeom.feature.create('co1', 'Compose');
compose1.set('formula', 'sq1-r1-r2-r3-r4-r5-pol1');
fillet1 = ucellgeom.feature.create('fil1', 'Fillet');
fillet1.set('radius', r1);
fillet1.selection('point').set('co1(1)', [3 5 8 11 12 14]);
fillet2 = ucellgeom.feature.create('fil2', 'Fillet');
fillet2.set('radius', r2);
fillet2.selection('point').set('fil1(1)', [4 8 9 12 14 17 21]);

holeplane = ucellgeom.feature.create('sq2', 'Square');
holeplane.set('pos', [-a/2 -a/2]);
holeplane.set('size', a);
holeplane.label('Air plane');
ucellgeom.runAll;

%% make selections
mphgeom(model);
P.xEnd1 = 1;
P.xEnd2 = 17;

P.yEnd1 = 2;
P.yEnd2 = 3;

% % debugging
% mphplot(model);
disp(P) % debugging

out = model;