function [model,P] = buildBoomerangUnitCell_2D(model,P)
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
ucellgeom = model.geom.create(ucellname, 2);
ucellgeom.label(ucelllabel);

%% create unit cell with boomerang geometry
% ucellWP = ucellgeom.feature.create('wp1', 'WorkPlane');
% ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
% ucellWP.set('unite', true);

ucellplane = ucellgeom.feature.create('pol1', 'Polygon');
ucellplane.label('Base plane');
ucellplane.set('source', 'table');
ucellplane.set('table', [0 0; a/2 (a/2)*sqrt(3); a*(3/2) (a/2)*sqrt(3); a 0; 0 0]);
rec_1 = ucellgeom.feature.create('r1', 'Rectangle');
rec_1.set('pos', [a/2+w/2 r/2+(w/4)*sqrt(3)+(a/2)/sqrt(3)]);
rec_1.set('base', 'center');
rec_1.set('size', [w r]);
rec_2 = ucellgeom.feature.create('r2', 'Rectangle');
rec_2.set('pos', [-a*(3/4)+a/2+a*(3/4)-a/2+w/2+a/2 -(w/4)*sqrt(3)+(a/2)/sqrt(3)]);
rec_2.set('rot', 120);
rec_2.set('size', [w r]);
rec_3 = ucellgeom.feature.create('r3', 'Rectangle');
rec_3.set('pos', [-a*(3/4)+a/2+w/2+a*(3/4)-a/2+w/2+a/2 (w/4)*sqrt(3)+(a/2)/sqrt(3)]);
rec_3.set('rot', 240);
rec_3.set('size', [w r]);

ucellgeom.create('fil1', 'Fillet');
ucellgeom.feature('fil1').set('radius', r1);
ucellgeom.feature('fil1').selection('point').set('r1(1)', [3 4]);
ucellgeom.feature('fil1').selection('point').set('r3(1)', [3 4]);
ucellgeom.feature('fil1').selection('point').set('r2(1)', [3 4]);
centerTriangle = ucellgeom.feature.create('pol2', 'Polygon');
centerTriangle.set('source', 'table');
centerTriangle.set('table', [a/2 (w/2)*sqrt(3)/2+(a/2)/sqrt(3); w+a/2 (w/2)*sqrt(3)/2+(a/2)/sqrt(3); a/2+w/2 -(w/2)*sqrt(3)/2+(a/2)/sqrt(3); a/2 (w/2)*sqrt(3)/2+(a/2)/sqrt(3)]);
composit_geom = ucellgeom.feature.create('co1', 'Compose');
composit_geom.set('formula', 'pol1-fil1(1)-fil1(2)-fil1(3)-pol2');
ucellgeom.create('fil2', 'Fillet');
ucellgeom.feature('fil2').set('radius', r2);
ucellgeom.feature('fil2').selection('point').set('co1(1)', [6 10 12]);
ucellgeom.runAll;

holeplane = ucellgeom.feature.create('pol3', 'Polygon');
holeplane.label('Air plane');
holeplane.set('source', 'table');
holeplane.set('table', [0 0; a/2 (a/2)*sqrt(3); a*(3/2) (a/2)*sqrt(3); a 0; 0 0]);
ucellgeom.runAll;

%% Making selections (manual)
mphgeom(model);
P.xEnd1 = 1;
P.xEnd2 = 12;

P.yEnd1 = 2;
P.yEnd2 = 7;
% Note that this will return no indices if there is no boundary at z=0

% % debugging
% mphplot(model);
disp(P) % debugging
out = model;
end