function [model,P] = buildHoleUnitCell_2D(model,P)
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

ucellplane = ucellgeom.feature.create('pol1', 'Polygon');
ucellplane.label('Base plane');
ucellplane.set('source', 'table');
ucellplane.set('table', [0 0; a/2 (a/2)*sqrt(3); a*(3/2) (a/2)*sqrt(3); a 0; 0 0]);

airHole = ucellgeom.feature.create('c1','Circle');
airHole.set('r',r);
airHole.set('pos', [a*(1/2+1/4) a*sqrt(3)/4+b]);
airHole.set('base', 'center');

holeplane = ucellgeom.feature.create('pol3', 'Polygon');
holeplane.label('Air plane');
holeplane.set('source', 'table');
holeplane.set('table', [0 0; a/2 (a/2)*sqrt(3); a*(3/2) (a/2)*sqrt(3); a 0; 0 0]);
ucellgeom.runAll;

%% Making selections (with box select)
mphgeom(model);
P.xEnd1 = 1;
P.xEnd2 = 4;

P.yEnd1 = 2;
P.yEnd2 = 3;
% Note that this will return no indices if there is no boundary at z=0

% % debugging
% mphplot(model);
disp(P) % debugging
out = model;
end