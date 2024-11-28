function [model,P] = buildHoleStrip_3D(model,P)
%
% buildBoomerangUnitCell.m
%
% Model exported on Jul 21 2024, 16:57 by COMSOL 6.1.0.357.

%% read the input parameters 
a = P.a;        % lattice constant 
b = P.b;        % the width of the waveguide (cavity region)
b_wvg = P.b_wvg;    % modulating the wvg width of the unit cell 
th = P.th;        % the height of the hole
r = P.r;        % the height of the hole 
r1 = P.r1;      % the fillet radius of the edges of the hole 
r2 = P.r2;      % the fillet radius of the center of the hole 
abssym = abs(P.mbevenz);    % symmetry in the z direction 
airDiskH = P.airDiskH; % the height the of air disk

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

ucellplane = ucellWP.geom.feature.create('r1', 'Rectangle');
ucellplane.label('Base plane');
ucellplane.set('pos', [-a/2 0]);
ucellplane.set('base','corner');
ucellplane.set('size',[a sqrt(3)*(4+1/2)*a+b_wvg]);

airHole1 = ucellWP.geom.feature.create('c1','Circle');
airHole1.set('r',r);
airHole1.set('pos', [0 a*sqrt(3)/2+b+b_wvg]);
airHole1.set('base', 'center');

airHole2 = ucellWP.geom.feature.create('c2','Circle');
airHole2.set('r',r);
airHole2.set('pos', [0 a*sqrt(3)/2+sqrt(3)*a+b_wvg]);
airHole2.set('base', 'center');

airHole3 = ucellWP.geom.feature.create('c3','Circle');
airHole3.set('r',r);
airHole3.set('pos', [0 a*sqrt(3)/2+2*sqrt(3)*a+b_wvg]);
airHole3.set('base', 'center');

airHole4 = ucellWP.geom.feature.create('c4','Circle');
airHole4.set('r',r);
airHole4.set('pos', [a/2 sqrt(3)*a+b_wvg]);
airHole4.set('base', 'center');

airHole5 = ucellWP.geom.feature.create('c5','Circle');
airHole5.set('r',r);
airHole5.set('pos', [-a/2 sqrt(3)*a+b_wvg]);
airHole5.set('base', 'center');

airHole6 = ucellWP.geom.feature.create('c6','Circle');
airHole6.set('r',r);
airHole6.set('pos', [a/2 sqrt(3)*a+a*sqrt(3)+b_wvg]);
airHole6.set('base', 'center');

airHole7 = ucellWP.geom.feature.create('c7','Circle');
airHole7.set('r',r);
airHole7.set('pos', [-a/2 sqrt(3)*a+a*sqrt(3)+b_wvg]);
airHole7.set('base', 'center');

airHole8 = ucellWP.geom.feature.create('c8','Circle');
airHole8.set('r',r);
airHole8.set('pos', [a/2 sqrt(3)*a+a*sqrt(3)+a*sqrt(3)+b_wvg]);
airHole8.set('base', 'center');

airHole9 = ucellWP.geom.feature.create('c9','Circle');
airHole9.set('r',r);
airHole9.set('pos', [-a/2 sqrt(3)*a+a*sqrt(3)+a*sqrt(3)+b_wvg]);
airHole9.set('base', 'center');

airHole10 = ucellWP.geom.feature.create('c10','Circle');
airHole10.set('r',r);
airHole10.set('pos', [0 sqrt(3)*a+a*sqrt(3)+a*sqrt(3)+a*sqrt(3)/2+b_wvg]);
airHole10.set('base', 'center');

airHole11 = ucellWP.geom.feature.create('c11','Circle');
airHole11.set('r',r);
airHole11.set('pos', [-a/2 sqrt(3)*a+a*sqrt(3)+a*sqrt(3)+a*sqrt(3)+b_wvg]);
airHole11.set('base', 'center');

airHole12 = ucellWP.geom.feature.create('c12','Circle');
airHole12.set('r',r);
airHole12.set('pos', [a/2 sqrt(3)*a+a*sqrt(3)+a*sqrt(3)+a*sqrt(3)+b_wvg]);
airHole12.set('base', 'center');

airHole13 = ucellWP.geom.feature.create('c13','Circle');
airHole13.set('r',r);
airHole13.set('pos', [0 sqrt(3)*a+a*sqrt(3)+a*sqrt(3)+a*sqrt(3)+a*sqrt(3)/2+b_wvg]);
airHole13.set('base', 'center');

compose1 = ucellWP.geom.feature.create('co1', 'Compose');
compose1.set('formula', 'r1-c1-c2-c3-c4-c5-c6-c7-c8-c9-c10-c11-c12-c13');

extrude = ucellgeom.feature.create('ext1', 'Extrude');
extrude.setIndex('distance', th, 0);
extrude.selection('input').set({'wp1'});

% add the air disk 
airPlane = ucellgeom.feature.create('wpair', 'WorkPlane');
airPlane.set('quickplane', 'zx');
airPlane.set('unite', true);
diskBase = airPlane.geom.feature.create('r1', 'Rectangle');
diskBase.set('pos', [0 -a/2]);
diskBase.set('size', [airDiskH a]);
revolve = ucellgeom.feature.create('rev1', 'Revolve');
revolve.set('angle2', -90);
revolve.selection('input').set({'wpair'});

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
    symZPlane.set('table', [-a/2 0; a/2 0; a/2 sqrt(3)*(4+1/2)*a+b_wvg; -a/2 sqrt(3)*(4+1/2)*a+b_wvg; -a/2 0]);

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
    ZsymSel.set('ymin', -delta).set('ymax', P.airDiskH+delta);
    ZsymSel.set('zmin', -delta).set('zmax', delta);
    ZsymSel.set('entitydim', 2).set('condition', 'allvertices'); % only want beam surface, exclude air holes
    ucellgeom.runCurrent;
    inds = model.selection([ucellname,'_ZsymSel']).inputEntities();
    P.bndSel.Zsym = inds';

end 


% holeplane = ucellWP.geom.feature.create('r2', 'Rectangle');
% holeplane.label('Air plane');
% holeplane.set('pos', [0 0]);
% holeplane.set('base','center');
% holeplane.set('size',[a sqrt(3)*5*a+2*b_wvg]);

% 
% PML_length = 1e-6;
% top_PML = ucellgeom.feature.create('r_PML_top', 'Rectangle');
% top_PML.label('PML_plane_top');
% top_PML.set('pos', [0 (PML_length+sqrt(3)*5*a)/2]);
% top_PML.set('base','center');
% top_PML.set('size',[a PML_length]);
% 
% bottom_PML = ucellgeom.feature.create('r_PML_bottom', 'Rectangle');
% bottom_PML.label('PML_plane_bottom');
% bottom_PML.set('pos', [0 -(PML_length+sqrt(3)*5*a)/2]);
% bottom_PML.set('base','center');
% bottom_PML.set('size',[a PML_length]);

ucellgeom.runAll;
% 
% %% add the PML regions
% model.component('comp1').coordSystem.create('pml1', 'PML');
% model.component('comp1').geom('geom1').run;
% model.component('comp1').coordSystem('pml1').selection.set([1 7]);

% %% Making selections (with box select)
% mphgeom(model);
selection_width = 10e-9;

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

y_boundary_symmetry = ucellgeom.feature.create('y_boundary_symmetry', 'BoxSelection');
y_boundary_symmetry.set('entitydim', 2);
y_boundary_symmetry.set('ymin', -selection_width/2);
y_boundary_symmetry.set('ymax', selection_width/2);
y_boundary_symmetry.set('inputent', 'all');
y_boundary_symmetry.set('condition', 'inside');
% 
% y_boundary_disksel_b = ucellgeom.feature.create('y_boundary_boxsel_b', 'BoxSelection');
% y_boundary_disksel_b.set('entitydim', 1);
% y_boundary_disksel_b.set('ymin', -sqrt(3)*5*a/2-b_wvg-selection_width/2);
% y_boundary_disksel_b.set('ymax', -sqrt(3)*5*a/2-b_wvg+selection_width/2);
% y_boundary_disksel_b.set('inputent', 'all');
% y_boundary_disksel_b.set('condition', 'inside');
% 
ucellgeom.selection.create('yboundaries','CumulativeSelection');
ucellgeom.selection('yboundaries').label('Cumulative Selection y boundaries');
y_boundary_symmetry.set('contributeto','yboundaries');

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