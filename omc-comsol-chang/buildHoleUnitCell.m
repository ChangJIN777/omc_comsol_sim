function [model,P] = buildHoleUnitCell(model,P)
%
% buildBoomerangUnitCell.m
%
% Model exported on Jul 21 2024, 16:57 by COMSOL 6.1.0.357.

%% read the input parameters 
a = P.a;        % lattice constant 
hx = P.hx;        % the diameter of the hole in x 
hy = P.hy;        % the diameter of the hole in y 
th = P.th;        % the height of the hole
beam_width = P.beam_width; % the beam width 
d_in = P.d_in; % the sidewall angle for the inside
d_out = P.d_out; % the sidewall angle for the outside
abssym = abs(P.mbevenz);    % symmetry in the z direction 
airDiskH = P.airDiskH; % the height the of air disk

%% Create component 
ucellcomp = model.modelNode.create('comp1');
ucellcomp.label('Unit cell FEM simulation');
ucelllabel = 'Unit Cell';
ucellname = 'geom1';
ucellgeom = model.geom.create(ucellname, 3);
ucellgeom.label(ucelllabel);

%% create unit cell
ucellWP = ucellgeom.feature.create('wp1', 'WorkPlane');
ucellWP.set('quickplane', 'yz');
ucellWP.set('quickx', -a/2);
ucellWP.set('unite', true);
basePolygon = ucellWP.geom.feature.create('pol1', 'Polygon');
basePolygon.set('source', 'table');
basePolygon.set('table', [-beam_width/2 -th/2; beam_width/2 -th/2; beam_width/2-tan(d_out)*th th/2 ; -beam_width/2+tan(d_out)*th th/2]);
ext1 = ucellgeom.feature.create('ext1', 'Extrude');
ext1.setIndex('distance', a, 0);
ext1.selection('input').set({'wp1'});
econ1 = ucellgeom.feature.create('econ1', 'ECone');
econ1.set('pos', [0 0 -th/2]);
econ1.set('axis', [0 0]);
hx_top = hx*(1-th*tan(d_in)/hy); hy_top = hy*(1-th*tan(d_in)/hy); 
econ1.set('semiaxes', [hy_top/2 hx_top/2]);
econ1.set('h', th);
econ1.set('rat', 1/(1-th*tan(d_in)/hy));
compose = ucellgeom.feature.create('co1', 'Compose');
compose.set('formula', 'ext1-econ1');
ucellgeom.run;
ucellgeom.run('fin');
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

ucellgeom.runAll;

model.component('comp1').geom('geom1').create('wp7', 'WorkPlane');
if ~P.run_optical
    model.component('comp1').geom('geom1').feature('wp7').active(false);
end
model.component('comp1').geom('geom1').feature('wp7').set('quickplane', 'xz');
model.component('comp1').geom('geom1').feature('wp7').set('unite', true);
model.component('comp1').geom('geom1').feature('wp7').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp7').geom.feature('r1').set('size', [a th]);
model.component('comp1').geom('geom1').feature('wp7').geom.feature('r1').set('pos', [-a/2 -th/2]);
model.component('comp1').geom('geom1').create('ext8', 'Extrude');
if ~P.run_optical
    model.component('comp1').geom('geom1').feature('ext8').active(false);
end
model.component('comp1').geom('geom1').feature('ext8').setIndex('distance', beam_width/2, 0);
model.component('comp1').geom('geom1').feature('ext8').selection('input').set({'wp7'});
model.component('comp1').geom('geom1').create('co2', 'Compose');
if ~P.run_optical
    model.component('comp1').geom('geom1').feature('co2').active(false);
end
model.component('comp1').geom('geom1').feature('co2').set('formula', 'co1-ext8');
model.component('comp1').geom('geom1').create('wp8', 'WorkPlane');
if ~P.run_optical
    model.component('comp1').geom('geom1').feature('wp8').active(false);
end
model.component('comp1').geom('geom1').feature('wp8').set('quickplane', 'xz');
model.component('comp1').geom('geom1').feature('wp8').set('unite', true);
model.component('comp1').geom('geom1').feature('wp8').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp8').geom.feature('r1').set('size', [a airDiskH]);
model.component('comp1').geom('geom1').feature('wp8').geom.feature('r1').set('pos', [-a/2 0]);
model.component('comp1').geom('geom1').create('rev1', 'Revolve');
if ~P.run_optical
    model.component('comp1').geom('geom1').feature('rev1').active(false);
end
model.component('comp1').geom('geom1').feature('rev1').set('angle2', 180);
model.component('comp1').geom('geom1').feature('rev1').set('axis', [-1 0]);
model.component('comp1').geom('geom1').feature('rev1').selection('input').set({'wp8'});
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');


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

if P.run_optical    
    y_boundary_symmetry = ucellgeom.feature.create('y_boundary_symmetry', 'BoxSelection');
    y_boundary_symmetry.set('entitydim', 2);
    y_boundary_symmetry.set('ymin', -selection_width/2);
    y_boundary_symmetry.set('ymax', selection_width/2);
    y_boundary_symmetry.set('inputent', 'all');
    y_boundary_symmetry.set('condition', 'inside');
    ucellgeom.selection.create('yboundaries','CumulativeSelection');
    ucellgeom.selection('yboundaries').label('Cumulative Selection y boundaries');
    y_boundary_symmetry.set('contributeto','yboundaries');
end

model.component('comp1').geom('geom1').run;

% Note that this will return no indices if there is no boundary at z=0

% % debugging
% mphplot(model);
disp(P) % debugging
out = model;
end