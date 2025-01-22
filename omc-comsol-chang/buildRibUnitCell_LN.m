function [model,P] = buildRibUnitCell_LN(model,P)

%% read the input parameters 
a = P.a;    % lattice constants 
s = P.s;  % spine width 
w = P.w;    % the beam width 
t = P.t;    % the rib width
th = P.th;  % thickness of the cavity 
d = P.d;    % side wall angle 

%% set the parameters 
model.param.set('a', a, 'lattice constant');
model.param.set('s', s, 'spine width');
model.param.set('w', w, 'beam width');
model.param.set('t', t, 'rib width');
model.param.set('th', th, 'thickness');
model.param.set('d', d, 'side wall angle');
model.param.set('airrad', '5*a');
model.param.set('kx', 'k*pi/a');
model.param.set('k', '1');
model.param.set('lbd0', '1.55[um]', 'wavelength band edge');

%% Create the components 
ucellcomp = model.modelNode.create('comp1');
ucellcomp.label('Unit cell FEM simulation');
ucellname = 'geom1';
ucellgeom = model.geom.create(ucellname,3);

%% create the unit cell
model.result.table.create('tbl1', 'Table');
model.result.table.create('tbl2', 'Table');
model.result.table.create('tbl3', 'Table');
model.result.table.create('tbl4', 'Table');

model.component('comp1').geom('geom1').geomRep('comsol');
model.component('comp1').geom('geom1').create('wp1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp1').set('quickz', '-th/2');
model.component('comp1').geom('geom1').feature('wp1').set('unite', true);
model.component('comp1').geom('geom1').feature('wp1').geom.create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').active(false);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').set('table', {'-a-a/2+t/2' 's/2';  ...
'-a/2-t/2' 's/2';  ...
'-a/2-t/2' 'w/2';  ...
'-a/2+t/2' 'w/2';  ...
'-a/2+t/2' 's/2';  ...
'a/2-t/2' 's/2';  ...
'a/2-t/2' 'w/2';  ...
'a/2+t/2' 'w/2';  ...
'a/2+t/2' 's/2';  ...
'a+a/2-t/2' 's/2';  ...
'a+a/2-t/2' '-s/2';  ...
'a/2+t/2' '-s/2';  ...
'a/2+t/2' '-w/2';  ...
'a/2-t/2' '-w/2';  ...
'a/2-t/2' '-s/2';  ...
'-a/2+t/2' '-s/2';  ...
'-a/2+t/2' '-w/2';  ...
'-a/2-t/2' '-w/2';  ...
'-a/2-t/2' '-s/2';  ...
'-a-a/2+t/2' '-s/2'});
model.component('comp1').geom('geom1').feature('wp1').geom.create('fil1', 'Fillet');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('fil1').active(false);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('fil1').set('radius', 'r1');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('fil1').selection('point').set('pol1(1)', [3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18]);
model.component('comp1').geom('geom1').feature('wp1').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('size', {'t' 'w'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('pos', {'-a/2' '0'});
model.component('comp1').geom('geom1').feature('wp1').geom.create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').set('size', {'t' 'w'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').set('pos', {'a/2' '0'});
model.component('comp1').geom('geom1').feature('wp1').geom.create('r3', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r3').set('size', {'a*2' 's'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r3').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r3').set('pos', [0 0]);
model.component('comp1').geom('geom1').create('ext1', 'Extrude');
model.component('comp1').geom('geom1').feature('ext1').set('extrudefrom', 'faces');
model.component('comp1').geom('geom1').feature('ext1').set('inputhandling', 'keep');
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', 'th', 0);
model.component('comp1').geom('geom1').feature('ext1').set('crossfaces', false);
model.component('comp1').geom('geom1').feature('ext1').selection('inputface').set('wp1.uni', [2 3 4 6 7 8]);
model.component('comp1').geom('geom1').create('wp3', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp3').set('quickplane', 'yz');
model.component('comp1').geom('geom1').feature('wp3').set('quickx', '-a');
model.component('comp1').geom('geom1').feature('wp3').set('unite', true);
model.component('comp1').geom('geom1').feature('wp3').geom.create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature('wp3').geom.feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature('wp3').geom.feature('pol1').set('table', {'w/2' '-th/2'; '(w-2*tan(d)*th)/2' 'th/2'; '-(w-2*tan(d)*th)/2' 'th/2'; '-w/2' '-th/2'});
model.component('comp1').geom('geom1').create('ext3', 'Extrude');
model.component('comp1').geom('geom1').feature('ext3').setIndex('distance', '2*a', 0);
model.component('comp1').geom('geom1').feature('ext3').selection('input').set({'wp3'});
model.component('comp1').geom('geom1').create('int1', 'Intersection');
model.component('comp1').geom('geom1').feature('int1').selection('input').set({'ext1' 'ext3'});
model.component('comp1').geom('geom1').create('wp4', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp4').set('quickplane', 'xz');
model.component('comp1').geom('geom1').feature('wp4').set('quicky', '-w/2');
model.component('comp1').geom('geom1').feature('wp4').set('unite', true);
model.component('comp1').geom('geom1').feature('wp4').geom.create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature('wp4').geom.feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature('wp4').geom.feature('pol1').set('table', {'-a/2-t/2' '-th/2'; '(-a/2-(1-2*tan(d)*th/t)*t/2)' 'th/2'; '-a/2+(1-2*tan(d)*th/t)*t/2' 'th/2'; '-a/2+t/2' '-th/2'});
model.component('comp1').geom('geom1').feature('wp4').geom.create('pol2', 'Polygon');
model.component('comp1').geom('geom1').feature('wp4').geom.feature('pol2').set('source', 'table');
model.component('comp1').geom('geom1').feature('wp4').geom.feature('pol2').set('table', {'a/2+t/2' '-th/2'; '-(-a/2-(1-2*tan(d)*th/t)*t/2)' 'th/2'; '-(-a/2+(1-2*tan(d)*th/t)*t/2)' 'th/2'; '-(-a/2+t/2)' '-th/2'});
model.component('comp1').geom('geom1').create('ext4', 'Extrude');
model.component('comp1').geom('geom1').feature('ext4').setIndex('distance', 'w', 0);
model.component('comp1').geom('geom1').feature('ext4').set('reverse', true);
model.component('comp1').geom('geom1').feature('ext4').selection('input').set({'wp4'});
model.component('comp1').geom('geom1').create('int2', 'Intersection');
model.component('comp1').geom('geom1').feature('int2').selection('input').set({'ext4' 'int1'});
model.component('comp1').geom('geom1').create('ext5', 'Extrude');
model.component('comp1').geom('geom1').feature('ext5').set('extrudefrom', 'faces');
model.component('comp1').geom('geom1').feature('ext5').setIndex('distance', 'th', 0);
model.component('comp1').geom('geom1').feature('ext5').selection('inputface').set('wp1.uni', [1 3 5 7 9]);
model.component('comp1').geom('geom1').create('wp2', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp2').set('quickplane', 'yz');
model.component('comp1').geom('geom1').feature('wp2').set('quickx', '-a');
model.component('comp1').geom('geom1').feature('wp2').set('unite', true);
model.component('comp1').geom('geom1').feature('wp2').geom.create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature('wp2').geom.feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature('wp2').geom.feature('pol1').set('table', {'-s/2' '-th/2'; 's/2' '-th/2'; '(s-2*tan(d)*th)/2' 'th/2'; '-(s-2*tan(d)*th)/2' 'th/2'});
model.component('comp1').geom('geom1').create('ext2', 'Extrude');
model.component('comp1').geom('geom1').feature('ext2').setIndex('distance', '2*a', 0);
model.component('comp1').geom('geom1').feature('ext2').selection('input').set({'wp2'});
model.component('comp1').geom('geom1').create('int3', 'Intersection');
model.component('comp1').geom('geom1').feature('int3').selection('input').set({'ext2' 'ext5'});
model.component('comp1').geom('geom1').create('wp5', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp5').set('quickplane', 'zy');
model.component('comp1').geom('geom1').feature('wp5').set('quickx', '-a/2');
model.component('comp1').geom('geom1').feature('wp5').set('unite', true);
model.component('comp1').geom('geom1').feature('wp5').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp5').geom.feature('r1').set('size', {'th' 'w'});
model.component('comp1').geom('geom1').feature('wp5').geom.feature('r1').set('pos', {'-th/2' '-w/2'});
model.component('comp1').geom('geom1').create('ext6', 'Extrude');
model.component('comp1').geom('geom1').feature('ext6').setIndex('distance', 'a/2', 0);
model.component('comp1').geom('geom1').feature('ext6').selection('input').set({'wp5'});
model.component('comp1').geom('geom1').create('wp6', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp6').set('quickplane', 'zy');
model.component('comp1').geom('geom1').feature('wp6').set('quickx', 'a/2');
model.component('comp1').geom('geom1').feature('wp6').set('unite', true);
model.component('comp1').geom('geom1').feature('wp6').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp6').geom.feature('r1').set('size', {'th' 'w'});
model.component('comp1').geom('geom1').feature('wp6').geom.feature('r1').set('pos', {'-th/2' '-w/2'});
model.component('comp1').geom('geom1').create('ext7', 'Extrude');
model.component('comp1').geom('geom1').feature('ext7').setIndex('distance', 'a/2', 0);
model.component('comp1').geom('geom1').feature('ext7').set('reverse', true);
model.component('comp1').geom('geom1').feature('ext7').selection('input').set({'wp6'});
model.component('comp1').geom('geom1').create('co1', 'Compose');
model.component('comp1').geom('geom1').feature('co1').set('formula', 'int2+int3-ext6-ext7');
model.component('comp1').geom('geom1').create('wp7', 'WorkPlane');
if ~P.run_optical
    model.component('comp1').geom('geom1').feature('wp7').active(false);
end
model.component('comp1').geom('geom1').feature('wp7').set('quickplane', 'xz');
model.component('comp1').geom('geom1').feature('wp7').set('unite', true);
model.component('comp1').geom('geom1').feature('wp7').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp7').geom.feature('r1').set('size', {'a' 'th'});
model.component('comp1').geom('geom1').feature('wp7').geom.feature('r1').set('pos', {'-a/2' '-th/2'});
model.component('comp1').geom('geom1').create('ext8', 'Extrude');
if ~P.run_optical
    model.component('comp1').geom('geom1').feature('ext8').active(false);
end
model.component('comp1').geom('geom1').feature('ext8').setIndex('distance', 'w/2', 0);
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
model.component('comp1').geom('geom1').feature('wp8').geom.feature('r1').set('size', {'a' '5*a'});
model.component('comp1').geom('geom1').feature('wp8').geom.feature('r1').set('pos', {'-a/2' '0'});
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