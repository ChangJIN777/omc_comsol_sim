function [model,P] = buildBoomerangStrip_3D(model,P)
%
% buildBoomerangUnitCell.m
%
% Model exported on Jul 21 2024, 16:57 by COMSOL 6.1.0.357.

%% read the input parameters 
% the parameters for the boomerang unit cell
a = P.a;        % lattice constant 
w = P.w;        % the width of the boomerang unit cell 
r = P.r;        % the length of the boomerang unit cell 
th = P.th;        % the thickness of the cavity 
% the parameters for the center unit cell 
wo = P.wo;
wi = P.wi;
ho = P.ho;
hi = P.hi; 
% P.d = sqrt(3)*a/4;
b_wvg = 0;
% fillet radius
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

ucellplane = ucellWP.geom.feature.create('r1', 'Rectangle');
ucellplane.label('Base plane');
ucellplane.set('pos', [-a/2 0]);
ucellplane.set('base','corner');
ucellplane.set('size',[a sqrt(3)*(2+1/2)*a+b_wvg]);
hole_label_list = {};

% hole1_labelList = add_BoomerangHole_center(ucellWP,P,1);
% hole_label_list = [hole_label_list,hole1_labelList];
% hole1_center = [0,P.d];

% adding center unit cell v2 
hole1_center = [0,P.d];
holeCenter_lableList = add_BoomerangHole_center_v2(ucellWP,hole1_center,P,1);
hole_label_list = [hole_label_list,holeCenter_lableList];

hole2_center = [a/2 sqrt(3)*a+b_wvg];
hole2_labelList = add_BoomerangHole(ucellWP,hole2_center,P,2);
hole_label_list = [hole_label_list,hole2_labelList];

hole3_center = [-a/2 sqrt(3)*a+b_wvg];
hole3_labelList = add_BoomerangHole(ucellWP,hole3_center,P,3);
hole_label_list = [hole_label_list,hole3_labelList];

hole4_center = [0 (3/2)*sqrt(3)*a+b_wvg];
hole4_labelList = add_BoomerangHole(ucellWP,hole4_center,P,4);
hole_label_list = [hole_label_list,hole4_labelList];

hole5_center = [-a/2 2*sqrt(3)*a+b_wvg];
hole5_labelList = add_BoomerangHole(ucellWP,hole5_center,P,5);
hole_label_list = [hole_label_list,hole5_labelList];

hole6_center = [a/2 2*sqrt(3)*a+b_wvg];
hole6_labelList = add_BoomerangHole(ucellWP,hole6_center,P,6);
hole_label_list = [hole_label_list,hole6_labelList];

hole7_center = [0 sqrt(3)*a*(5/2)+b_wvg];
hole7_labelList = add_BoomerangHole(ucellWP,hole7_center,P,7);
hole_label_list = [hole_label_list,hole7_labelList];

hole8_center = [a/2 sqrt(3)*a*3+b_wvg];
hole8_labelList = add_BoomerangHole(ucellWP,hole8_center,P,8);
hole_label_list = [hole_label_list,hole8_labelList];

hole9_center = [-a/2 sqrt(3)*a*3+b_wvg];
hole9_labelList = add_BoomerangHole(ucellWP,hole9_center,P,9);
hole_label_list = [hole_label_list,hole9_labelList];

hole10_center = [0 sqrt(3)*a*(7/2)+b_wvg];
hole10_labelList = add_BoomerangHole(ucellWP,hole10_center,P,10);
hole_label_list = [hole_label_list,hole10_labelList];

hole11_center = [-a/2 sqrt(3)*a*4+b_wvg];
hole11_labelList = add_BoomerangHole(ucellWP,hole11_center,P,11);
hole_label_list = [hole_label_list,hole11_labelList];

hole12_center = [a/2 sqrt(3)*a*4+b_wvg];
hole12_labelList = add_BoomerangHole(ucellWP,hole12_center,P,12);
hole_label_list = [hole_label_list,hole12_labelList];

hole13_center = [0 sqrt(3)*a*(9/2)+b_wvg];
hole13_labelList = add_BoomerangHole(ucellWP,hole13_center,P,13);
hole_label_list = [hole_label_list,hole13_labelList];

compose1 = ucellWP.geom.feature.create('co1', 'Compose');
compose_string = 'r1';
for i=1:length(hole_label_list)
    compose_string = strcat(compose_string,'-');
    compose_string = strcat(compose_string,hole_label_list(i));
end
    
compose1.set('formula', compose_string);

selection_width = 10e-9;
% add fillets 
addFillet(P,ucellWP.geom,1,hole1_center,selection_width*2)
addFillet(P,ucellWP.geom,2,hole2_center,selection_width)
addFillet(P,ucellWP.geom,3,hole3_center,selection_width)
addFillet(P,ucellWP.geom,4,hole4_center,selection_width)
addFillet(P,ucellWP.geom,5,hole5_center,selection_width)
addFillet(P,ucellWP.geom,6,hole6_center,selection_width)
addFillet(P,ucellWP.geom,7,hole7_center,selection_width)
addFillet(P,ucellWP.geom,8,hole8_center,selection_width)
addFillet(P,ucellWP.geom,9,hole9_center,selection_width)
addFillet(P,ucellWP.geom,10,hole10_center,selection_width)
addFillet(P,ucellWP.geom,11,hole11_center,selection_width)
addFillet(P,ucellWP.geom,12,hole12_center,selection_width)
addFillet(P,ucellWP.geom,13,hole13_center,selection_width)

model.component('comp1').geom('geom1').run('wp1');

model.component('comp1').geom('geom1').feature.create('ext1', 'Extrude');
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', th, 0);
model.component('comp1').geom('geom1').run('ext1');
% extrude.setIndex('distance', th, 0);
% extrude.selection('input').set({'wp1'});

% add the air disk 
if P.add_airDisk
    airDiskH = P.airDiskH; % the height the of air disk
    airPlane = ucellgeom.feature.create('wpair', 'WorkPlane');
    airPlane.set('quickplane', 'zx');
    airPlane.set('unite', true);
    diskBase = airPlane.geom.feature.create('r1', 'Rectangle');
    diskBase.set('pos', [0 -a/2]);
    diskBase.set('size', [airDiskH a]);
    revolve = ucellgeom.feature.create('rev1', 'Revolve');
    revolve.set('angle2', -180);
    revolve.selection('input').set({'wpair'});
end

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
    ZsymSel.set('ymin', -delta).set('ymax', (sqrt(3)*(4+1/2)*a+b_wvg)+delta);
    ZsymSel.set('zmin', -delta).set('zmax', delta);
    ZsymSel.set('entitydim', 2).set('condition', 'allvertices'); % only want beam surface, exclude air holes
    ucellgeom.runCurrent;
    inds = model.selection([ucellname,'_ZsymSel']).inputEntities();
    P.bndSel.Zsym = inds';

end 

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

function hole_rec3_label_list = add_BoomerangHole(ucellWP,hole_center,P,holeIndex)
    hole_rec1_label = sprintf('h%dr1',holeIndex);
    hole_rec2_label = sprintf('h%dr2',holeIndex);
    hole_rec3_label = sprintf('h%dr3',holeIndex);
    hole_rec3_label_list = {hole_rec1_label,hole_rec2_label,hole_rec3_label};
    a = P.a;
    r = P.r;
    w = P.w;
    hole_rec_1 = ucellWP.geom.feature.create(hole_rec1_label, 'Rectangle');
    hole_rec_1.set('pos', [hole_center(1) r/2+hole_center(2)]);
    hole_rec_1.set('base', 'center');
    hole_rec_1.set('size', [w r]);
    hole_rec_2 = ucellWP.geom.feature.create(hole_rec2_label, 'Rectangle');
    hole_rec_2.set('pos', [-sqrt(3)*r/4+hole_center(1) -r/4+hole_center(2)]);
    hole_rec_2.set('base', 'center');
    hole_rec_2.set('rot', 120);
    hole_rec_2.set('size', [w r]);
    hole_rec_3 = ucellWP.geom.feature.create(hole_rec3_label, 'Rectangle');
    hole_rec_3.set('pos', [sqrt(3)*r/4+hole_center(1) -r/4+hole_center(2)]);
    hole_rec_3.set('base', 'center');
    hole_rec_3.set('rot', 240);
    hole_rec_3.set('size', [w r]);
end

function hole_rec3_label_list = add_BoomerangHole_center(ucellWP,P,holeIndex)
    hole_rec1_label = sprintf('h%dr1',holeIndex);
    hole_rec2_label = sprintf('h%dr2',holeIndex);
    hole_rec3_label = sprintf('h%dr3',holeIndex);
    hole_rec3_label_list = {hole_rec1_label,hole_rec2_label,hole_rec3_label};
    wo = P.wo;
    wi = P.wi;
    ho = P.ho;
    hi = P.hi;
    d = P.d;
    hole_rec_1 = ucellWP.geom.feature.create(hole_rec1_label, 'Rectangle');
    hole_rec_1.set('pos', [-wo/2 hi+d]);
    hole_rec_1.set('base', 'corner');
    hole_rec_1.set('size', [wo ho-hi]);
    hole_rec_2 = ucellWP.geom.feature.create(hole_rec2_label, 'Rectangle');
    hole_rec_2.set('pos', [-wo/2 d]);
    hole_rec_2.set('base', 'corner');
    hole_rec_2.set('size', [(wo-wi)/2 ho]);
    hole_rec_3 = ucellWP.geom.feature.create(hole_rec3_label, 'Rectangle');
    hole_rec_3.set('pos', [wi/2 d]);
    hole_rec_3.set('base', 'corner');
    hole_rec_3.set('size', [(wo-wi)/2 ho]);
end

function hole_rec3_label_list = add_BoomerangHole_center_v2(ucellWP,hole_center,P,holeIndex)
    hole_rec1_label = sprintf('h%dr1',holeIndex);
    hole_rec2_label = sprintf('h%dr2',holeIndex);
    hole_rec3_label = sprintf('h%dr3',holeIndex);
    hole_rec4_label = sprintf('h%dr4',holeIndex);
    hole_rec5_label = sprintf('h%dr5',holeIndex);
    hole_rec3_label_list = {hole_rec1_label,hole_rec2_label,hole_rec3_label,hole_rec4_label,hole_rec5_label};
    a = P.a;
    r = P.r;
    w = P.w;
    h = P.h;
    d1 = P.d1;
    hole_rec_1 = ucellWP.geom.feature.create(hole_rec1_label, 'Rectangle');
    hole_rec_1.set('pos', [hole_center(1) r/2+hole_center(2)]);
    hole_rec_1.set('base', 'center');
    hole_rec_1.set('size', [w r]);
    hole_rec_2 = ucellWP.geom.feature.create(hole_rec2_label, 'Rectangle');
    hole_rec_2.set('pos', [-sqrt(3)*r/4+hole_center(1) -r/4+hole_center(2)]);
    hole_rec_2.set('base', 'center');
    hole_rec_2.set('rot', 120);
    hole_rec_2.set('size', [w r]);
    hole_rec_3 = ucellWP.geom.feature.create(hole_rec3_label, 'Rectangle');
    hole_rec_3.set('pos', [sqrt(3)*r/4+hole_center(1) -r/4+hole_center(2)]);
    hole_rec_3.set('base', 'center');
    hole_rec_3.set('rot', 240);
    hole_rec_3.set('size', [w r]);
    hole_rec_4 = ucellWP.geom.feature.create(hole_rec4_label, 'Rectangle');
    hole_rec_4.set('pos', [sqrt(3)*r/2+hole_center(1)+w/4 -r/2+hole_center(2)+sqrt(3)*w/4]);
    hole_rec_4.set('base', 'corner');
    hole_rec_4.set('rot', 180);
    hole_rec_4.set('size', [d1 h]);
    hole_rec_5 = ucellWP.geom.feature.create(hole_rec5_label, 'Rectangle');
    hole_rec_5.set('pos', [-sqrt(3)*r/2+hole_center(1)-w/4 -r/2+hole_center(2)+sqrt(3)*w/4-h]);
    hole_rec_5.set('base', 'corner');
    hole_rec_5.set('size', [d1 h]);
    
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
    disksel2.set('r', r+selection_width);
    disksel2.set('rin', r-selection_width);
    disksel2.set('condition', 'allvertices');
    fil2 = ucellgeom.feature.create(fillet2_label, 'Fillet');
    fil2.set('radius', r2);
    fil2.selection('point').named(disksel2_label);
end

function addFillet_center(P,ucellgeom,hole_indx,hole_pos,selection_width)
    wi = P.wi;
    wo = P.wo;
    hi = P.hi;
    ho = P.ho;
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
    disksel1.set('r', sqrt(hi^2+(wi/2)^2)+selection_width);
    disksel1.set('rin', wi/2-selection_width);
    disksel1.set('condition', 'allvertices');
    fil1 = ucellgeom.feature.create(fillet1_label, 'Fillet');
    fil1.set('radius', r1);
    fil1.selection('point').named(disksel1_label);
    disksel2 = ucellgeom.feature.create(disksel2_label, 'DiskSelection');
    disksel2.set('entitydim', 0);
    disksel2.set('posy', hole_pos(2));
    disksel2.set('posx', hole_pos(1));
    disksel2.set('r', sqrt(ho^2+(wo/2)^2)+selection_width);
    disksel2.set('rin', wo/2-selection_width);
    disksel2.set('condition', 'allvertices');
    fil2 = ucellgeom.feature.create(fillet2_label, 'Fillet');
    fil2.set('radius', r2);
    fil2.selection('point').named(disksel2_label);
end