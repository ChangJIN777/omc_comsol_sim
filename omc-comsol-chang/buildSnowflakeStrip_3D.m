function [model,P] = buildSnowflakeStrip_3D(model,P)
%
% buildBoomerangUnitCell.m
%
% Model exported on Jul 21 2024, 16:57 by COMSOL 6.1.0.357.

%% read the input parameters 
% the parameters for the boomerang unit cell
a = P.a;        % lattice constant 
w = P.w;        % the width of the boomerang unit cell 
r = P.r;        % the length of the boomerang unit cell 
b = P.b;
th = P.th;        % the thickness of the cavity 
% the parameters for the center unit cell 
wo = P.wo;
wi = P.wi;
ho = P.ho;
hi = P.hi; 
d = P.d;
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

%% create unit cell with snowflake geometry
ucellWP = ucellgeom.feature.create('wp1', 'WorkPlane');
ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
% ucellWP.set('unite', true);

ucellplane = ucellWP.geom.feature.create('pol1', 'Polygon');
ucellplane.label('Base plane');
ucellplane.set('source', 'table');
ucellplane.set('table', [-a/2 0; a/2 0; a/2 b+sqrt(3)*a+d/2; -a/2 b+sqrt(3)*a+d/2; -a/2 0]);
hole_label_list = {};

% add the holes for the strip
% adding center unit cell v2 
hole1_center = [-a/2 b+d/2];
holeCenter_lableList = add_SnowflakeHole(ucellWP,hole1_center,P,1);
hole_label_list = [hole_label_list,holeCenter_lableList];

hole2_center = [a/2 b+d/2];
hole2_labelList = add_SnowflakeHole(ucellWP,hole2_center,P,2);
hole_label_list = [hole_label_list,hole2_labelList];

hole3_center = [0 b+sqrt(3)*a/2+d/2];
hole3_labelList = add_SnowflakeHole(ucellWP,hole3_center,P,3);
hole_label_list = [hole_label_list,hole3_labelList];

hole4_center = [-a/2 b+sqrt(3)*a+d/2];
hole4_labelList = add_SnowflakeHole(ucellWP,hole4_center,P,4);
hole_label_list = [hole_label_list,hole4_labelList];

hole5_center = [a/2 b+sqrt(3)*a+d/2];
hole5_labelList = add_SnowflakeHole(ucellWP,hole5_center,P,5);
hole_label_list = [hole_label_list,hole5_labelList];

% hole6_center = [0 b+sqrt(3)*a*(3/2)+d/2];
% hole6_labelList = add_SnowflakeHole(ucellWP,hole6_center,P,6);
% hole_label_list = [hole_label_list,hole6_labelList];
% 
% hole7_center = [a/2 b+sqrt(3)*a*2+d/2];
% hole7_labelList = add_SnowflakeHole(ucellWP,hole7_center,P,7);
% hole_label_list = [hole_label_list,hole7_labelList];
% 
% hole8_center = [-a/2 b+sqrt(3)*a*2+d/2];
% hole8_labelList = add_SnowflakeHole(ucellWP,hole8_center,P,8);
% hole_label_list = [hole_label_list,hole8_labelList];

% hole9_center = [a/2-(a/2+w/2) b+a*sqrt(3)/2+sqrt(3)*a/2];
% hole9_labelList = add_SnowflakeHole(ucellWP,hole9_center,P,9);
% hole_label_list = [hole_label_list,hole9_labelList];
% 
% hole10_center = [a/2-a-(a/2+w/2) b+a*sqrt(3)/2+sqrt(3)*a/2];
% hole10_labelList = add_SnowflakeHole(ucellWP,hole10_center,P,10);
% hole_label_list = [hole_label_list,hole10_labelList];




%% create unit cell with boomerang geometry in the lower cavity
% ucellWP_lower = ucellgeom.create('wp_lower', 'WorkPlane');
% ucellWP_lower.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
% ucellplane_lower = ucellWP_lower.geom.feature.create('r_base', 'Rectangle');
% ucellplane_lower.set('pos', [-a/2 0]);
% ucellplane_lower.set('size', [a b]);

rec1 = ucellWP.geom.feature.create('r1_base', 'Rectangle');
rec1.set('pos', [0 d/2+(ho-hi)/2+hi]);
rec1.set('base', 'center');
rec1.set('size', [wo ho-hi]);
hole_label_list = [hole_label_list,'r1_base'];

rec2 = ucellWP.geom.feature.create('r2_base', 'Rectangle');
rec2.set('pos', [-wo/2+(wo-wi)/4 d/2+ho/2]);
rec2.set('base', 'center');
rec2.set('size', [(wo-wi)/2 ho]);
hole_label_list = [hole_label_list,'r2_base'];

rec3 = ucellWP.geom.feature.create('r3_base', 'Rectangle');
rec3.set('pos', [wo/2-(wo-wi)/4 d/2+ho/2]);
rec3.set('base', 'center');
rec3.set('size', [(wo-wi)/2 ho]);
hole_label_list = [hole_label_list,'r3_base'];

%% compose the 2D geometry
compose1 = ucellWP.geom.feature.create('co1', 'Compose');
compose_string = 'pol1';
for i=1:length(hole_label_list)
    compose_string = strcat(compose_string,'-');
    compose_string = strcat(compose_string,hole_label_list(i));
end
    
compose1.set('formula', compose_string);

% compose1 = ucellWP.geom.feature.create('co1', 'Compose');
% compose1.set('formula', 'r_base-r1-r2-r3');

%% add fillets 
selection_width = w/2;
% add fillets 
addFillet(P,ucellWP.geom,1,hole1_center,selection_width*2);
addFillet(P,ucellWP.geom,2,hole2_center,selection_width);
addFillet(P,ucellWP.geom,3,hole3_center,selection_width);
addFillet(P,ucellWP.geom,4,hole4_center,selection_width);
addFillet(P,ucellWP.geom,5,hole5_center,selection_width);
% addFillet(P,ucellWP.geom,6,hole6_center,selection_width);
% addFillet(P,ucellWP.geom,7,hole7_center,selection_width);
% addFillet(P,ucellWP.geom,8,hole8_center,selection_width);
% addFillet(P,ucellWP.geom,9,hole9_center,selection_width);
% addFillet(P,ucellWP.geom,10,hole10_center,selection_width);

% fillet1 = ucellWP.geom.feature.create('fil1', 'Fillet');
% fillet1.set('radius', r1);
% fillet1.selection('point').set('co1(1)', [40 42 61 62 128 129 149 151]);
% fillet2 = ucellWP.geom.feature.create('fil2', 'Fillet');
% fillet2.set('radius', r2);
% fillet2.selection('point').set('fil1(1)', [10 13]);
% ucellgeom.runAll;

%% extrude 
model.component('comp1').geom('geom1').run('wp1');

model.component('comp1').geom('geom1').feature.create('ext1', 'Extrude');
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', th, 0);
model.component('comp1').geom('geom1').run('ext1');

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
    symZPlane.set('table', [-a/2 0; a/2 0; a/2 b+sqrt(3)*a+d/2; -a/2 b+sqrt(3)*a+d/2; -a/2 0]);

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
    ZsymSel.set('ymin', -delta).set('ymax', b+sqrt(3)*a+d/2+delta);
    ZsymSel.set('zmin', -delta).set('zmax', delta);
    ZsymSel.set('entitydim', 2).set('condition', 'allvertices'); % only want beam surface, exclude air holes
    ZsymSel.set('condition', 'inside');
    ucellgeom.runCurrent;
    inds = model.selection([ucellname,'_ZsymSel']).inputEntities();
    P.bndSel.Zsym = inds';

end

%% Making selections (with box select)
mphgeom(model);
selection_width = 1e-9; 

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

function hole_rec_label_list = add_SnowflakeHole(ucellWP,hole_center,P,holeIndex)
    %% add a snowflake hole with specified location and labels 
    hole_rec1_label = sprintf('h%dr1',holeIndex);
    hole_rec2_label = sprintf('h%dr2',holeIndex);
    hole_rec3_label = sprintf('h%dr3',holeIndex);
    
    hole_rec_label_list = {hole_rec1_label,hole_rec2_label,hole_rec3_label};

    a = P.a;
    r = P.r;
    w = P.w;

    c1_rec_1 = ucellWP.geom.feature.create(hole_rec1_label, 'Rectangle');
    c1_rec_1.set('pos', [hole_center(1) hole_center(2)]);
    c1_rec_1.set('base', 'center');
    c1_rec_1.set('size', [r*2 w]);
    
    c1_rec_2 = ucellWP.geom.feature.create(hole_rec2_label, 'Rectangle');
    c1_rec_2.set('pos', [hole_center(1) hole_center(2)]);
    c1_rec_2.set('rot', 60);
    c1_rec_2.set('base', 'center');
    c1_rec_2.set('size', [r*2 w]);
    
    c1_rec_3 = ucellWP.geom.feature.create(hole_rec3_label, 'Rectangle');
    c1_rec_3.set('pos', [hole_center(1) hole_center(2)]);
    c1_rec_3.set('rot', 120);
    c1_rec_3.set('base', 'center');
    c1_rec_3.set('size', [r*2 w]);

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
