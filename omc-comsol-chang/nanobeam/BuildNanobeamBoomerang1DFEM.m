% Function to construct geometry (nanobeam, air/mech PML etc.) for FEM
% simulations
% Cleaven Chia, 11/21/18

function [model,P] = BuildNanobeamBoomerang1DFEM(model,P)

%% extract geometry parameters from P
a = P.a;        % lattice constant 
w = P.w;        % the width of the hole 
th = P.th;        % the height of the hole
r = P.r;        % the height of the hole 
wo = P.wo;        % the height of the hole in the lower portion
wi = P.wi;      
ho = P.ho;      
hi = P.hi;
d = P.d;        % the width of the hole in the lower portion
b = P.b;        % 
r1 = P.r1;      % the fillet radius of the edges of the hole 
r2 = P.r2;      % the fillet radius of the center of the hole 
% abssym = abs(P.mbevenz);    % symmetry in the z direction 

% asym = P.asym;

% if isfield(P,'asymCav') && P.asymCav
%     geom = P.geom;      % the geometric parameters associated with the cavity
%     beamLen = P.beamLen;
%     P.xc = P.beamLenHalfL;
% else
%     geom = P.geomHalf;
%     beamLen = P.beamLenHalf;
%     P.xc = 0;
% end
geom = P.geom;

% cavity geometries 
MN_left = P.MN_left;
MN_right = P.MN_right;
TN = P.TN;
% the parameters in the defect region 
a_ctr = P.a_ctr;                 % for taperTo = 'custom': lattice constant of center hole
ho_ctr = P.ho_ctr;                    % for taperTo = 'custom': hole height of center hole
hi_ctr = P.hi_ctr;                    % for taperTo = 'custom': hole height of center hole
wo_ctr = P.wo_ctr;                    % for taperTo = 'custom': hole height of center hole
wi_ctr = P.wi_ctr;                    % for taperTo = 'custom': hole height of center hole
% get the parameter lists of the cavity 
nholes = size(geom,1); % the number of unit cells 
wo_list = geom(:,1)';
wi_list = geom(:,2)';
ho_list = geom(:,3)';
hi_list = geom(:,4)';
a_list = geom(:,5)';
% calculate the location of the center of each super unit cell
xpos = zeros(1,nholes);
cci = MN_left + TN;
for index = cci+1:nholes
    xpos(index) = xpos(index-1)+a_list(index-1)/2+a_list(index)/2;
end
for i = 2:cci
    index = cci+1-i;
    xpos(index) = xpos(index+1)-a_list(index+1)/2-a_list(index)/2;
end
beamLen = abs(xpos(nholes)-xpos(1));
P.beamLen = beamLen;
% if using rectangular cross section, use half of height due to symmetry
if strcmp(P.xsect,'rect')
    thi = P.th/2;
end

%% Create component
comp = model.modelNode.create('comp');
comp.label('Nanobeam FEM simulation');
% if P.solveOpt
%     beamlabel = 'Nanobeam with air cylinder';
%     beamname = 'beamCyl';
% else
%     beamlabel = 'Nanobeam';
%     beamname = 'beam';
% end

% do not worry about the optical mode for now 
P.geomname = 'slab';
slabgeom = model.geom.create(P.geomname, 3);
slabgeom.label('slab');

%% Create beam with rectangular cross-section
slabWP = slabgeom.feature.create('slabWP', 'WorkPlane');
slabWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', 0);

slabplane = slabWP.geom.feature.create('slabplane', 'Rectangle');
slabplane.set('type', 'solid').set('base', 'corner');
slabplane.set('pos', [-beamLen/2 0]).set('size', [beamLen sqrt(3)*2*a]);
slabgeom.runCurrent;
    
holeList = {};
holeFormula = [];
extrude_labels = {};
% add the supercells of the cavities one by one
for k = 1:nholes
    % add the boomerang unit cells 
    xloc = xpos(k);
    if P.addshield
        label_list = buildBoomerangCells(P,slabgeom,xloc,k);
    else 
        label_list = {};
    end
    holeList = [holeList,label_list];
    % add the lower cavity region 
    P.a = a_list(k);
    P.ho = ho_list(k);
    P.hi = hi_list(k);
    P.wo = wo_list(k);
    P.wi = wi_list(k);
    lowerCavity = buildLowerCell(slabgeom,P,k,xloc);
    holeList = [holeList,{lowerCavity}];
end
    
% % compose workplane
% beamComp = beamWP.geom.feature.create('beamComp', 'Compose');
% beamComp.selection('input').set(['beamplane',holeList]);
% beamComp.set('formula', ['beamplane',holeFormula]).set('intbnd',false);
% beamgeom.runCurrent;
% 
% 
% % extrude
% finBeamTag = 'beamHoles';
% beamHoles = beamgeom.feature.create(finBeamTag, 'Extrude');
% beamHoles.set('distance', thi);
% beamgeom.runCurrent;  
displayBeamStr = 'Nanobeam, rectangular cross-section';

% % track max dimensions for selections
% totLen = beamLen;
% maxWid = wid;
% maxThi = thi;

% track max dimensions for selections
totLen = beamLen;
maxWid = a;
maxThi = thi;



%% create the intersect 
compose = slabgeom.feature.create('col','Compose');
% compose the formula for making the structure 
compose_formula = ['('];
for i=1:length(holeList)
    if i == length(holeList)
        compose_formula = [compose_formula,holeList{i},')'];
    else
        compose_formula = [compose_formula,holeList{i},'+'];
    end
end
holeList = [{'slabWP'}, holeList];
compose_formula = ['slabWP*',compose_formula];
compose.selection('input').set(holeList);
compose.set('formula', compose_formula);

% use box selection for extruding
% box selection for the boundary condition 
base_Selection = slabgeom.feature.create('2d_structures', 'BoxSelection');
base_Selection.set('entitydim', -1);
% x_boundary_disksel_r.set('xmin', a/2-selection_width/2);
% x_boundary_disksel_r.set('xmax', a/2+selection_width/2);
base_Selection.set('inputent', 'all');
base_Selection.set('condition', 'inside');

% extrude 
extrude = slabgeom.feature.create('ext1', 'Extrude');
extrude.set('extrudefrom', 'faces');
extrude.selection('inputface').named('2d_structures');
extrude.setIndex('distance', th/2, 0);

%% Run geometry
slabgeom.run;

%% to implement: adding phononic mirrors
% use createNanobeamGeom
% then update total length

%% Making selections (with box select) for boundary conditions 
% mphgeom(model);
selection_width = 2.5e-9;

% box selection for the boundary condition 
x_boundary_disksel_r = slabgeom.feature.create('x_boundary_boxsel_r', 'BoxSelection');
x_boundary_disksel_r.set('entitydim', 2);
x_boundary_disksel_r.set('xmin', beamLen/2-selection_width/2);
x_boundary_disksel_r.set('xmax', beamLen/2+selection_width/2);
x_boundary_disksel_r.set('inputent', 'all');
x_boundary_disksel_r.set('condition', 'inside');

x_boundary_disksel_l = slabgeom.feature.create('x_boundary_boxsel_l', 'BoxSelection');
x_boundary_disksel_l.set('entitydim', 2);
x_boundary_disksel_l.set('xmin', -beamLen/2-selection_width/2);
x_boundary_disksel_l.set('xmax', -beamLen/2+selection_width/2);
x_boundary_disksel_l.set('inputent', 'all');
x_boundary_disksel_l.set('condition', 'inside');

slabgeom.selection.create('xboundaries','CumulativeSelection');
slabgeom.selection('xboundaries').label('Cumulative Selection x boundaries');
x_boundary_disksel_r.set('contributeto','xboundaries');
x_boundary_disksel_l.set('contributeto','xboundaries');

y_boundary_symmetry = slabgeom.feature.create('y_boundary_symmetry', 'BoxSelection');
y_boundary_symmetry.set('entitydim', 2);
y_boundary_symmetry.set('ymin', -selection_width/2);
y_boundary_symmetry.set('ymax', selection_width/2);
y_boundary_symmetry.set('inputent', 'all');
y_boundary_symmetry.set('condition', 'inside');

slabgeom.selection.create('yboundaries','CumulativeSelection');
slabgeom.selection('yboundaries').label('Cumulative Selection y boundaries');
y_boundary_symmetry.set('contributeto','yboundaries');

z_boundary_symmetry = slabgeom.feature.create('z_boundary_symmetry', 'BoxSelection');
z_boundary_symmetry.set('entitydim', 2);
z_boundary_symmetry.set('zmin', -selection_width/2);
z_boundary_symmetry.set('zmax', selection_width/2);
z_boundary_symmetry.set('inputent', 'all');
z_boundary_symmetry.set('condition', 'inside');

slabgeom.selection.create('zboundaries','CumulativeSelection');
slabgeom.selection('zboundaries').label('Cumulative Selection z boundaries');
z_boundary_symmetry.set('contributeto','zboundaries');

% if we dont add 2D shield 
if P.addshield~=1
    y_boundary_fixed = slabgeom.feature.create('y_boundary_fixed', 'BoxSelection');
    y_boundary_fixed.set('entitydim', 2);
    y_boundary_fixed.set('ymin', b-selection_width/2);
    y_boundary_fixed.set('ymax', b+selection_width/2);
    y_boundary_fixed.set('inputent', 'all');
    y_boundary_fixed.set('condition', 'inside');
    
    slabgeom.selection.create('yboundaries_fixed','CumulativeSelection');
    slabgeom.selection('yboundaries_fixed').label('Cumulative Selection y fixed boundaries');
    y_boundary_fixed.set('contributeto','yboundaries_fixed');
end

%% Create air cylinder around beam
displayCylStr = '';
if isfield(P,'airrad') && P.solveOpt
    cyl_cut_wp = slabgeom.feature.create('cyl_cut_wp', 'WorkPlane');
    cyl_cut_wp.set('planetype', 'quick').set('quickplane', 'xz');
    cyl_cut = cyl_cut_wp.geom.feature.create('cyl_cut','Rectangle');
    cyl_cut.set('type','solid').set('pos',[0,0]).set('size',[totLen,P.airrad]);

    air_cyl = slabgeom.feature.create('air_cyl', 'Revolve');
    air_cyl.selection('input').set('cyl_cut_wp');

    air_cyl.set('angle1','-180');
    air_cyl.set('axis',[1,0]).set('pos',[0,0]).set('angle2','0');

    % Compose final geometry
    beamHolesAir = slabgeom.feature.create('beamHolesAir', 'Compose');
    beamHolesAir.selection('input').set('symZComp');
    beamHolesAir.selection('input').set('air_cyl');
    beamHolesAir.set('formula', ['air_cyl-','symZComp']);
    slabgeom.run;
    displayCylStr = ', air cylinder';
    finBeamTag = 'beamHolesAir';

    % track max dimensions for selections
    totLen = beamLen;
    maxWid = max(maxWid,2*P.airrad);
    maxThi = max(maxThi,2*P.airrad);
end

% display(['Geometry created - ',displayBeamStr,displayCylStr,displayPMLStr]);

%% Create domain selections after full geometry is constructed

% Create selection with beam only
% define box slightly larger and fully containing full beam volume
delta = 10e-9; 
beamSel = slabgeom.create('beamSel', 'BoxSelection');
beamSel.set('xmin', -beamLen/2-delta).set('xmax', beamLen/2 + delta);
beamSel.set('ymin', -delta).set('ymax', sqrt(3)*P.a*4 + delta);
beamSel.set('zmin', -delta).set('zmax', thi + delta);
beamSel.set('entitydim', 3).set('condition', 'inside');
slabgeom.runCurrent;
inds = model.selection([P.geomname,'_beamSel']).inputEntities();
P.domSel.beam = inds';

% Create selection with cylinder only
% - define box slightly larger and fully containing full cylinder volume
% then take difference selection with beam
if isfield(P,'airrad') && P.solveOpt
    delta = 10e-9; 
    beamCylSel = slabgeom.create('beamCylSel', 'BoxSelection');
    beamCylSel.set('xmin', -delta).set('xmax', beamLen + delta);
    beamCylSel.set('ymin', -delta).set('ymax', P.airrad + delta);
    beamCylSel.set('zmin', -P.airrad-delta).set('zmax', P.airrad + delta);
    beamCylSel.set('condition', 'inside');
    slabgeom.runCurrent;
    
    cylSel = slabgeom.create('cylSel', 'DifferenceSelection');
    cylSel.set('entitydim',3).set('add','beamCylSel').set('subtract','beamSel');
    slabgeom.runCurrent;
    inds = model.selection([P.geomname,'_cylSel']).inputEntities();
    P.domSel.cyl = inds';
end

% return the model 
out = model;
end 

%% define the subfunction for make variout unit cell geometries 
function ucellWP = buildBoomerangCell(P,ucellgeom,workPlaneName, loc)
    %% create unit cell with boomerang geometry
    ucellWP = ucellgeom.feature.create(workPlaneName, 'WorkPlane');
    ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', 0);
    a = P.a;
    w = P.w;
    r = P.r;
    ucellplane = ucellWP.geom.feature.create('pol1', 'Polygon');
    ucellplane.label('Base plane');
    ucellplane.set('source', 'table');
    ucellplane.set('table', [0 0; a/2 (a/2)*sqrt(3); a*(3/2) (a/2)*sqrt(3); a 0; 0 0]);
    rec_1 = ucellWP.geom.feature.create('r1', 'Rectangle');
    rec_1.set('pos', [a/2+w/2 r/2+(w/4)*sqrt(3)+(a/2)/sqrt(3)]);
    rec_1.set('base', 'center');
    rec_1.set('size', [w r]);
    rec_2 = ucellWP.geom.feature.create('r2', 'Rectangle');
    rec_2.set('pos', [-a*(3/4)+a/2+a*(3/4)-a/2+w/2+a/2 -(w/4)*sqrt(3)+(a/2)/sqrt(3)]);
    rec_2.set('rot', 120);
    rec_2.set('size', [w r]);
    rec_3 = ucellWP.geom.feature.create('r3', 'Rectangle');
    rec_3.set('pos', [-a*(3/4)+a/2+w/2+a*(3/4)-a/2+w/2+a/2 (w/4)*sqrt(3)+(a/2)/sqrt(3)]);
    rec_3.set('rot', 240);
    rec_3.set('size', [w r]);
    
    % implement functions to add fillets to the unit cells 
    addFillets(P,ucellWP)
    
    % set the displacement of the unit cell 
    ucellWP.set('displ', loc);
end

function ucellWP_dup = cellDuplicate(ucellgeom,workPlaneName,workPlaneName_dup,loc)
    % this function duplicates the given unit cell geometry and displace
    % them by loc 
    ucellWP_dup = ucellgeom.feature.duplicate(workPlaneName_dup, workPlaneName);
    ucellWP_dup.set('displ', loc);
end

function addFillets(P,ucellWP)
    r1= P.r1;
    r2 = P.r2;
    a = P.a;
    w = P.w;
    ucellWP.geom.create('fil1', 'Fillet');
    ucellWP.geom.feature('fil1').set('radius', r1);
    ucellWP.geom.feature('fil1').selection('point').set('r1(1)', [3 4]);
    ucellWP.geom.feature('fil1').selection('point').set('r3(1)', [3 4]);
    ucellWP.geom.feature('fil1').selection('point').set('r2(1)', [3 4]);
    centerTriangle = ucellWP.geom.feature.create('pol2', 'Polygon');
    centerTriangle.set('source', 'table');
    centerTriangle.set('table', [a/2 (w/2)*sqrt(3)/2+(a/2)/sqrt(3); w+a/2 (w/2)*sqrt(3)/2+(a/2)/sqrt(3); a/2+w/2 -(w/2)*sqrt(3)/2+(a/2)/sqrt(3); a/2 (w/2)*sqrt(3)/2+(a/2)/sqrt(3)]);
    composit_geom = ucellWP.geom.feature.create('co1', 'Compose');
    composit_geom.set('formula', 'pol1-fil1(1)-fil1(2)-fil1(3)-pol2');
    ucellWP.geom.create('fil2', 'Fillet');
    ucellWP.geom.feature('fil2').set('radius', r2);
    ucellWP.geom.feature('fil2').selection('point').set('co1(1)', [6 10 12]);
end 

function label_list = buildBoomerangCells(P,ucellgeom,xloc,cell_num)
    % buildSnowFlakeRegion: this function builds the boomerang cells that compose one
    % super unit cell
    %  xloc - specify the x location of the unit cell 
    %  cell_num - keeping track of which unit cell in the cavity sequence we are adding 
    a = P.a;
    b = P.b;
    w = P.w;
    cell_num = num2str(cell_num);
    label_list = {}; % list of labels corresponding to each boomerang unit cell
    % adding the first unit cell 
    loc1 = [xloc+a/2-(a/2+w/2) b];
    workPlaneName = ['wp_',cell_num,'_cell_',num2str(1)];
    ucellWP = buildBoomerangCell(P,ucellgeom,workPlaneName, loc1);
    label_list = [label_list,{workPlaneName}];
    % duplicate the second unit cell 
    loc2 = [xloc+a/2-a-(a/2+w/2) b];
    workPlaneName2 = ['wp_',cell_num,'_cell_',num2str(2)];
    ucellWP2 = cellDuplicate(ucellgeom,workPlaneName,workPlaneName2, loc2);
    label_list = [label_list,{workPlaneName2}];
    
    loc3 = [xloc-(a/2+w/2) b+sqrt(3)*a/2];
    workPlaneName3 = ['wp_',cell_num,'_cell_',num2str(3)];
    ucellWP3 = cellDuplicate(ucellgeom,workPlaneName,workPlaneName3, loc3);
    label_list = [label_list,{workPlaneName3}];
    
    loc4 = [xloc+a-(a/2+w/2) b+sqrt(3)*a/2];
    workPlaneName4 = ['wp_',cell_num,'_cell_',num2str(4)];
    ucellWP4 = cellDuplicate(ucellgeom,workPlaneName,workPlaneName4, loc4);
    label_list = [label_list,{workPlaneName4}];
    
    loc5 = [xloc-a-(a/2+w/2) b+sqrt(3)*a/2];
    workPlaneName5 = ['wp_',cell_num,'_cell_',num2str(5)];
    ucellWP5 = cellDuplicate(ucellgeom,workPlaneName,workPlaneName5, loc5);
    label_list = [label_list,{workPlaneName5}];
    
    loc6 = [xloc+a/2-(a/2+w/2) b+a*sqrt(3)/2+sqrt(3)*a/2];
    workPlaneName6 = ['wp_',cell_num,'_cell_',num2str(6)];
    ucellWP6 = cellDuplicate(ucellgeom,workPlaneName,workPlaneName6, loc6);
    label_list = [label_list,{workPlaneName6}];
    
    loc7 = [xloc+a/2-a-(a/2+w/2) b+a*sqrt(3)/2+sqrt(3)*a/2];
    workPlaneName7 = ['wp_',cell_num,'_cell_',num2str(7)];
    ucellWP7 = cellDuplicate(ucellgeom,workPlaneName,workPlaneName7, loc7);
    label_list = [label_list,{workPlaneName7}];
    
end

function workPlaneName = buildLowerCell(ucellgeom,P,cell_num,xloc)
%% this function is used to build lower portion of the boomerang supercell 
cell_num = num2str(cell_num);
workPlaneName = ['lowerCell_',cell_num];
a = P.a;
b = P.b;
d = P.d;
ho = P.ho;
hi = P.hi;
wo = P.wo;
wi = P.wi;
r1 = P.r1;
r2 = P.r2;
ucellWP_lower = ucellgeom.feature.create(workPlaneName, 'WorkPlane');
ucellWP_lower.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', 0);
ucellplane_lower = ucellWP_lower.geom.feature.create('r_base', 'Rectangle');
ucellplane_lower.set('pos', [-a/2 0]);
ucellplane_lower.set('size', [a b]);

rec1 = ucellWP_lower.geom.feature.create('r1', 'Rectangle');
rec1.set('pos', [0 d/2+(ho-hi)/2+hi]);
rec1.set('base', 'center');
rec1.set('size', [wo ho-hi]);

rec2 = ucellWP_lower.geom.feature.create('r2', 'Rectangle');
rec2.set('pos', [-wo/2+(wo-wi)/4 d/2+ho/2]);
rec2.set('base', 'center');
rec2.set('size', [(wo-wi)/2 ho]);

rec3 = ucellWP_lower.geom.feature.create('r3', 'Rectangle');
rec3.set('pos', [wo/2-(wo-wi)/4 d/2+ho/2]);
rec3.set('base', 'center');
rec3.set('size', [(wo-wi)/2 ho]);

compose1 = ucellWP_lower.geom.feature.create('co1_lower', 'Compose');
compose1.set('formula', 'r_base-r1-r2-r3');

% fillet1 = ucellWP_lower.geom.feature.create('fil1', 'Fillet');
% fillet1.set('radius', r1);
% fillet1.selection('point').set('co1(1)', [3 4 6 9 13 14]);
% fillet2 = ucellWP_lower.geom.feature.create('fil2', 'Fillet');
% fillet2.set('radius', r2);
% fillet2.selection('point').set('fil1(1)', [10 13]);

% set the displacement of the unit cell 
ucellWP_lower.set('displ', [xloc,0]);

end




















