function [model,P] = buildCrossUnitCell(model,P)
% build the unit cell geometry for rib OMC
% buildRibUnitCell.m
%
% Model exported on Jul 4 2024, 16:12 by COMSOL 6.1.0.357.

nperiod = P.nperiod;
holeatedge = P.holeatedge; % 1/0 if unit cell terminates in middle of hole/dielectric
nholes = nperiod + holeatedge;

%% read the parameters from the input 
a = P.a;
w = P.w;
h = P.h;
th = P.th;
r1 = P.r1;
r2 = P.r2;
evenz = P.mbevenz;

%% Create component
ucellcomp = model.modelNode.create('comp1');
% model.component.create('comp1', true);
ucellcomp.label('Unit cell FEM simulation');
ucelllabel = 'Unit Cell';
ucellname = 'geom1';

% ucellgeom = model.component('comp1').geom.create(ucellname, 3);
ucellgeom = model.geom.create(ucellname, 3);
ucellgeom.label(ucelllabel);

%% create unit cell with rectangular cross-section and cross 
ucellWP = ucellgeom.feature.create('wp1', 'WorkPlane');
ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);

ucellPlane = ucellWP.geom.feature.create('r_ucell', 'Rectangle');
ucellPlane.set('pos', [0 0]);
ucellPlane.set('base', 'center');
ucellPlane.set('size', [a a]);

rec_1 = ucellWP.geom.feature.create('r1', 'Rectangle');
rec_1.set('pos', [0 0]);
rec_1.set('base', 'center');
rec_1.set('size', [h w]);

rec_2 = ucellWP.geom.feature.create('r2', 'Rectangle');
rec_2.set('pos', [0 0]);
rec_2.set('base', 'center');
rec_2.set('size', [w h]);

ucellgeom.runCurrent;

composit_geom = ucellWP.geom.feature.create('co1', 'Compose');
composit_geom.set('formula', 'r_ucell-r1-r2');

ucellgeom.runCurrent;

% add the fillet to the structures 
selection_width = 10e-9;
addFillet(P,ucellWP,selection_width);
ucellgeom.runCurrent;

% extrude
finBeamTag = 'ext1';
ucellHoles = ucellgeom.feature.create(finBeamTag,'Extrude');
ucellHoles.set('distance',th);
ucellgeom.runCurrent;
% 
% debugging 
figure; 
mphgeom(model);
ucellgeom.runAll;
%% Symmetry in z 
if abs(P.mbevenz)
    symZth = th/2;
    symW = a;
    % create symmetry block
    symZWP = ucellgeom.feature.create('symZWP', 'WorkPlane');
    symZWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -symZth);
    symZPlane = symZWP.geom.feature.create('symZPlane', 'Rectangle');
    symZPlane.set('type', 'solid').set('base', 'center');
    symZPlane.set('pos', [0 0]).set('size', [P.a symW]);

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
    ZsymSel.set('xmin', -nperiod*a/2-delta).set('xmax', nperiod*a/2+delta);
    ZsymSel.set('ymin', -nperiod*a/2-delta).set('ymax', nperiod*a/2+delta);
    ZsymSel.set('zmin', -delta).set('zmax', delta);
    ZsymSel.set('entitydim', 2).set('condition', 'allvertices'); % only want beam surface, exclude air holes
    ucellgeom.runCurrent;
    inds = model.selection([ucellname,'_ZsymSel']).inputEntities();
    P.bndSel.Zsym = inds';

end

%% Making selections
P.xEnd1 =  bndindex(ucellgeom, [-a/2 0 0], [1 0 0]);
P.xEnd2 = bndindex(ucellgeom, [ a/2 0 0], [1 0 0]);

P.yEnd1 =  bndindex(ucellgeom, [0 -a/2 0], [0 1 0]);
P.yEnd2 = bndindex(ucellgeom, [0 a/2 0], [0 1 0]);
% Note that this will return no indices if there is no boundary at z=0
P.zEnd = bndindex(ucellgeom, [h/4 h/4 -th/2], [0 0 1]);
out = model;
end 

function addFillet(P,ucellWP,selection_width)
    w = P.w;
    h = P.h;
    r1 = P.r1;
    r2 = P.r2;
    disksel1_label = sprintf('h_disksel1');
    disksel2_label = sprintf('h_disksel2');
    fillet1_label = sprintf('h_fil1');
    fillet2_label = sprintf('h_fil2');
    % hole 
    disksel1 = ucellWP.geom.feature.create(disksel1_label, 'DiskSelection');
    disksel1.set('entitydim', 0);
    disksel1.set('posx', 0);
    disksel1.set('posy', 0);
    disksel1.set('r', w*(3/2));
    disksel1.set('rin', w/2);
    disksel1.set('condition', 'allvertices');
    fil1 = ucellWP.geom.feature.create(fillet1_label, 'Fillet');
    fil1.set('radius', r1);
    fil1.selection('point').named(disksel1_label);
    disksel2 = ucellWP.geom.feature.create(disksel2_label, 'DiskSelection');
    disksel2.set('entitydim', 0);
    disksel2.set('posy', 0);
    disksel2.set('posx', 0);
    disksel2.set('r', h/2+selection_width/2);
    disksel2.set('rin', h/2-selection_width/2);
    disksel2.set('condition', 'allvertices');
    fil2 = ucellWP.geom.feature.create(fillet2_label, 'Fillet');
    fil2.set('radius', r2);
    fil2.selection('point').named(disksel2_label);
end