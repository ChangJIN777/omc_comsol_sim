function [model,P] = buildRibUnitCell(model,P)
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
wi = P.wi;
hi = P.hi;
ho = P.ho;
wo = P.wo;
ai = P.ai;
th = P.th;
P.beamLen = nperiod*a;
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

%% create unit cell with rectangular cross-section and rib geoemtry
ucellWP = ucellgeom.feature.create('wp1', 'WorkPlane');
ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
ucellplane = ucellWP.geom.feature.create('ucellplane', 'Rectangle');
ucellplane.set('type', 'solid').set('base', 'corner');
% ucellplane.set('pos', [-P.beamLen/2 -a/2]).set('size', [P.beamLen a]); % no y symmetry 
ucellplane.set('pos', [-P.beamLen/2 0]).set('size', [P.beamLen w/2]); % with y symmetry 
ucellgeom.runCurrent;

% create holes 
holeList = {};
holeFormula = [];
for k = 1:nholes
    holeID = ['hole_' num2str(k)];
    holeID_upper = [holeID '_upper'];
    % holeID_lower = [holeID '_lower'];
    hole_upper = ucellWP.geom.feature.create(holeID_upper, 'Polygon');
    % hole_lower = ucellWP.geom.feature.create(holeID_lower, 'Polygon');
    hole_upper.set('source', 'table');
    % hole_lower.set('source', 'table');
    hole_upper.set('table', [-wo/2 ai/2; -wi/2 ai/2; -wi/2 ai/2+hi; wi/2 ai/2+hi; wi/2 ai/2; wo/2 ai/2; wo/2 ai/2+ho; -wo/2 ai/2+ho; -wo/2 ai/2]);
    % hole_lower.set('table', [-wo -ai; -wi -ai; -wi -(ai+hi); wi -(ai+hi); wi -ai; wo -ai; wo -(ai+ho); -wo -(ai+ho); -wo -ai]);
    holeList = [holeList,holeID_upper];
    holeFormula = [holeFormula,'-',holeID,'_upper'];
end
ucellgeom.runCurrent;

% compose workplane 
ucellComp = ucellWP.geom.feature.create('ucellComp','Compose');
ucellComp.selection('input').set(['ucellplane',holeList]);
ucellComp.set('formula', ['ucellplane',holeFormula]);
ucellgeom.runCurrent;

% extrude
finBeamTag = 'ucellholes';
ucellHoles = ucellgeom.feature.create(finBeamTag,'Extrude');
ucellHoles.set('distance',th);
ucellgeom.runCurrent;
displayBeamStr = 'Unit cell, rib geometry, rectangular cross-section';
% 
% debugging 
figure; 
mphgeom(model);
ucellgeom.runAll;
%% Symmetry in z 
if abs(P.mbevenz)
    symZth = th/2;
    symW = w;
    % create symmetry block
    symZWP = ucellgeom.feature.create('symZWP', 'WorkPlane');
    symZWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -symZth);
    symZPlane = symZWP.geom.feature.create('symZPlane', 'Rectangle');
    symZPlane.set('type', 'solid').set('base', 'corner');
    symZPlane.set('pos', [-P.beamLen/2 -symW/2]).set('size', [P.beamLen symW]);

    % extrude symmetry block
    ucellgeom.runCurrent;
    symZPlaneExt = ucellgeom.feature.create('symZPlaneExt', 'Extrude');
    symZPlaneExt.set('distance', symZth);

    % compose: unit cell - symmetry block
    symZComp = ucellgeom.feature.create('symZComp', 'Compose');
    symZComp.selection('input').set(finBeamTag);
    symZComp.selection('input').set('symZPlaneExt');
    symZComp.set('formula', [finBeamTag,' - symZPlaneExt']);
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

%% Create selections for boundary conditions
delta = 10e-9;
% beam x-symmetry plane
XsymSel = ucellgeom.create('XsymSel', 'BoxSelection');
XsymSel.set('xmin', -delta).set('xmax', delta);
XsymSel.set('ymin', -delta).set('ymax', a/2+delta);
XsymSel.set('zmin', -th/2*(1-evenz)-delta).set('zmax', th/2+delta);
XsymSel.set('entitydim', 2).set('condition', 'allvertices');
ucellgeom.runCurrent;
inds = model.selection([ucellname,'_XsymSel']).inputEntities();
P.bndSel.Xsym = inds';

% beam x-end plane
XendSel = ucellgeom.create('XendSel', 'BoxSelection');
XendSel.set('xmin', nperiod*a-delta).set('xmax', nperiod*a+delta);
XendSel.set('ymin', -delta).set('ymax', a/2+delta);
XendSel.set('zmin', -th/2*(1-evenz)-delta).set('zmax', th/2+delta);
XendSel.set('entitydim', 2).set('condition', 'allvertices');
ucellgeom.runCurrent;
inds = model.selection([ucellname,'_XendSel']).inputEntities();
P.bndSel.Xend = inds';

% beam y-symmetry plane 
YsymSel = ucellgeom.create('YsymSel', 'BoxSelection');
YsymSel.set('xmin', -delta).set('xmax', nperiod*a+delta);
YsymSel.set('ymin', -delta).set('ymax', delta);
YsymSel.set('zmin', -th/2*(1-evenz)-delta).set('zmax', th/2+delta);
YsymSel.set('entitydim', 2).set('condition', 'allvertices');
ucellgeom.runCurrent;
inds = model.selection([ucellname,'_YsymSel']).inputEntities();
P.bndSel.Ysym = inds';

% debugging 
mphplot(ucellgeom);

%% Making selections
P.xEnd1 =  bndindex(ucellgeom, [-a/2 0 0], [1 0 0]);
P.xEnd2 = bndindex(ucellgeom, [ a/2 0 0], [1 0 0]);

P.yEnd1 =  bndindex(ucellgeom, [0 0 0], [0 1 0]);
P.yEnd2 = bndindex(ucellgeom, [0 w/2 0], [0 1 0]);
% Note that this will return no indices if there is no boundary at z=0
P.zEnd = bndindex(ucellgeom, [0 0 0], [0 0 1]);
disp(P.xEnd1)
out = model;
