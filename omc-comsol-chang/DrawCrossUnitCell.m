
function [model,P] = DrawCrossUnitCell(model,P)

%DRAWCROSS Summary of this function goes here
%   Detailed explanation goes here
a = P.a;
hc = P.hc;
wc = P.wc;
th = P.th;
r1 = P.r1;
r2 = P.r2;
abssym = abs(P.mbevenz);

%% Create component
ucellcomp = model.modelNode.create('comp1');
% model.component.create('comp1', true);
ucellcomp.label('Unit cell FEM simulation');
ucelllabel = 'Unit Cell';
ucellname = 'geom1';
ucellgeom = model.geom.create(ucellname, 3);
ucellgeom.label(ucelllabel);

wp = ucellgeom.create('wp1', 'WorkPlane');
wp.set('unite', true);
if ~abssym
    wp.set('quickz',-th/2);
end

wp.geom.create('sq1', 'Square');
wp.geom.feature('sq1').set('base', 'center');
wp.geom.feature('sq1').set('size', a);
wp.geom.create('r1', 'Rectangle');
wp.geom.feature('r1').set('base', 'center');
wp.geom.feature('r1').set('size', [hc wc]);
wp.geom.create('r2', 'Rectangle');
wp.geom.feature('r2').set('base', 'center');
wp.geom.feature('r2').set('size', [wc hc]);
wp.geom.create('dif1', 'Difference');
wp.geom.feature('dif1').selection('input').set({'sq1'});
wp.geom.feature('dif1').selection('input2').set({'r1' 'r2'});
wp.geom.create('fil1', 'Fillet');
wp.geom.feature('fil1').set('radius', r1);
wp.geom.feature('fil1').selection('point').set('dif1(1)', [3 4 5 8 9 12 13 14]);
wp.geom.create('fil2', 'Fillet');
wp.geom.feature('fil2').set('radius', r2);
wp.geom.feature('fil2').selection('point').set('fil1(1)', [8 9 16 17]);

% extrude
finBeamTag = 'ucellholes';
ucellHoles = ucellgeom.feature.create(finBeamTag, 'Extrude');
ucellHoles.selection('input').set({'wp1'});
if abssym
ucellHoles.setIndex('distance', th/2, 0);
else
ucellHoles.setIndex('distance', th, 0);
end
ucellgeom.run;
ucellgeom.run('fin');

%% Symmetry in z 
if abs(P.mbevenz)
    symZth = th/2;
    symW = a;
    % create symmetry block
    symZWP = ucellgeom.feature.create('symZWP', 'WorkPlane');
    symZWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -symZth);
    symZPlane = symZWP.geom.feature.create('symZPlane', 'Rectangle');
    symZPlane.set('type', 'solid').set('base', 'corner');
    symZPlane.set('pos', [-P.a/2 -symW/2]).set('size', [P.a symW]);

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
    ZsymSel.set('xmin', -a/2-delta).set('xmax', a/2+delta);
    ZsymSel.set('ymin', -a/2-delta).set('ymax', a/2+delta);
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
P.zEnd = bndindex(ucellgeom, [0 0 0], [0 0 1]);
disp(P) % debugging

end