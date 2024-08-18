function [model,P] = buildBoomerangUnitCellStrip(model,P)
%
% buildBoomerangUnitCell.m
%
% Model exported on Jul 21 2024, 16:57 by COMSOL 6.1.0.357.

%% read the input parameters 
a = P.a;        % lattice constant 
w = P.w;        % the width of the hole 
th = P.th;        % the height of the hole
r = P.r;        % the height of the hole 
h = P.h;        % the height of the hole in the lower portion
d = P.d;        % the width of the hole in the lower portion
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


%% duplicate the unit cells 
ucellWP.set('displ', [-(a/2+w/2) a*sqrt(3)/4]);

ucellWP2 = ucellgeom.feature.duplicate('wp2', 'wp1');
ucellWP2.set('displ', [a-(a/2+w/2) a*sqrt(3)/4]);

ucellWP3 = ucellgeom.feature.duplicate('wp3', 'wp1');
ucellWP3.set('displ', [-a-(a/2+w/2) a*sqrt(3)/4]);

ucellWP4 = ucellgeom.feature.duplicate('wp4', 'wp1');
ucellWP4.set('displ', [a/2-(a/2+w/2) a*sqrt(3)/4+a*sqrt(3)/2]);

ucellWP5 = ucellgeom.feature.duplicate('wp5', 'wp1');
ucellWP5.set('displ', [a/2-a-(a/2+w/2) a*sqrt(3)/4+a*sqrt(3)/2]);

ucellWP6 = ucellgeom.feature.duplicate('wp6', 'wp1');
ucellWP6.set('displ', [-(a/2+w/2) a*sqrt(3)/4+sqrt(3)*a]);

ucellWP7 = ucellgeom.feature.duplicate('wp7', 'wp1');
ucellWP7.set('displ', [a-(a/2+w/2) a*sqrt(3)/4+sqrt(3)*a]);

ucellWP8 = ucellgeom.feature.duplicate('wp8', 'wp1');
ucellWP8.set('displ', [-a-(a/2+w/2) a*sqrt(3)/4+sqrt(3)*a]);

ucellWP9 = ucellgeom.feature.duplicate('wp9', 'wp1');
ucellWP9.set('displ', [a/2-(a/2+w/2) a*sqrt(3)/4+a*sqrt(3)/2+sqrt(3)*a]);

ucellWP10 = ucellgeom.feature.duplicate('wp10', 'wp1');
ucellWP10.set('displ', [a/2-a-(a/2+w/2) a*sqrt(3)/4+a*sqrt(3)/2+sqrt(3)*a]);

% deactivate unused planes 
ucellWP.active(false);
ucellWP2.active(false);
ucellWP3.active(false);


%% create unit cell with boomerang geometry in the lower cavity
ucellWP_lower = ucellgeom.create('wp_lower', 'WorkPlane');
ucellWP_lower.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
ucellplane_lower = ucellWP_lower.geom.feature.create('r_base', 'Rectangle');
ucellplane_lower.set('base', 'center');
ucellplane_lower.set('pos', [0 -sqrt(3)*a*(3/8)+a/2]);
ucellplane_lower.set('size', [a a*sqrt(3)*(3/4)]);
rec1 = ucellWP_lower.geom.feature.create('r1', 'Rectangle');
rec1.set('pos', [0 r/2+sqrt(3)*w/4]);
rec1.set('base', 'center');
rec1.set('size', [w r]);
rec2 = ucellWP_lower.geom.feature.create('r2', 'Rectangle');
rec2.set('pos', [0 -sqrt(3)*w/2+sqrt(3)*w/4]);
rec2.set('rot', 120);
rec2.set('size', [w r]);
rec3 = ucellWP_lower.geom.feature.create('r3', 'Rectangle');
rec3.set('pos', [w/2 sqrt(3)*w/4]);
rec3.set('rot', 240);
rec3.set('size', [w r]);
rec4 = ucellWP_lower.geom.feature.create('r4', 'Rectangle');
rec4.set('pos', [w/2+sqrt(3)*r/2 -r/2+sqrt(3)*w/4]);
rec4.set('rot', 180);
rec4.set('size', [d h]);
rec5 = ucellWP_lower.geom.feature.create('r5', 'Rectangle');
rec5.set('pos', [-w/2-sqrt(3)*r/2+a/2+d-a/2 a/2-r/2+sqrt(3)*w/4-a/2]);
rec5.set('rot', 180);
rec5.set('size', [d h]);
pol1 = ucellWP_lower.geom.feature.create('pol1', 'Polygon');
pol1.set('tableconstr', {'off' 'off'});
pol1.set('source', 'table');
pol1.set('table', [-w/2 sqrt(3)*w/4; w/2 sqrt(3)*w/4; 0 -sqrt(3)*w/4; -w/2 sqrt(3)*w/4]);
compose1 = ucellWP_lower.geom.feature.create('co1', 'Compose');
compose1.set('formula', 'r_base-r1-r2-r3-r4-r5-pol1');
fillet1 = ucellWP_lower.geom.feature.create('fil1', 'Fillet');
fillet1.set('radius', r1);
fillet1.selection('point').set('co1(1)', [3 5 8 11 12 14]);
fillet2 = ucellWP_lower.geom.feature.create('fil2', 'Fillet');
fillet2.set('radius', r2);
fillet2.selection('point').set('fil1(1)', [4 8 9 12 14 17 21]);

ucellWP_lower.set('displ', [0 a*sqrt(3)*(3/4)-a/2]);

%% the strip work plane 
ucellWP11 = ucellgeom.feature.create('wp11', 'WorkPlane');
ucellWP11.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
rec_base = ucellWP11.geom.feature.create('r_base', 'Rectangle');
rec_base.set('size',[a sqrt(3)*2*a]);
rec_base.set('pos',[-a/2 0]);

ucellWP12 = ucellgeom.feature.create('wp12', 'WorkPlane');
ucellWP12.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
rec_base2 = ucellWP12.geom.feature.create('r_base2', 'Rectangle');
rec_base2.set('size',[a sqrt(3)*a/4]);
rec_base2.set('pos',[-a/2 0]);
ucellgeom.runAll;


%% create the intersect
compose2 = ucellgeom.feature.create('symZComp1', 'Compose');
compose2.selection('input').set({'wp10' 'wp11' 'wp4' 'wp5' 'wp6' 'wp7' 'wp8' ...
'wp9' 'wp10' 'wp11' 'wp12','wp_lower'});
compose2.set('formula', 'wp11*(wp_lower+wp4+wp5+wp6+wp7+wp8+wp9+wp10+wp12)');
ucellgeom.runAll;

% extrude 
extrude = ucellgeom.feature.create('ext1', 'Extrude');
extrude.set('extrudefrom', 'faces');
extrude.setIndex('distance', th, 0);
extrude.selection('inputface').set('symZComp1(1)',[1 2 3 4 5 6 7 8]);
ucellgeom.runAll;

%% Symmetry in z 
if abs(P.mbevenz)
    symZth = th/2;
    symW = a;
    % create symmetry block
    symZWP = ucellgeom.feature.create('symZWP', 'WorkPlane');
    symZWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -symZth);
%     symZPlane = symZWP.geom.feature.create('pol2', 'Polygon');
%     symZPlane.set('source', 'table');
%     symZPlane.set('table', [0 0; symW/2 (symW/2)*sqrt(3); symW*(3/2) (symW/2)*sqrt(3); symW 0; 0 0]);
    symZPlane = symZWP.geom.feature.create('pol2', 'Rectangle');
    symZPlane.set('size',[a sqrt(3)*2*a]);
    symZPlane.set('pos',[-a/2 0]);
    
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
    ZsymSel.set('ymin', -a/2-delta).set('ymax', a/2+delta);
    ZsymSel.set('zmin', -delta).set('zmax', delta);
    ZsymSel.set('entitydim', 2).set('condition', 'allvertices'); % only want beam surface, exclude air holes
    ucellgeom.runCurrent;
    inds = model.selection([ucellname,'_ZsymSel']).inputEntities();
    P.bndSel.Zsym = inds';

end
%% Making selections (manual)
mphgeom(model);
P.xEnd1 =  bndindex(ucellgeom, [-a/2 0 0], [-1 0 0]);
P.xEnd2 = bndindex(ucellgeom, [a/2 0 0], [1 0 0]);

P.yEnd1 =  bndindex(ucellgeom, [0 0 0], [0 1 0]);
P.yEnd2 = bndindex(ucellgeom, [0 sqrt(3)*2*a 0], [0 1 0]);
% Note that this will return no indices if there is no boundary at z=0
P.zEnd = bndindex(ucellgeom, [0 0 0], [0 0 1]);
% % debugging
% mphplot(model);
disp(P) % debugging
out = model;
end



