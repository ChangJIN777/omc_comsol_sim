% This script creates a single unit cell of triangular nanobeam
% Modified from previous version (Sean Meenehan, 01/26/11) for
% compatibility with COMSOL 5.x.
% Input arguments: model, P
% - model: COMSOL model
% - P: data structure consisting of following geometry parameters
%      a: lattice constant;  w: beam width; th: beam height;
%      hx: hole height; hy: hole width; 
%      hole: hole type ('rect' for rectangular, or 'elps' for elliptical); 
%      optBands: if included, create unit cell for optical bandstructure
%      calculations;
%      airrad: radius of air cylinder for optical bandstructure
%      calculations
% Output: model
% - model: updated COMSOL model with built nanobeam geometry 
%
% Cleaven Chia, 09/06/16

function [model,P] = DrawBlockTet(model,P)

nperiod = P.nperiod;
holeatedge = P.holeatedge; % 1/0 if unit cell terminates in middle of hole/dielectric
nholes = nperiod + holeatedge;

a = P.a;
w = P.w;
th = P.th;
hx = P.hx; % length of block along x-axis
% cy = P.cy;
% tx = P.tx;
hy = P.hy; % width of tether
filRad = P.filRad;
evenz = P.mbevenz;

% lithographic constraints
if hx >= a
    error('block dimensions larger than unit cell dimensions')
end
if hy >= w
    error('tether width larger than block width')
end

%% Create component
ucellcomp = model.modelNode.create('ucell');
ucellcomp.label('Unit cell FEM simulation');
ucelllabel = 'Unit Cell';
ucellname = 'ucell';

ucellgeom = model.geom.create(ucellname, 3);
ucellgeom.label(ucelllabel);

%% Create unit cell with rectangular cross-section
ucellWP = ucellgeom.feature.create('ucellWP', 'WorkPlane');
ucellWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
% ucellplane = ucellWP.geom.feature.create('beamplane', 'Rectangle');
% ucellplane.set('type', 'solid').set('base', 'corner');
% ucellplane.set('pos', [0 0]).set('size', [nperiod*a w/2]);
% ucellgeom.runCurrent;

% solid block
block = ucellWP.geom.create('block', 'Rectangle');
block.set('type', 'solid').set('base', 'corner');
block.set('pos', [(a-hx)/2 0]).set('size', [hx w/2]);
ucellgeom.runCurrent;

% tether
tether = ucellWP.geom.create('tether', 'Rectangle');
tether.set('type', 'solid').set('base', 'center');
tether.set('pos', [a/2 hy/4]).set('size', [a hy/2]);
ucellgeom.runCurrent;

% compose
ucellComp = ucellWP.geom.feature.create('ucellComp', 'Compose');
ucellComp.selection('input').set({'block' 'tether'});
ucellComp.set('formula', 'block + tether').set('intbnd', false);
ucellgeom.runCurrent;

%% If fillet: choose vertices for fillet
if filRad > 0
    
    
    % create selection
    delta = 10e-9; 
    ucellFilSel = ucellWP.geom.create('ucellFilSel', 'BoxSelection');
    ucellFilSel.set('xmin', (a-hx)/2-delta).set('xmax', (a+hx)/2+delta);
    ucellFilSel.set('ymin', hy/2-delta).set('ymax', w/2+delta);
    ucellFilSel.set('condition', 'inside').set('entitydim', 0);
    ucellgeom.runCurrent;
    
    ucellFillet = ucellWP.geom.create('ucellFillet', 'Fillet');
    ucellFillet.set('radius', filRad);
    ucellFillet.selection('point').named('ucellFilSel');
    ucellgeom.runCurrent;
end

%% Extrude unit cell
ucellPlaneExt = ucellgeom.feature.create('ucellPlaneExt', 'Extrude');
ucellPlaneExt.set('distance', th);
ucellgeom.runCurrent;
ucellName = 'ucellPlaneExt';

%% Symmetries in y and z
if evenz
    % create symmetry block
    symZ = abs(evenz);
    ucsymZWP = ucellgeom.feature.create('ucsymZWP', 'WorkPlane');
    ucsymZWP.set('planetype', 'quick').set('quickplane', 'xy').set('quickz', -th/2);
    ucsymZPlane = ucsymZWP.geom.feature.create('ucsymZPlane', 'Rectangle');
    ucsymZPlane.set('type', 'solid').set('base', 'corner');
    ucsymZPlane.set('pos', [0 0]).set('size', [a w/2]);

    % extrude symmetry block
    ucellgeom.runCurrent;
    ucsymZPlaneExt = ucellgeom.feature.create('ucsymZPlaneExt', 'Extrude');
    ucsymZPlaneExt.set('distance', th/2^symZ);
    
    % compose: unit cell - symmetry block
    ucsymZComp = ucellgeom.feature.create('ucsymZComp', 'Compose');
    ucsymZComp.selection('input').set(ucellName);
    ucsymZComp.selection('input').set('ucsymZPlaneExt');
    ucsymZComp.set('formula', [ucellName,' - ucsymZPlaneExt']);
    ucellgeom.run;
    ucellName = 'ucsymZComp';
end

%% Create selections for boundary conditions
delta = 10e-9;
% beam x-symmetry plane
XsymSel = ucellgeom.create('XsymSel', 'BoxSelection');
XsymSel.set('xmin', -delta).set('xmax', delta);
XsymSel.set('ymin', -delta).set('ymax', w/2+delta);
XsymSel.set('zmin', -th/2*(1-evenz)-delta).set('zmax', th/2+delta);
XsymSel.set('entitydim', 2).set('condition', 'allvertices');
ucellgeom.runCurrent;
inds = model.selection([ucellname,'_XsymSel']).inputEntities();
P.bndSel.Xsym = inds';

% beam x-end plane
XendSel = ucellgeom.create('XendSel', 'BoxSelection');
XendSel.set('xmin', nperiod*a-delta).set('xmax', nperiod*a+delta);
XendSel.set('ymin', -delta).set('ymax', w/2+delta);
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

% beam z-symmetry plane
ZsymSel = ucellgeom.create('ZsymSel', 'BoxSelection');
ZsymSel.set('xmin', -delta).set('xmax', nperiod*a+delta);
ZsymSel.set('ymin', -delta).set('ymax', w/2+delta);
ZsymSel.set('zmin', -delta).set('zmax', delta);
ZsymSel.set('entitydim', 2).set('condition', 'allvertices'); % only want beam surface, exclude air holes
ucellgeom.runCurrent;
inds = model.selection([ucellname,'_ZsymSel']).inputEntities();
P.bndSel.Zsym = inds';
