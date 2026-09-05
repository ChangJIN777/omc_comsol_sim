function out = model
%
% BoomerangWithAirDisk.m
%
% Model exported on Aug 13 2026, 12:02 by COMSOL 6.3.0.290.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/Users/changjin/Documents/GitHub/omc_comsol_sim/comsol_templates');

model.component.create('comp1', false);

model.component('comp1').geom.create('geom1', 3);

model.component('comp1').label('Unit cell FEM simulation');

model.component('comp1').geom('geom1').label('Unit Cell');
model.component('comp1').geom('geom1').geomRep('cadps');
model.component('comp1').geom('geom1').designBooleans(false);
model.component('comp1').geom('geom1').create('wp1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp1').set('quickz', -1.1E-7);
model.component('comp1').geom('geom1').feature('wp1').geom.create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').label('Base plane');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').set('table', [0 0; 2.4E-7 4.156921938165305E-7; 7.199999999999999E-7 4.156921938165305E-7; 4.8E-7 0; 0 0]);
model.component('comp1').geom('geom1').feature('wp1').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('size', [1.4E-7 1.77E-7]);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('pos', [3.5999999999999994E-7 2.9634609690826523E-7]);
model.component('comp1').geom('geom1').feature('wp1').geom.create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').set('size', [1.4E-7 1.77E-7]);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').set('pos', [2.8335675176507715E-7 1.6359609690826526E-7]);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').set('rot', 120);
model.component('comp1').geom('geom1').feature('wp1').geom.create('r3', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r3').set('size', [1.4E-7 1.77E-7]);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r3').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r3').set('pos', [4.3664324823492274E-7 1.6359609690826526E-7]);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r3').set('rot', 240);
model.component('comp1').geom('geom1').feature('wp1').geom.create('co1', 'Compose');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('co1').set('formula', 'pol1-r1-r2-r3');
model.component('comp1').geom('geom1').feature('wp1').geom.create('h_disksel1', 'DiskSelection');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel1').set('entitydim', 0);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel1').set('posx', 3.5999999999999994E-7);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel1').set('posy', 2.0784609690826525E-7);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel1').set('r', 1.4E-7);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel1').set('rin', 4.949747468305833E-8);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel1').set('condition', 'allvertices');
model.component('comp1').geom('geom1').feature('wp1').geom.create('h_fil1', 'Fillet');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_fil1').set('radius', 1.0E-8);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_fil1').selection('point').named('h_disksel1');
model.component('comp1').geom('geom1').feature('wp1').geom.create('h_disksel2', 'DiskSelection');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel2').set('entitydim', 0);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel2').set('posx', 3.5999999999999994E-7);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel2').set('posy', 2.0784609690826525E-7);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel2').set('r', 2.02E-7);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel2').set('rin', 1.52E-7);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel2').set('condition', 'allvertices');
model.component('comp1').geom('geom1').feature('wp1').geom.create('h_fil2', 'Fillet');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_fil2').set('radius', 1.0E-8);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_fil2').selection('point').named('h_disksel2');
model.component('comp1').geom('geom1').create('ext1', 'Extrude');
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', '2.2E-7', 0);
model.component('comp1').geom('geom1').feature('ext1').selection('input').set({'wp1'});
model.component('comp1').geom('geom1').create('symZWP', 'WorkPlane');
model.component('comp1').geom('geom1').feature('symZWP').set('quickz', -1.1E-7);
model.component('comp1').geom('geom1').feature('symZWP').geom.create('pol2', 'Polygon');
model.component('comp1').geom('geom1').feature('symZWP').geom.feature('pol2').set('source', 'table');
model.component('comp1').geom('geom1').feature('symZWP').geom.feature('pol2').set('table', [0 0; 2.4E-7 4.156921938165305E-7; 7.199999999999999E-7 4.156921938165305E-7; 4.8E-7 0; 0 0]);
model.component('comp1').geom('geom1').create('symZPlaneExt', 'Extrude');
model.component('comp1').geom('geom1').feature('symZPlaneExt').setIndex('distance', '1.1E-7', 0);
model.component('comp1').geom('geom1').feature('symZPlaneExt').selection('input').set({'symZWP'});
model.component('comp1').geom('geom1').create('symZComp', 'Compose');
model.component('comp1').geom('geom1').feature('symZComp').set('formula', 'ext1 - symZPlaneExt');
model.component('comp1').geom('geom1').create('symZWP1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('symZWP1').label('airDisk');
model.component('comp1').geom('geom1').feature('symZWP1').geom.create('pol2', 'Polygon');
model.component('comp1').geom('geom1').feature('symZWP1').geom.feature('pol2').set('source', 'table');
model.component('comp1').geom('geom1').feature('symZWP1').geom.feature('pol2').set('table', [0 0; 2.4E-7 4.156921938165305E-7; 7.199999999999999E-7 4.156921938165305E-7; 4.8E-7 0; 0 0]);
model.component('comp1').geom('geom1').create('ext2', 'Extrude');
model.component('comp1').geom('geom1').feature('ext2').setIndex('distance', '5000[um]', 0);
model.component('comp1').geom('geom1').feature('ext2').selection('input').set({'symZWP1'});
model.component('comp1').geom('geom1').create('ZsymSel', 'BoxSelection');
model.component('comp1').geom('geom1').feature('ZsymSel').set('entitydim', 2);
model.component('comp1').geom('geom1').feature('ZsymSel').set('xmin', '-inf');
model.component('comp1').geom('geom1').feature('ZsymSel').set('xmax', 'inf');
model.component('comp1').geom('geom1').feature('ZsymSel').set('ymin', '-inf');
model.component('comp1').geom('geom1').feature('ZsymSel').set('ymax', 'inf');
model.component('comp1').geom('geom1').feature('ZsymSel').set('zmin', -1.0E-8);
model.component('comp1').geom('geom1').feature('ZsymSel').set('zmax', 1.0E-8);
model.component('comp1').geom('geom1').feature('ZsymSel').set('condition', 'allvertices');
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.label('BoomerangWithAirDisk.mph');

out = model;
