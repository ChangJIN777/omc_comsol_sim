function out = model
%
% testblocktet.m
%
% Model exported on Apr 29 2019, 18:33 by COMSOL 5.3.1.348.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('G:\My Drive\OMC\COMSOL\UnitCell1D');

model.component.create('comp');

model.component('comp').label('Unit cell FEM simulation');

model.component('comp').geom.create('ucell', 3);
model.component('comp').geom('ucell').label('Unit Cell');
model.component('comp').geom('ucell').feature.create('ucellWP', 'WorkPlane');
model.component('comp').geom('ucell').feature('ucellWP').set('planetype', 'quick');
model.component('comp').geom('ucell').feature('ucellWP').set('quickplane', 'xy');
model.component('comp').geom('ucell').feature('ucellWP').set('quickz', -2.0E-7);
model.component('comp').geom('ucell').feature('ucellWP').geom.create('block', 'Rectangle');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('block').set('type', 'solid');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('block').set('base', 'corner');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('block').set('pos', [1.5E-7 0]);
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('block').set('size', [3.0E-7 4.0E-7]);
model.component('comp').geom('ucell').run('ucellWP');
model.component('comp').geom('ucell').feature('ucellWP').geom.create('tether', 'Rectangle');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('tether').set('type', 'solid');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('tether').set('base', 'center');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('tether').set('pos', [3.0E-7 5.0E-8]);
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('tether').set('size', [6.0E-7 1.0E-7]);
model.component('comp').geom('ucell').run('ucellWP');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature.create('ucellComp', 'Compose');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('ucellComp').selection('input').set({'block' 'tether'});
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('ucellComp').set('formula', 'block + tether');
model.component('comp').geom('ucell').run('ucellWP');
model.component('comp').geom('ucell').feature.create('ucellPlaneExt', 'Extrude');
model.component('comp').geom('ucell').feature('ucellPlaneExt').set('distance', 4.0E-7);
model.component('comp').geom('ucell').run('ucellPlaneExt');
model.component('comp').geom('ucell').feature.create('ucsymZWP', 'WorkPlane');
model.component('comp').geom('ucell').feature('ucsymZWP').set('planetype', 'quick');
model.component('comp').geom('ucell').feature('ucsymZWP').set('quickplane', 'xy');
model.component('comp').geom('ucell').feature('ucsymZWP').set('quickz', -2.0E-7);
model.component('comp').geom('ucell').feature('ucsymZWP').geom.feature.create('ucsymZPlane', 'Rectangle');
model.component('comp').geom('ucell').feature('ucsymZWP').geom.feature('ucsymZPlane').set('type', 'solid');
model.component('comp').geom('ucell').feature('ucsymZWP').geom.feature('ucsymZPlane').set('base', 'corner');
model.component('comp').geom('ucell').feature('ucsymZWP').geom.feature('ucsymZPlane').set('pos', [0 0]);
model.component('comp').geom('ucell').feature('ucsymZWP').geom.feature('ucsymZPlane').set('size', [6.0E-7 4.0E-7]);
model.component('comp').geom('ucell').run('ucsymZWP');
model.component('comp').geom('ucell').feature.create('ucsymZPlaneExt', 'Extrude');
model.component('comp').geom('ucell').feature('ucsymZPlaneExt').set('distance', 2.0E-7);
model.component('comp').geom('ucell').feature.create('ucsymZComp', 'Compose');
model.component('comp').geom('ucell').feature('ucsymZComp').selection('input').set({'ucsymZPlaneExt'});
model.component('comp').geom('ucell').feature('ucsymZComp').set('formula', 'ucellPlaneExt - ucsymZPlaneExt');
model.component('comp').geom('ucell').run;

model.label('testblocktet.mph');

model.component('comp').geom('ucell').run('ucsymZComp');
model.component('comp').geom('ucell').create('boxsel1', 'BoxSelection');
model.component('comp').geom('ucell').feature('boxsel1').set('xmin', '1.5e-7');
model.component('comp').geom('ucell').feature('boxsel1').set('xmax', '4.5e-7');
model.component('comp').geom('ucell').feature('boxsel1').set('ymin', 0);
model.component('comp').geom('ucell').feature('boxsel1').set('ymax', '4e-7');
model.component('comp').geom('ucell').feature('boxsel1').set('zmin', 0);
model.component('comp').geom('ucell').feature('boxsel1').set('zmax', '2e-7');
model.component('comp').geom('ucell').run('boxsel1');
model.component('comp').geom('ucell').feature.move('boxsel1', 4);
model.component('comp').geom('ucell').feature.move('boxsel1', 3);
model.component('comp').geom('ucell').feature.move('boxsel1', 2);
model.component('comp').geom('ucell').feature.move('boxsel1', 1);
model.component('comp').geom('ucell').feature.remove('boxsel1');
model.component('comp').geom('ucell').feature('ucellWP').geom.run('ucellComp');
model.component('comp').geom('ucell').feature('ucellWP').geom.create('boxsel1', 'BoxSelection');
model.component('comp').geom('ucell').feature('ucellWP').geom.run('ucellComp');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('boxsel1').set('entitydim', 0);
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('boxsel1').set('xmin', '1.5e-7');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('boxsel1').set('xmax', '4.5e-7');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('boxsel1').set('ymin', 0);
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('boxsel1').set('ymax', '4e-7');
model.component('comp').geom('ucell').feature('ucellWP').geom.run('boxsel1');
model.component('comp').geom('ucell').feature('ucellWP').geom.run('boxsel1');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('boxsel1').set('ymin', '1e-7');
model.component('comp').geom('ucell').feature('ucellWP').geom.run('boxsel1');
model.component('comp').geom('ucell').feature('ucellWP').geom.run('boxsel1');
model.component('comp').geom('ucell').feature('ucellWP').geom.create('fil1', 'Fillet');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('fil1').selection('point').named('boxsel1');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('fil1').set('radius', '1e-8');
model.component('comp').geom('ucell').feature('ucellWP').geom.run('fil1');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('fil1').set('radius', '3e-8');
model.component('comp').geom('ucell').feature('ucellWP').geom.run('fil1');
model.component('comp').geom('ucell').feature('ucellWP').geom.runPre('fil1');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('fil1').selection('point').set('ucellComp', [4 5 7 8]);
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('fil1').selection('point').named('boxsel1');
model.component('comp').geom('ucell').feature('ucellWP').geom.run('fil1');
model.component('comp').geom('ucell').feature('ucellWP').geom.runPre('fil1');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('fil1').selection('point').set('ucellComp', [4 7]);
model.component('comp').geom('ucell').feature('ucellWP').geom.run('fil1');
model.component('comp').geom('ucell').feature('ucellWP').geom.run('boxsel1');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('ucellComp').set('intbnd', false);
model.component('comp').geom('ucell').feature('ucellWP').geom.run('ucellComp');
model.component('comp').geom('ucell').feature('ucellWP').geom.runPre('fil1');
model.component('comp').geom('ucell').feature('ucellWP').geom.feature('fil1').selection('point').named('boxsel1');
model.component('comp').geom('ucell').feature('ucellWP').geom.run('fil1');

out = model;
