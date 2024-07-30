function out = model
%
% ribBuildUnitCellCode.m
%
% Model exported on Jul 29 2024, 18:32 by COMSOL 6.1.0.357.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\LoncarGroup\Downloads');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('wp1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp1').set('unite', true);
model.component('comp1').geom('geom1').feature('wp1').set('quickz', '-th');

model.param.set('th', '180[nm]');
model.param.set('a', '500[nm]');
model.param.set('wo', '400[nm]');
model.param.set('ho', '200[nm]');
model.param.set('wi', '300[nm]');
model.param.set('hi', '100[nm]');
model.param.set('r1', '40[nm]');
model.param.set('r2', '45[nm]');

model.component('comp1').geom('geom1').feature('wp1').geom.create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').setIndex('table', 0, 0, 0);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').setIndex('table', 0, 0, 1);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').setIndex('table', 'a', 1, 0);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').setIndex('table', 0, 1, 1);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').setIndex('table', 'a*(3/2)', 2, 0);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').setIndex('table', 'a*sqrt(3)/2', 2, 1);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').setIndex('table', 'a*sqrt(3)/2', 3, 1);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').setIndex('table', 'a/2', 3, 0);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').setIndex('table', 0, 4, 0);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').setIndex('table', 0, 4, 1);
model.component('comp1').geom('geom1').feature('wp1').geom.run('pol1');
model.component('comp1').geom('geom1').feature('wp1').geom.feature.remove('pol1');
model.component('comp1').geom('geom1').feature('wp1').geom.create('sq1', 'Square');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('sq1').set('size', 'a');
model.component('comp1').geom('geom1').feature('wp1').geom.run('sq1');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('sq1').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.run('sq1');
model.component('comp1').geom('geom1').feature('wp1').geom.run('sq1');
model.component('comp1').geom('geom1').feature('wp1').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('size', {'wo' 'd'});

model.param.set('d', '80[nm]');
model.param.set('ai', '200[nm]');

model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('base', 'corner');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('pos', {'-wo/2' '0'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('pos', {'0' 'd/2+ai+hi'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('r1');
model.component('comp1').geom('geom1').feature('wp1').geom.run('r1');

model.param.set('hi', '50[nm]');

model.component('comp1').geom('geom1').run('wp1');
model.component('comp1').geom('geom1').runPre('fin');

model.param.set('ai', '100[nm]');

model.component('comp1').geom('geom1').feature('wp1').geom.run('r1');

model.param.set('ai', '50[nm]');

model.component('comp1').geom('geom1').feature('wp1').geom.run('r1');
model.component('comp1').geom('geom1').feature('wp1').geom.feature.duplicate('r2', 'r1');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').set('pos', {'0' '-(d/2+ai+hi)'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('r2');
model.component('comp1').geom('geom1').feature('wp1').geom.feature.duplicate('r3', 'r2');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r3').set('size', {'d2' 'd'});

model.param.set('d2', '80[nm]');

model.component('comp1').geom('geom1').feature('wp1').geom.feature('r3').set('size', {'d2' 'ho'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r3').set('pos', {'-(wo/2-d2)' '-(ho/2+ai)'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('r3');

model.param.set('ho', '100[nm]');

model.component('comp1').geom('geom1').feature('wp1').geom.run('r3');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r3').set('pos', {'-wi/2-d2/2' '-(ho/2+ai)'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('r3');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r3').set('pos', {'-wi/2+d2/2' '-(ho/2+ai)'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('r3');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r3').set('pos', {'-wo/2+d2/2' '-(ho/2+ai)'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('r3');
model.component('comp1').geom('geom1').feature('wp1').geom.feature.duplicate('r4', 'r3');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r4').set('pos', {'-wo/2+d2/2' '(ho/2+ai)'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('r4');
model.component('comp1').geom('geom1').feature('wp1').geom.feature.duplicate('r5', 'r4');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r5').set('pos', {'wo/2-d2/2' '(ho/2+ai)'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('r5');
model.component('comp1').geom('geom1').feature('wp1').geom.feature.duplicate('r6', 'r5');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r6').set('pos', {'wo/2-d2/2' '-(ho/2+ai)'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('r6');
model.component('comp1').geom('geom1').run('wp1');
model.component('comp1').geom('geom1').create('co1', 'Compose');
model.component('comp1').geom('geom1').feature.remove('co1');
model.component('comp1').geom('geom1').run('wp1');
model.component('comp1').geom('geom1').create('co1', 'Compose');
model.component('comp1').geom('geom1').feature('co1').selection('input').set({'wp1'});
model.component('comp1').geom('geom1').feature.remove('co1');
model.component('comp1').geom('geom1').feature('wp1').geom.run('r6');
model.component('comp1').geom('geom1').feature('wp1').geom.create('co1', 'Compose');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('co1').selection('input').set({'r1' 'r2' 'r3' 'r4' 'r5' 'r6' 'sq1'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('co1').set('formula', 'sq1-r1-r2-r3-r4-r5-r6');
model.component('comp1').geom('geom1').feature('wp1').geom.run('co1');
model.component('comp1').geom('geom1').feature('wp1').geom.run('co1');
model.component('comp1').geom('geom1').feature('wp1').geom.create('fil1', 'Fillet');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('fil1').set('radius', 'r1');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('fil1').selection('point').set('co1', [6 7 12 13 16 17 22 23]);
model.component('comp1').geom('geom1').feature('wp1').geom.run('fil1');

model.param.set('r1', '30[nm]');

model.component('comp1').geom('geom1').feature('wp1').geom.run('fil1');
model.component('comp1').geom('geom1').feature('wp1').geom.run('fil1');
model.component('comp1').geom('geom1').feature('wp1').geom.create('fil2', 'Fillet');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('fil2').set('radius', 'r2');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('fil2').selection('point').set('fil1', [3 10 15 18 19 22 27 34]);
model.component('comp1').geom('geom1').feature('wp1').geom.run('fil2');

model.param.set('r2', '40[nm]');

model.component('comp1').geom('geom1').feature('wp1').geom.feature.copyTo('fil3', 'geom1/wp1/fil2', 'fil2');
model.component('comp1').geom('geom1').feature('wp1').geom.feature.remove('fil2');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('fil3').tag('fil2');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('fil2').label('Fillet 2');
model.component('comp1').geom('geom1').feature('wp1').geom.run('fil1');

model.param.set('r2', '30[nm]');

model.component('comp1').geom('geom1').feature('wp1').geom.run('fil1');

model.param.set('d', '50[nm]');

model.component('comp1').geom('geom1').feature('wp1').geom.run('fil1');

model.param.set('d', '80[nm]');

model.component('comp1').geom('geom1').feature('wp1').geom.run('fil1');

model.param.set('hi', '100[nm]');

model.component('comp1').geom('geom1').feature('wp1').geom.run('fil1');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('fil2').selection('point').set('fil1', [3 8 13 16 17 20 25 30]);
model.component('comp1').geom('geom1').feature('wp1').geom.run('fil2');
model.component('comp1').geom('geom1').run('wp1');
model.component('comp1').geom('geom1').feature.create('ext1', 'Extrude');
model.component('comp1').geom('geom1').feature('ext1').set('workplane', 'wp1');
model.component('comp1').geom('geom1').feature('ext1').selection('input').set({'wp1'});
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', 'th', 0);
model.component('comp1').geom('geom1').feature('wp1').set('quickz', '-th/2');
model.component('comp1').geom('geom1').run('ext1');

out = model;
