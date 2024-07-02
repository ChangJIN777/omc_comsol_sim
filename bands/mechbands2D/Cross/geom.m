function out = model
%
% geom.m
%
% Model exported on Jun 29 2024, 14:08 by COMSOL 6.2.0.415.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\mhaas\Documents\GitHub\device_design\bands\mechbands2D\Cross');

model.label('test_geom.mph');

model.param.set('k', '0');
model.param.set('a', '4e-07[m]');
model.param.set('kx', 'if(k<1,pi/a*k,if(k<2,pi/a,(3-k)*pi/a))');
model.param.set('ky', 'if(k<1,0,if(k<2,(k-1)*pi/a,(3-k)*pi/a))');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);
model.component('comp1').geom('geom1').create('wp1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp1').set('unite', true);
model.component('comp1').geom('geom1').feature('wp1').geom.create('sq1', 'Square');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('sq1').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('sq1').set('size', 4.0E-7);
model.component('comp1').geom('geom1').feature('wp1').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('size', [3.0E-7 8.0E-8]);
model.component('comp1').geom('geom1').feature('wp1').geom.create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').set('size', [8.0E-8 3.0E-7]);
model.component('comp1').geom('geom1').feature('wp1').geom.create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('dif1').selection('input').set({'sq1'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('dif1').selection('input2').set({'r1' 'r2'});
model.component('comp1').geom('geom1').feature('wp1').geom.create('fil1', 'Fillet');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('fil1').set('radius', 2.0E-8);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('fil1').selection('point').set('dif1(1)', [3 4 5 8 9 12 13 14]);
model.component('comp1').geom('geom1').feature('wp1').geom.create('fil2', 'Fillet');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('fil2').set('radius', 2.0E-8);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('fil2').selection('point').set('fil1(1)', [8 9 16 17]);
model.component('comp1').geom('geom1').create('ext1', 'Extrude');
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', '8.0E-8', 0);
model.component('comp1').geom('geom1').feature('ext1').selection('input').set({'wp1'});
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').material.create('diamond', 'Common');

model.component('comp1').physics.create('smech', 'SolidMechanics', 'geom1');
model.component('comp1').physics('smech').create('pbcX', 'PeriodicCondition', 2);
model.component('comp1').physics('smech').feature('pbcX').set('manualDestinationSelection', true);
model.component('comp1').physics('smech').feature('pbcX').selection.set([1 30]);
model.component('comp1').physics('smech').feature('pbcX').selection('destinationDomains').set([30]);
model.component('comp1').physics('smech').create('pbcY', 'PeriodicCondition', 2);
model.component('comp1').physics('smech').feature('pbcY').set('manualDestinationSelection', true);
model.component('comp1').physics('smech').feature('pbcY').selection.set([2 5]);
model.component('comp1').physics('smech').feature('pbcY').selection('destinationDomains').set([5]);
model.component('comp1').physics('smech').create('symBCs', 'SymmetrySolid', 2);
model.component('comp1').physics('smech').feature('symBCs').selection.set([3]);
model.component('comp1').physics('smech').create('asymBCs', 'Antisymmetry', 2);

model.component('comp1').material('diamond').label('diamond');
model.component('comp1').material('diamond').propertyGroup('def').set('density', '3500');
model.component('comp1').material('diamond').propertyGroup('def').set('youngsmodulus', '1.05E12');
model.component('comp1').material('diamond').propertyGroup('def').set('poissonsratio', '0.2');

model.component('comp1').physics('smech').feature('pbcX').set('PeriodicType', 'Floquet');
model.component('comp1').physics('smech').feature('pbcX').set('kFloquet', {'kx'; 'ky'; '0'});
model.component('comp1').physics('smech').feature('pbcX').label('Periodic BC, x-direction');
model.component('comp1').physics('smech').feature('pbcY').set('PeriodicType', 'Floquet');
model.component('comp1').physics('smech').feature('pbcY').set('kFloquet', {'kx'; 'ky'; '0'});
model.component('comp1').physics('smech').feature('pbcY').label('Periodic BC, y-direction');
model.component('comp1').physics('smech').feature('symBCs').label('Symmetric BC');
model.component('comp1').physics('smech').feature('asymBCs').label('Anti-symmetric BC');

out = model;
