function out = model
%
% Untitled.m
%
% Model exported on Jul 7 2024, 15:26 by COMSOL 6.1.0.357.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\LoncarGroup\Documents\comsol_sim\omc_comsol_sim\omc-comsol-chang\test\rib\070724');

model.param.set('k', '0');
model.param.set('a', '2e-07[m]');
model.param.set('kx', 'if(k<1,pi/a*k,if(k<2,pi/a,(3-k)*pi/a))');
model.param.set('ky', 'if(k<1,0,if(k<2,(k-1)*pi/a,(3-k)*pi/a))');

model.component.create('comp1');

model.component('comp1').label('Unit cell FEM simulation');

model.component('comp1').geom.create('geom1', 3);
model.component('comp1').geom('geom1').label('Unit Cell');
model.component('comp1').geom('geom1').feature.create('wp1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp1').set('planetype', 'quick');
model.component('comp1').geom('geom1').feature('wp1').set('quickplane', 'xy');
model.component('comp1').geom('geom1').feature('wp1').set('quickz', -8.0E-8);
model.component('comp1').geom('geom1').feature('wp1').geom.feature.create('ucellplane', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('ucellplane').set('type', 'solid');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('ucellplane').set('base', 'corner');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('ucellplane').set('pos', [-1.0E-7 -1.0E-7]);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('ucellplane').set('size', [2.0E-7 2.0E-7]);
model.component('comp1').geom('geom1').run('wp1');
model.component('comp1').geom('geom1').feature('wp1').geom.feature.create('hole_1_upper', 'Polygon');
model.component('comp1').geom('geom1').feature('wp1').geom.feature.create('hole_1_lower', 'Polygon');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('hole_1_upper').set('source', 'table');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('hole_1_lower').set('source', 'table');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('hole_1_upper').set('table', [-6.0E-8 5.0E-8; -3.0E-8 5.0E-8; -3.0E-8 6.999999999999999E-8; 3.0E-8 6.999999999999999E-8; 3.0E-8 5.0E-8; 6.0E-8 5.0E-8; 6.0E-8 8.0E-8; -6.0E-8 8.0E-8; -6.0E-8 5.0E-8]);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('hole_1_lower').set('table', [-6.0E-8 -5.0E-8; -3.0E-8 -5.0E-8; -3.0E-8 -6.999999999999999E-8; 3.0E-8 -6.999999999999999E-8; 3.0E-8 -5.0E-8; 6.0E-8 -5.0E-8; 6.0E-8 -8.0E-8; -6.0E-8 -8.0E-8; -6.0E-8 -5.0E-8]);
model.component('comp1').geom('geom1').run('wp1');
model.component('comp1').geom('geom1').feature('wp1').geom.feature.create('ucellComp', 'Compose');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('ucellComp').selection('input').set({'ucellplane'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('ucellComp').set('formula', 'ucellplane-hole_1_upper-hole_1_lower');
model.component('comp1').geom('geom1').run('wp1');
model.component('comp1').geom('geom1').feature.create('ucellholes', 'Extrude');
model.component('comp1').geom('geom1').feature('ucellholes').set('distance', 1.6E-7);
model.component('comp1').geom('geom1').run('ucellholes');
model.component('comp1').geom('geom1').feature.create('symZWP', 'WorkPlane');
model.component('comp1').geom('geom1').feature('symZWP').set('planetype', 'quick');
model.component('comp1').geom('geom1').feature('symZWP').set('quickplane', 'xy');
model.component('comp1').geom('geom1').feature('symZWP').set('quickz', -8.0E-8);
model.component('comp1').geom('geom1').feature('symZWP').geom.feature.create('symZPlane', 'Rectangle');
model.component('comp1').geom('geom1').feature('symZWP').geom.feature('symZPlane').set('type', 'solid');
model.component('comp1').geom('geom1').feature('symZWP').geom.feature('symZPlane').set('base', 'corner');
model.component('comp1').geom('geom1').feature('symZWP').geom.feature('symZPlane').set('pos', [-1.0E-7 -1.0E-7]);
model.component('comp1').geom('geom1').feature('symZWP').geom.feature('symZPlane').set('size', [2.0E-7 2.0E-7]);
model.component('comp1').geom('geom1').run('symZWP');
model.component('comp1').geom('geom1').feature.create('symZPlaneExt', 'Extrude');
model.component('comp1').geom('geom1').feature('symZPlaneExt').set('distance', 8.0E-8);
model.component('comp1').geom('geom1').feature.create('symZComp', 'Compose');
model.component('comp1').geom('geom1').feature('symZComp').selection('input').set({'symZPlaneExt'});
model.component('comp1').geom('geom1').feature('symZComp').set('formula', 'ucellholes - symZPlaneExt');
model.component('comp1').geom('geom1').run('symZComp');
model.component('comp1').geom('geom1').create('ZsymSel', 'BoxSelection');
model.component('comp1').geom('geom1').feature('ZsymSel').set('xmin', -1.0999999999999999E-7);
model.component('comp1').geom('geom1').feature('ZsymSel').set('xmax', 1.0999999999999999E-7);
model.component('comp1').geom('geom1').feature('ZsymSel').set('ymin', -1.0999999999999999E-7);
model.component('comp1').geom('geom1').feature('ZsymSel').set('ymax', 1.0999999999999999E-7);
model.component('comp1').geom('geom1').feature('ZsymSel').set('zmin', -1.0E-8);
model.component('comp1').geom('geom1').feature('ZsymSel').set('zmax', 1.0E-8);
model.component('comp1').geom('geom1').feature('ZsymSel').set('entitydim', 2);
model.component('comp1').geom('geom1').feature('ZsymSel').set('condition', 'allvertices');
model.component('comp1').geom('geom1').run('ZsymSel');

model.component('comp1').material.create('diamond');
model.component('comp1').material('diamond').label('diamond');
model.component('comp1').material('diamond').propertyGroup('def').set('density', 3500);
model.component('comp1').material('diamond').propertyGroup.create('AnisotropicVoGrp', 'Anisotropic, Voigt notation');
model.component('comp1').material('diamond').propertyGroup('AnisotropicVoGrp').set('DVo', [1.1785E12 2.25E10 1.1785E12 1.25E11 1.25E11 1.076E12 0 0 0 5.78E11 0 0 0 -1.5362188012202615E-5 5.78E11 6.271469146668096E-6 -6.271469146668096E-6 0 0 0 4.755E11]);
model.component('comp1').material('diamond').selection.geom('geom1', 3);

model.component('comp1').geom('geom1').run;

model.component('comp1').material('diamond').selection.set([1]);

out = model;
