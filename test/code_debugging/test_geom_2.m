function out = model
%
% test_geom_2.m
%
% Model exported on Jul 13 2024, 15:44 by COMSOL 6.1.0.357.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\LoncarGroup\Documents\comsol_sim\omc_comsol_sim\test\code_debugging');

model.label('test_geom.mph');

model.param.set('kx', '0');

model.component.create('comp1', false);

model.component('comp1').geom.create('geom1', 3);

model.component('comp1').label('Unit cell FEM simulation');

model.component('comp1').mesh.create('mesh');

model.component('comp1').geom('geom1').label('Unit Cell');
model.component('comp1').geom('geom1').create('wp1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp1').set('quickz', -1.15E-7);
model.component('comp1').geom('geom1').feature('wp1').geom.create('ucellplane', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('ucellplane').set('pos', [-2.0E-6 0]);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('ucellplane').set('size', [4.0E-6 2.0E-6]);
model.component('comp1').geom('geom1').feature('wp1').geom.create('hole_1_upper', 'Polygon');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('hole_1_upper').set('source', 'table');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('hole_1_upper').set('table', [-5.0E-7 7.5E-8; -3.35E-7 7.5E-8; -3.35E-7 7.5E-7; 3.35E-7 7.5E-7; 3.35E-7 7.5E-8; 5.0E-7 7.5E-8; 5.0E-7 1.4999999999999998E-6; -5.0E-7 1.4999999999999998E-6; -5.0E-7 7.5E-8]);
model.component('comp1').geom('geom1').feature('wp1').geom.create('ucellComp', 'Compose');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('ucellComp').set('formula', 'ucellplane-hole_1_upper');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('ucellComp').selection('input').set({'ucellplane' 'hole_1_upper'});
model.component('comp1').geom('geom1').create('ucellholes', 'Extrude');
model.component('comp1').geom('geom1').feature('ucellholes').setIndex('distance', '2.3E-7', 0);
model.component('comp1').geom('geom1').feature('ucellholes').selection('input').set({'wp1'});
model.component('comp1').geom('geom1').create('symZWP', 'WorkPlane');
model.component('comp1').geom('geom1').feature('symZWP').set('quickz', -1.15E-7);
model.component('comp1').geom('geom1').feature('symZWP').geom.create('symZPlane', 'Rectangle');
model.component('comp1').geom('geom1').feature('symZWP').geom.feature('symZPlane').set('pos', [-2.0E-6 -2.0E-6]);
model.component('comp1').geom('geom1').feature('symZWP').geom.feature('symZPlane').set('size', [4.0E-6 4.0E-6]);
model.component('comp1').geom('geom1').create('symZPlaneExt', 'Extrude');
model.component('comp1').geom('geom1').feature('symZPlaneExt').setIndex('distance', '1.15E-7', 0);
model.component('comp1').geom('geom1').feature('symZPlaneExt').selection('input').set({'symZWP'});
model.component('comp1').geom('geom1').create('symZComp', 'Compose');
model.component('comp1').geom('geom1').feature('symZComp').set('formula', 'ucellholes - symZPlaneExt');
model.component('comp1').geom('geom1').create('ZsymSel', 'BoxSelection');
model.component('comp1').geom('geom1').feature('ZsymSel').set('entitydim', 2);
model.component('comp1').geom('geom1').feature('ZsymSel').set('xmin', -2.01E-6);
model.component('comp1').geom('geom1').feature('ZsymSel').set('xmax', 2.01E-6);
model.component('comp1').geom('geom1').feature('ZsymSel').set('ymin', -2.01E-6);
model.component('comp1').geom('geom1').feature('ZsymSel').set('ymax', 2.01E-6);
model.component('comp1').geom('geom1').feature('ZsymSel').set('zmin', -1.0E-8);
model.component('comp1').geom('geom1').feature('ZsymSel').set('zmax', 1.0E-8);
model.component('comp1').geom('geom1').feature('ZsymSel').set('condition', 'allvertices');
model.component('comp1').geom('geom1').create('XsymSel', 'BoxSelection');
model.component('comp1').geom('geom1').feature('XsymSel').set('entitydim', 2);
model.component('comp1').geom('geom1').feature('XsymSel').set('xmin', -1.0E-8);
model.component('comp1').geom('geom1').feature('XsymSel').set('xmax', 1.0E-8);
model.component('comp1').geom('geom1').feature('XsymSel').set('ymin', -1.0E-8);
model.component('comp1').geom('geom1').feature('XsymSel').set('ymax', 2.01E-6);
model.component('comp1').geom('geom1').feature('XsymSel').set('zmin', -1.0E-8);
model.component('comp1').geom('geom1').feature('XsymSel').set('zmax', 1.25E-7);
model.component('comp1').geom('geom1').feature('XsymSel').set('condition', 'allvertices');
model.component('comp1').geom('geom1').create('XendSel', 'BoxSelection');
model.component('comp1').geom('geom1').feature('XendSel').set('entitydim', 2);
model.component('comp1').geom('geom1').feature('XendSel').set('xmin', 3.99E-6);
model.component('comp1').geom('geom1').feature('XendSel').set('xmax', 4.01E-6);
model.component('comp1').geom('geom1').feature('XendSel').set('ymin', -1.0E-8);
model.component('comp1').geom('geom1').feature('XendSel').set('ymax', 2.01E-6);
model.component('comp1').geom('geom1').feature('XendSel').set('zmin', -1.0E-8);
model.component('comp1').geom('geom1').feature('XendSel').set('zmax', 1.25E-7);
model.component('comp1').geom('geom1').feature('XendSel').set('condition', 'allvertices');
model.component('comp1').geom('geom1').create('YsymSel', 'BoxSelection');
model.component('comp1').geom('geom1').feature('YsymSel').set('entitydim', 2);
model.component('comp1').geom('geom1').feature('YsymSel').set('xmin', -1.0E-8);
model.component('comp1').geom('geom1').feature('YsymSel').set('xmax', 4.01E-6);
model.component('comp1').geom('geom1').feature('YsymSel').set('ymin', -1.0E-8);
model.component('comp1').geom('geom1').feature('YsymSel').set('ymax', 1.0E-8);
model.component('comp1').geom('geom1').feature('YsymSel').set('zmin', -1.0E-8);
model.component('comp1').geom('geom1').feature('YsymSel').set('zmax', 1.25E-7);
model.component('comp1').geom('geom1').feature('YsymSel').set('condition', 'allvertices');
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').material.create('diamond', 'Common');
model.component('comp1').material('diamond').propertyGroup.create('AnisotropicVoGrp', 'Anisotropic, Voigt notation');

model.component('comp1').physics.create('smech', 'SolidMechanics', 'geom1');
model.component('comp1').physics('smech').create('pbcX', 'PeriodicCondition', 2);
model.component('comp1').physics('smech').feature('pbcX').set('manualDestinationSelection', true);
model.component('comp1').physics('smech').feature('pbcX').selection.set([1 14]);
model.component('comp1').physics('smech').feature('pbcX').selection('destinationDomains').set([14]);
model.component('comp1').physics('smech').create('symBCs', 'SymmetrySolid', 2);
model.component('comp1').physics('smech').create('asymBCs', 'Antisymmetry', 2);
model.component('comp1').physics('smech').create('fix1', 'Fixed', 2);
model.component('comp1').physics('smech').feature('fix1').selection.set([1 5 14]);

model.component('comp1').mesh('mesh').autoMeshSize(3);

model.component('comp1').material('diamond').label('diamond');
model.component('comp1').material('diamond').propertyGroup('def').set('density', '3500');
model.component('comp1').material('diamond').propertyGroup('AnisotropicVoGrp').set('DVo', {'1.1785E12' '2.25E10' '1.25E11' '0' '0' '6.271469146668096E-6' '2.25E10' '1.1785E12' '1.25E11' '0'  ...
'0' '-6.271469146668096E-6' '1.25E11' '1.25E11' '1.076E12' '0' '0' '0' '0' '0'  ...
'0' '5.78E11' '-1.5362188012202615E-5' '0' '0' '0' '0' '-1.5362188012202615E-5' '5.78E11' '0'  ...
'6.271469146668096E-6' '-6.271469146668096E-6' '0' '0' '0' '4.755E11'});

model.component('comp1').physics('smech').feature('lemm1').set('SolidModel', 'Anisotropic');
model.component('comp1').physics('smech').feature('lemm1').set('AnisotropicOption', 'AnisotropicVo');
model.component('comp1').physics('smech').feature('pbcX').set('PeriodicType', 'Floquet');
model.component('comp1').physics('smech').feature('pbcX').set('kFloquet', {'kx'; '0'; '0'});
model.component('comp1').physics('smech').feature('pbcX').label('Periodic BC, x-direction');
model.component('comp1').physics('smech').feature('symBCs').label('Symmetric BC');
model.component('comp1').physics('smech').feature('asymBCs').label('Anti-symmetric BC');

model.study.create('study');
model.study('study').create('std_param', 'Parametric');
model.study('study').create('std_eigv', 'Eigenfrequency');

model.sol.create('solv');
model.sol('solv').study('study');
model.sol('solv').attach('study');
model.sol('solv').create('solv_stdstep', 'StudyStep');
model.sol('solv').create('solv_vars', 'Variables');
model.sol('solv').create('solv_eigv', 'Eigenvalue');
model.sol.create('psolv');
model.sol('psolv').study('study');

model.batch.create('pbatch', 'Parametric');
model.batch('pbatch').create('pbatch_solseq', 'Solutionseq');
model.batch('pbatch').study('study');

model.study('study').feature('std_param').set('pname', {'kx'});
model.study('study').feature('std_param').set('plistarr', {'0, 87266.4626, 174532.9252, 261799.3878, 349065.8504, 436332.313, 523598.7756, 610865.2382, 698131.7008, 785398.1634'});
model.study('study').feature('std_param').set('punit', {''});
model.study('study').feature('std_eigv').set('neigs', 9);
model.study('study').feature('std_eigv').set('neigsactive', true);
model.study('study').feature('std_eigv').set('shift', '0');

model.sol('solv').attach('study');
model.sol('solv').feature('solv_stdstep').label('Compile Equations: Eigenfrequency');
model.sol('solv').feature('solv_vars').label('Dependent Variables 1');
model.sol('solv').feature('solv_eigv').label('Eigenvalue Solver 1');
model.sol('solv').feature('solv_eigv').set('transform', 'eigenfrequency');
model.sol('solv').feature('solv_eigv').set('neigs', 9);
model.sol('solv').feature('solv_eigv').feature('dDef').label('Direct 1');
model.sol('solv').feature('solv_eigv').feature('dDef').set('linsolver', 'spooles');
model.sol('solv').feature('solv_eigv').feature('aDef').label('Advanced 1');

model.batch('pbatch').set('control', 'std_param');
model.batch('pbatch').set('pname', {'kx'});
model.batch('pbatch').set('plistarr', {'0, 87266.4626, 174532.9252, 261799.3878, 349065.8504, 436332.313, 523598.7756, 610865.2382, 698131.7008, 785398.1634'});
model.batch('pbatch').set('punit', {''});
model.batch('pbatch').set('err', true);
model.batch('pbatch').feature('pbatch_solseq').set('seq', 'solv');
model.batch('pbatch').feature('pbatch_solseq').set('psol', 'psolv');
model.batch('pbatch').feature('pbatch_solseq').set('param', {'"k", "0"' '"k", "87266.4626"' '"k", "174532.9252"' '"k", "261799.3878"' '"k", "349065.8504"' '"k", "436332.313"' '"k", "523598.7756"' '"k", "610865.2382"' '"k", "698131.7008"' '"k", "785398.1634"'});
model.batch('pbatch').attach('study');
model.batch('pbatch').run;

model.label('test_geom.mph');

model.component('comp1').physics('smech').feature('fix1').selection.set([5]);

out = model;
