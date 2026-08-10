function out = model
%
% trusty_boomerang_script.m
%
% Model exported on Aug 9 2026, 21:54 by COMSOL 6.3.0.290.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/Users/changjin/Documents/GitHub/omc_comsol_sim/python_script_2D/comsol');

model.label('trusty_boomerang.mph');

model.param.set('k', '0');
model.param.set('a', '4.8e-07[m]');
model.param.set('kx', 'if(k<1,(-1/sqrt(3))*k*(pi/a),if(k<2,(1/sqrt(3))*pi/a*(k-2),0))');
model.param.set('ky', 'if(k<1,(pi/a)*k,if(k<2,(k+2)*pi/(3*a),(3-k)*4*pi/(3*a)))');

model.component.create('comp1', false);

model.component('comp1').geom.create('geom1', 3);

model.component('comp1').label('Unit cell FEM simulation');

model.result.table.create('tbl1', 'Table');
model.result.evaluationGroup.create('std1EvgFrq', 'EvaluationGroup');
model.result.evaluationGroup.create('std1mpf1', 'EvaluationGroup');
model.result.evaluationGroup('std1EvgFrq').create('gev1', 'EvalGlobal');
model.result.evaluationGroup('std1mpf1').create('gev1', 'EvalGlobal');

model.component('comp1').mesh.create('mesh');

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
model.component('comp1').geom('geom1').create('ZsymSel', 'BoxSelection');
model.component('comp1').geom('geom1').feature('ZsymSel').set('entitydim', 2);
model.component('comp1').geom('geom1').feature('ZsymSel').set('xmin', -2.5E-7);
model.component('comp1').geom('geom1').feature('ZsymSel').set('xmax', 2.5E-7);
model.component('comp1').geom('geom1').feature('ZsymSel').set('ymin', -2.5E-7);
model.component('comp1').geom('geom1').feature('ZsymSel').set('ymax', 2.5E-7);
model.component('comp1').geom('geom1').feature('ZsymSel').set('zmin', -1.0E-8);
model.component('comp1').geom('geom1').feature('ZsymSel').set('zmax', 1.0E-8);
model.component('comp1').geom('geom1').feature('ZsymSel').set('condition', 'allvertices');
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').material.create('diamond', 'Common');
model.component('comp1').material('diamond').propertyGroup.create('AnisotropicVoGrp', 'AnisotropicVoGrp', 'Anisotropic, Voigt notation');

model.component('comp1').common.create('mpf1', 'ParticipationFactors');

model.component('comp1').physics.create('smech', 'SolidMechanics', 'geom1');
model.component('comp1').physics('smech').create('pbcX', 'PeriodicCondition', 2);
model.component('comp1').physics('smech').feature('pbcX').set('manualDestinationSelection', true);
model.component('comp1').physics('smech').feature('pbcX').selection.set([1 22]);
model.component('comp1').physics('smech').feature('pbcX').selection('destinationDomains').set([22]);
model.component('comp1').physics('smech').create('pbcY', 'PeriodicCondition', 2);
model.component('comp1').physics('smech').feature('pbcY').set('manualDestinationSelection', true);
model.component('comp1').physics('smech').feature('pbcY').selection.set([2 9]);
model.component('comp1').physics('smech').feature('pbcY').selection('destinationDomains').set([9]);
model.component('comp1').physics('smech').create('symBCs', 'SymmetrySolid', 2);
model.component('comp1').physics('smech').feature('symBCs').selection.set([3]);
model.component('comp1').physics('smech').create('asymBCs', 'Antisymmetry', 2);
model.component('comp1').physics('smech').feature('asymBCs').selection.set([3]);

model.component('comp1').mesh('mesh').autoMeshSize(3);

model.result.table('tbl1').label('Eigenvalues');

model.component('comp1').material('diamond').label('diamond');
model.component('comp1').material('diamond').propertyGroup('def').set('density', '3500');
model.component('comp1').material('diamond').propertyGroup('AnisotropicVoGrp').set('DVo', {'1.1785E12' '2.25E10' '1.25E11' '0' '0' '0' '2.25E10' '1.1785E12' '1.25E11' '0'  ...
'0' '0' '1.25E11' '1.25E11' '1.076E12' '0' '0' '0' '0' '0'  ...
'0' '5.78E11' '0' '0' '0' '0' '0' '0' '5.78E11' '0'  ...
'0' '0' '0' '0' '0' '4.755E11'});

model.component('comp1').physics('smech').feature('lemm1').set('SolidModel', 'Anisotropic');
model.component('comp1').physics('smech').feature('lemm1').set('AnisotropicOption', 'AnisotropicVo');
model.component('comp1').physics('smech').feature('pbcX').set('PeriodicType', 'Floquet');
model.component('comp1').physics('smech').feature('pbcX').set('kFloquet', {'kx'; 'ky'; '0'});
model.component('comp1').physics('smech').feature('pbcX').label('Periodic BC, x-direction');
model.component('comp1').physics('smech').feature('pbcY').set('PeriodicType', 'Floquet');
model.component('comp1').physics('smech').feature('pbcY').set('kFloquet', {'kx'; 'ky'; '0'});
model.component('comp1').physics('smech').feature('pbcY').label('Periodic BC, y-direction');
model.component('comp1').physics('smech').feature('symBCs').label('Symmetric BC');
model.component('comp1').physics('smech').feature('asymBCs').label('Anti-symmetric BC');

model.study.create('study');
model.study('study').create('std_param', 'Parametric');
model.study('study').create('std_eigv', 'Eigenfrequency');
model.study('study').feature('std_eigv').set('useadvanceddisable', true);
model.study('study').feature('std_eigv').set('disabledphysics', {'smech/asymBCs'});
model.study.create('std1');
model.study('std1').create('std_param', 'Parametric');
model.study('std1').create('std_eigv1', 'Eigenfrequency');
model.study('std1').feature('std_eigv1').set('useadvanceddisable', true);
model.study('std1').feature('std_eigv1').set('disabledphysics', {'smech/symBCs'});

model.batch.create('pbatch', 'Parametric');
model.batch('pbatch').create('pbatch_solseq', 'Solutionseq');
model.batch('pbatch').study('study');

model.sol.create('solv');
model.sol('solv').attach('study');
model.sol('solv').create('solv_stdstep', 'StudyStep');
model.sol('solv').create('solv_vars', 'Variables');
model.sol('solv').create('solv_eigv', 'Eigenvalue');
model.sol.create('psolv');
model.sol('psolv').study('study');
model.sol.create('sol1');
model.sol('sol1').attach('std1');
model.sol.create('sol2');
model.sol('sol2').study('std1');
model.sol('sol2').label('Parametric Solutions 1');

model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').set('data', 'dset4');
model.result('pg1').create('surf1', 'Surface');
model.result('pg1').feature('surf1').create('def', 'Deform');

model.study('study').label('Study_symmetric');
model.study('study').feature('std_param').set('pname', {'k'});
model.study('study').feature('std_param').set('plistarr', {'range(0,1/9,3-1/9)'});
model.study('study').feature('std_param').set('punit', {''});
model.study('study').feature('std_eigv').set('neigs', 10);
model.study('study').feature('std_eigv').set('neigsactive', true);
model.study('study').feature('std_eigv').set('shift', '0');
model.study('std1').feature('std_param').set('pname', {'k'});
model.study('std1').feature('std_param').set('plistarr', {'range(0,1/9,3-1/9)'});
model.study('std1').feature('std_param').set('punit', {''});
model.study('std1').feature('std_eigv1').set('neigs', 10);
model.study('std1').feature('std_eigv1').set('neigsactive', true);
model.study('std1').feature('std_eigv1').set('shift', '0');

model.batch('pbatch').set('control', 'std_param');
model.batch('pbatch').set('pname', {'k'});
model.batch('pbatch').set('plistarr', {'range(0,1/9,3-1/9)'});
model.batch('pbatch').set('punit', {''});
model.batch('pbatch').set('err', true);
model.batch('pbatch').feature('pbatch_solseq').set('seq', 'solv');
model.batch('pbatch').feature('pbatch_solseq').set('psol', 'psolv');
model.batch('pbatch').feature('pbatch_solseq').set('param', {'"k", "0"' '"k", "0.11111"' '"k", "0.22222"' '"k", "0.33333"' '"k", "0.44444"' '"k", "0.55556"' '"k", "0.66667"' '"k", "0.77778"' '"k", "0.88889"' '"k", "1"'  ...
'"k", "1.1111"' '"k", "1.2222"' '"k", "1.3333"' '"k", "1.4444"' '"k", "1.5556"' '"k", "1.6667"' '"k", "1.7778"' '"k", "1.8889"' '"k", "2"' '"k", "2.1111"'  ...
'"k", "2.2222"' '"k", "2.3333"' '"k", "2.4444"' '"k", "2.5556"' '"k", "2.6667"' '"k", "2.7778"' '"k", "2.8889"'});
model.batch('pbatch').attach('study');

model.study('std1').createAutoSequences('jobs');

model.batch('p1').feature('so1').set('psol', 'sol2');

model.sol('solv').feature('solv_stdstep').label('Compile Equations: Eigenfrequency');
model.sol('solv').feature('solv_vars').label('Dependent Variables 1');
model.sol('solv').feature('solv_eigv').label('Eigenvalue Solver 1');
model.sol('solv').feature('solv_eigv').set('neigs', 10);
model.sol('solv').feature('solv_eigv').set('shift', '0');
model.sol('solv').feature('solv_eigv').set('filtereigexpression', {'real(freq)+1e-6>0'});
model.sol('solv').feature('solv_eigv').set('filtereigdescription', {'Damped natural frequency'});
model.sol('solv').feature('solv_eigv').feature('dDef').label('Direct 1');
model.sol('solv').feature('solv_eigv').feature('dDef').set('linsolver', 'spooles');
model.sol('solv').feature('solv_eigv').feature('aDef').label('Advanced 1');
model.sol('sol1').createAutoSequence('std1');

model.study('std1').runNoGen;

model.result.evaluationGroup('std1EvgFrq').label('Eigenfrequencies (Study 1)');
model.result.evaluationGroup('std1EvgFrq').set('data', 'dset4');
model.result.evaluationGroup('std1EvgFrq').set('looplevelinput', {'all' 'all'});
model.result.evaluationGroup('std1EvgFrq').feature('gev1').set('expr', {'2*pi*freq' 'imag(freq)/abs(freq)' 'abs(freq)/imag(freq)/2'});
model.result.evaluationGroup('std1EvgFrq').feature('gev1').set('unit', {'rad/s' '1' '1'});
model.result.evaluationGroup('std1EvgFrq').feature('gev1').set('descr', {'Angular frequency' 'Damping ratio' 'Quality factor'});
model.result.evaluationGroup('std1EvgFrq').feature('gev1').set('const', {'solid.refpntx' '0' 'Reference point for moment computation, x-coordinate'; 'solid.refpnty' '0' 'Reference point for moment computation, y-coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z-coordinate'});
model.result.evaluationGroup('std1EvgFrq').run;
model.result.evaluationGroup('std1mpf1').label('Participation Factors (Study 1)');
model.result.evaluationGroup('std1mpf1').set('data', 'dset4');
model.result.evaluationGroup('std1mpf1').set('looplevelinput', {'all' 'all'});
model.result.evaluationGroup('std1mpf1').feature('gev1').set('expr', {'mpf1.pfLnormX' 'mpf1.pfLnormY' 'mpf1.pfLnormZ' 'mpf1.pfRnormX' 'mpf1.pfRnormY' 'mpf1.pfRnormZ' 'mpf1.mEffLX' 'mpf1.mEffLY' 'mpf1.mEffLZ' 'mpf1.mEffRX'  ...
'mpf1.mEffRY' 'mpf1.mEffRZ'});
model.result.evaluationGroup('std1mpf1').feature('gev1').set('unit', {'1' '1' '1' '1' '1' '1' 'kg' 'kg' 'kg' 'kg*m^2'  ...
'kg*m^2' 'kg*m^2'});
model.result.evaluationGroup('std1mpf1').feature('gev1').set('descr', {'Participation factor, normalized, X-translation' 'Participation factor, normalized, Y-translation' 'Participation factor, normalized, Z-translation' 'Participation factor, normalized, X-rotation' 'Participation factor, normalized, Y-rotation' 'Participation factor, normalized, Z-rotation' 'Effective modal mass, X-translation' 'Effective modal mass, Y-translation' 'Effective modal mass, Z-translation' 'Effective modal mass, X-rotation'  ...
'Effective modal mass, Y-rotation' 'Effective modal mass, Z-rotation'});
model.result.evaluationGroup('std1mpf1').feature('gev1').set('const', {'solid.refpntx' '0' 'Reference point for moment computation, x-coordinate'; 'solid.refpnty' '0' 'Reference point for moment computation, y-coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z-coordinate'});
model.result.evaluationGroup('std1mpf1').run;
model.result('pg1').label('Mode Shape (smech)');
model.result('pg1').set('looplevel', [9 27]);
model.result('pg1').set('showlegends', false);
model.result('pg1').feature('surf1').set('const', {'solid.refpntx' '0' 'Reference point for moment computation, x-coordinate'; 'solid.refpnty' '0' 'Reference point for moment computation, y-coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z-coordinate'});
model.result('pg1').feature('surf1').set('colortable', 'AuroraBorealis');
model.result('pg1').feature('surf1').set('threshold', 'manual');
model.result('pg1').feature('surf1').set('thresholdvalue', 0.2);
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg1').feature('surf1').feature('def').set('scale', 42628.30607276289);
model.result('pg1').feature('surf1').feature('def').set('scaleactive', false);

out = model;
