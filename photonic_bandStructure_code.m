function out = model
%
% photonic_bandStructure_code.m
%
% Model exported on Jul 26 2024, 17:44 by COMSOL 6.1.0.357.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\LoncarGroup\Documents\comsol_sim\omc_comsol_sim');

model.param.set('kx', 'if(k<1,(pi/a)*k*(sqrt(3)/2),if(k<2,(sqrt(3)/2)*pi/a,(sqrt(3)/2)*(3-k)*pi/a))');
model.param.set('ky', 'if(k<1,(pi/a)*k*(-1/2),if(k<2,(k-1-1/2)*pi/a,(1/2)*(3-k)*pi/a))');
model.param.set('a', '3.6e-07[m]');
model.param.set('k', '0');

model.component.create('comp1', false);

model.component('comp1').geom.create('geom1', 2);

model.component('comp1').label('Unit cell FEM simulation');

model.result.table.create('tbl1', 'Table');
model.result.table.create('tbl2', 'Table');
model.result.table.create('tbl3', 'Table');
model.result.table.create('tbl4', 'Table');
model.result.table.create('tbl5', 'Table');
model.result.table.create('tbl6', 'Table');
model.result.table.create('tbl7', 'Table');
model.result.table.create('tbl8', 'Table');
model.result.table.create('tbl9', 'Table');
model.result.table.create('tbl10', 'Table');
model.result.table.create('tbl11', 'Table');
model.result.table.create('tbl12', 'Table');
model.result.table.create('tbl13', 'Table');
model.result.table.create('tbl14', 'Table');

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').label('Unit Cell');
model.component('comp1').geom('geom1').create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature('pol1').label('Base plane');
model.component('comp1').geom('geom1').feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature('pol1').set('table', {'0' '0'; 'a/2' 'sqrt(3)*a/2'; 'a*(3/2)' 'sqrt(3)*a/2'; 'a' '0'; '0' '0'});
model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('pos', [2.1749999999999998E-7 2.063990010960491E-7]);
model.component('comp1').geom('geom1').feature('r1').set('base', 'center');
model.component('comp1').geom('geom1').feature('r1').set('size', [7.5E-8 1.4E-7]);
model.component('comp1').geom('geom1').create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('r2').set('pos', [2.1749999999999998E-7 7.144709581221619E-8]);
model.component('comp1').geom('geom1').feature('r2').set('rot', 120);
model.component('comp1').geom('geom1').feature('r2').set('size', [7.5E-8 1.4E-7]);
model.component('comp1').geom('geom1').create('r3', 'Rectangle');
model.component('comp1').geom('geom1').feature('r3').set('pos', [2.55E-7 1.363990010960491E-7]);
model.component('comp1').geom('geom1').feature('r3').set('rot', 240);
model.component('comp1').geom('geom1').feature('r3').set('size', [7.5E-8 1.4E-7]);
model.component('comp1').geom('geom1').create('fil1', 'Fillet');
model.component('comp1').geom('geom1').feature('fil1').set('radius', 3.5E-8);
model.component('comp1').geom('geom1').feature('fil1').selection('point').set('r1(1)', [3 4]);
model.component('comp1').geom('geom1').feature('fil1').selection('point').set('r3(1)', [3 4]);
model.component('comp1').geom('geom1').feature('fil1').selection('point').set('r2(1)', [3 4]);
model.component('comp1').geom('geom1').create('pol2', 'Polygon');
model.component('comp1').geom('geom1').feature('pol2').set('source', 'table');
model.component('comp1').geom('geom1').feature('pol2').set('table', [1.8E-7 1.363990010960491E-7; 2.55E-7 1.363990010960491E-7; 2.1749999999999998E-7 7.144709581221619E-8; 1.8E-7 1.363990010960491E-7]);
model.component('comp1').geom('geom1').create('co1', 'Compose');
model.component('comp1').geom('geom1').feature('co1').set('formula', 'pol1-(fil1(1)+fil1(2)+fil1(3)+pol2)');
model.component('comp1').geom('geom1').create('fil2', 'Fillet');
model.component('comp1').geom('geom1').feature('fil2').set('radius', 4.0E-8);
model.component('comp1').geom('geom1').feature('fil2').selection('point').set('co1(1)', [6 10 12]);
model.component('comp1').geom('geom1').create('pol3', 'Polygon');
model.component('comp1').geom('geom1').feature('pol3').label('Base plane 1');
model.component('comp1').geom('geom1').feature('pol3').set('source', 'table');
model.component('comp1').geom('geom1').feature('pol3').set('table', {'0' '0'; 'a/2' 'sqrt(3)*a/2'; 'a*(3/2)' 'sqrt(3)*a/2'; 'a' '0'; '0' '0'});
model.component('comp1').geom('geom1').run;

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material('mat1').selection.set([1]);
model.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat2').selection.set([2]);
model.component('comp1').material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');

model.component('comp1').physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom1');
model.component('comp1').physics('ewfd').create('pc1', 'PeriodicCondition', 1);
model.component('comp1').physics('ewfd').feature('pc1').selection.set([1 13]);
model.component('comp1').physics('ewfd').create('pc2', 'PeriodicCondition', 1);
model.component('comp1').physics('ewfd').feature('pc2').selection.set([2 7]);

model.component('comp1').view('view1').axis.set('xmin', -1.3499963813501381E-8);
model.component('comp1').view('view1').axis.set('xmax', 5.534999445444555E-7);
model.component('comp1').view('view1').axis.set('ymin', -6.898315518810705E-8);
model.component('comp1').view('view1').axis.set('ymax', 3.8075231145739963E-7);

model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'2.4' '0' '0' '0' '2.4' '0' '0' '0' '2.4'});
model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});

model.component('comp1').physics('ewfd').prop('components').set('components', 'inplane');
model.component('comp1').physics('ewfd').feature('pc1').set('PeriodicType', 'Floquet');
model.component('comp1').physics('ewfd').feature('pc1').set('kFloquet', {'kx'; 'ky'; '0'});
model.component('comp1').physics('ewfd').feature('pc1').label('Periodic Condition x direction');
model.component('comp1').physics('ewfd').feature('pc2').set('PeriodicType', 'Floquet');
model.component('comp1').physics('ewfd').feature('pc2').set('kFloquet', {'kx'; 'ky'; '0'});
model.component('comp1').physics('ewfd').feature('pc2').label('Periodic Condition y direction');

model.study.create('std1');
model.study('std1').create('param', 'Parametric');
model.study('std1').create('eig', 'Eigenfrequency');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('e1', 'Eigenvalue');
model.sol('sol1').feature('e1').create('d1', 'Direct');
model.sol.create('sol2');
model.sol('sol2').study('std1');
model.sol('sol2').label('Parametric Solutions 1');
model.sol.create('sol3');
model.sol('sol3').study('std1');
model.sol('sol3').label('Parametric Solutions 2');

model.batch.create('p1', 'Parametric');
model.batch('p1').create('so1', 'Solutionseq');
model.batch('p1').study('std1');

model.result.numerical.create('gev1', 'EvalGlobal');
model.result.numerical.create('gev2', 'EvalGlobal');
model.result.numerical('gev1').set('probetag', 'none');
model.result.numerical('gev2').set('data', 'dset3');
model.result.numerical('gev2').set('probetag', 'none');
model.result.create('pg1', 'PlotGroup2D');
model.result.create('pg2', 'PlotGroup2D');
model.result.create('pg3', 'PlotGroup1D');
model.result('pg1').create('surf1', 'Surface');
model.result('pg2').set('data', 'dset3');
model.result('pg2').create('surf1', 'Surface');
model.result('pg3').set('data', 'dset3');
model.result('pg3').create('glob1', 'Global');

model.study('std1').feature('param').set('pname', {'k'});
model.study('std1').feature('param').set('plistarr', {'range(0,1/9,3-1/9)'});
model.study('std1').feature('param').set('punit', {''});
model.study('std1').feature('eig').set('neigsactive', true);
model.study('std1').feature('eig').set('shift', '226');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Eigenfrequency');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('e1').label('Eigenvalue Solver 1.1');
model.sol('sol1').feature('e1').set('transform', 'eigenfrequency');
model.sol('sol1').feature('e1').set('shift', '226');
model.sol('sol1').feature('e1').feature('dDef').label('Direct 2');
model.sol('sol1').feature('e1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
model.sol('sol1').feature('e1').feature('d1').label('Suggested Direct Solver (ewfd)');
model.sol('sol1').runAll;

model.batch('p1').set('control', 'param');
model.batch('p1').set('pname', {'k'});
model.batch('p1').set('plistarr', {'range(0,1/9,3-1/9)'});
model.batch('p1').set('punit', {''});
model.batch('p1').set('err', true);
model.batch('p1').feature('so1').set('seq', 'sol1');
model.batch('p1').feature('so1').set('psol', 'sol3');
model.batch('p1').feature('so1').set('param', {'"k","0"' '"k","0.111111111111111"' '"k","0.222222222222222"' '"k","0.333333333333333"' '"k","0.444444444444444"' '"k","0.555555555555556"' '"k","0.666666666666667"' '"k","0.777777777777778"' '"k","0.888888888888889"' '"k","1"'  ...
'"k","1.11111111111111"' '"k","1.22222222222222"' '"k","1.33333333333333"' '"k","1.44444444444444"' '"k","1.55555555555556"' '"k","1.66666666666667"' '"k","1.77777777777778"' '"k","1.88888888888889"' '"k","2"' '"k","2.11111111111111"'  ...
'"k","2.22222222222222"' '"k","2.33333333333333"' '"k","2.44444444444444"' '"k","2.55555555555556"' '"k","2.66666666666667"' '"k","2.77777777777778"' '"k","2.88888888888889"'});
model.batch('p1').attach('std1');
model.batch('p1').run;

model.result.numerical('gev1').label('Eigenfrequencies (ewfd)');
model.result.numerical('gev1').set('table', 'tbl13');
model.result.numerical('gev1').set('expr', {'ewfd.freq' 'ewfd.Qfactor'});
model.result.numerical('gev1').set('unit', {'THz' '1'});
model.result.numerical('gev1').set('descr', {'Frequency' 'Quality factor'});
model.result.numerical('gev2').label('Eigenfrequencies (ewfd) 1');
model.result.numerical('gev2').set('table', 'tbl14');
model.result.numerical('gev2').set('expr', {'ewfd.freq' 'ewfd.Qfactor'});
model.result.numerical('gev2').set('unit', {'THz' '1'});
model.result.numerical('gev2').set('descr', {'Frequency' 'Quality factor'});
model.result.numerical('gev1').setResult;
model.result.numerical('gev2').setResult;
model.result('pg1').label('Electric Field (ewfd)');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature('surf1').set('smooth', 'internal');
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg2').label('Electric Field (ewfd) 1');
model.result('pg2').set('frametype', 'spatial');
model.result('pg2').feature('surf1').set('smooth', 'internal');
model.result('pg2').feature('surf1').set('resolution', 'normal');
model.result('pg3').set('xlabel', 'k');
model.result('pg3').set('ylabel', 'Frequency (THz)');
model.result('pg3').set('xlabelactive', false);
model.result('pg3').set('ylabelactive', false);
model.result('pg3').feature('glob1').set('expr', {'ewfd.freq'});
model.result('pg3').feature('glob1').set('unit', {'THz'});
model.result('pg3').feature('glob1').set('descr', {'Frequency'});
model.result('pg3').feature('glob1').set('xdatasolnumtype', 'outer');
model.result('pg3').feature('glob1').set('linewidth', 'preference');

out = model;
