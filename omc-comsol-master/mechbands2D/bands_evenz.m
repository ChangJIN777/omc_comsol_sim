function out = model
%
% bands_evenz.m
%
% Model exported on Aug 4 2019, 14:59 by COMSOL 5.3.1.348.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('F:\cross2D\Square tethered geometry\sqTets_a_135nm_w_135nm_t_75nm_cx_115nm_cy_115nm_tx_10nm_ty_10nm_Rfil_30.0nm');

model.param.set('k', '0');
model.param.set('a', '1.35e-07[m]');
model.param.set('w', '1.35e-07[m]');
model.param.set('kx', 'if(k<1,pi/a*k,if(k<2,pi/a,(3-k)*pi/a))');
model.param.set('ky', 'if(k<1,0,if(k<2,(k-1)*pi/w,(3-k)*pi/w))');

model.geom.create('ucell', 3);
model.component('mod1').geom('ucell').feature.create('ucellWP', 'WorkPlane');
model.component('mod1').geom('ucell').feature('ucellWP').set('planetype', 'quick');
model.component('mod1').geom('ucell').feature('ucellWP').set('quickplane', 'xy');
model.component('mod1').geom('ucell').feature('ucellWP').set('quickz', -3.75E-8);
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature.create('CSq', 'Polygon');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('CSq').set('x', [0 5.75E-8 0 -5.75E-8]);
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('CSq').set('y', [5.75E-8 0 -5.75E-8 0]);
model.component('mod1').geom('ucell').feature('ucellWP').geom.create('Llink', 'Rectangle');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('Llink').set('type', 'solid');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('Llink').set('base', 'corner');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('Llink').set('pos', [-6.75E-8 -5.0E-9]);
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('Llink').set('size', [6.75E-8 1.0E-8]);
model.component('mod1').geom('ucell').run('ucellWP');
model.component('mod1').geom('ucell').feature('ucellWP').geom.create('Rlink', 'Rectangle');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('Rlink').set('type', 'solid');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('Rlink').set('base', 'corner');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('Rlink').set('pos', [0 -5.0E-9]);
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('Rlink').set('size', [6.75E-8 1.0E-8]);
model.component('mod1').geom('ucell').run('ucellWP');
model.component('mod1').geom('ucell').feature('ucellWP').geom.create('Ulink', 'Rectangle');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('Ulink').set('type', 'solid');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('Ulink').set('base', 'corner');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('Ulink').set('pos', [-5.0E-9 0]);
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('Ulink').set('size', [1.0E-8 6.75E-8]);
model.component('mod1').geom('ucell').run('ucellWP');
model.component('mod1').geom('ucell').feature('ucellWP').geom.create('Blink', 'Rectangle');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('Blink').set('type', 'solid');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('Blink').set('base', 'corner');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('Blink').set('pos', [-5.0E-9 -6.75E-8]);
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('Blink').set('size', [1.0E-8 6.75E-8]);
model.component('mod1').geom('ucell').run('ucellWP');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature.create('ucellComp', 'Compose');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('ucellComp').selection('input').set({'CSq' 'Llink' 'Rlink' 'Ulink' 'Blink'});
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('ucellComp').set('formula', 'CSq + Llink + Rlink + Ulink + Blink');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('ucellComp').set('intbnd', 'off');
model.component('mod1').geom('ucell').run('ucellWP');
model.component('mod1').geom('ucell').feature('ucellWP').geom.create('ucellFillet', 'Fillet');
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('ucellFillet').set('radius', 3.0E-8);
model.component('mod1').geom('ucell').feature('ucellWP').geom.feature('ucellFillet').selection('point').set('ucellComp', [3 4 6 7 10 11 13 14]);
model.component('mod1').geom('ucell').run('ucellWP');
model.component('mod1').geom('ucell').feature.create('ucellPlaneExt', 'Extrude');
model.component('mod1').geom('ucell').feature('ucellPlaneExt').set('distance', 7.5E-8);
model.component('mod1').geom('ucell').run('ucellPlaneExt');
model.component('mod1').geom('ucell').feature.create('ucsymZWP', 'WorkPlane');
model.component('mod1').geom('ucell').feature('ucsymZWP').set('planetype', 'quick');
model.component('mod1').geom('ucell').feature('ucsymZWP').set('quickplane', 'xy');
model.component('mod1').geom('ucell').feature('ucsymZWP').set('quickz', -3.75E-8);
model.component('mod1').geom('ucell').feature('ucsymZWP').geom.feature.create('ucsymZPlane', 'Rectangle');
model.component('mod1').geom('ucell').feature('ucsymZWP').geom.feature('ucsymZPlane').set('type', 'solid');
model.component('mod1').geom('ucell').feature('ucsymZWP').geom.feature('ucsymZPlane').set('base', 'corner');
model.component('mod1').geom('ucell').feature('ucsymZWP').geom.feature('ucsymZPlane').set('pos', [-6.75E-8 -6.75E-8]);
model.component('mod1').geom('ucell').feature('ucsymZWP').geom.feature('ucsymZPlane').set('size', [1.35E-7 1.35E-7]);
model.component('mod1').geom('ucell').run('ucsymZWP');
model.component('mod1').geom('ucell').feature.create('ucsymZPlaneExt', 'Extrude');
model.component('mod1').geom('ucell').feature('ucsymZPlaneExt').set('distance', 3.75E-8);
model.component('mod1').geom('ucell').feature.create('ucsymZComp', 'Compose');
model.component('mod1').geom('ucell').feature('ucsymZComp').selection('input').set({'ucsymZPlaneExt'});
model.component('mod1').geom('ucell').feature('ucsymZComp').set('formula', 'ucellPlaneExt - ucsymZPlaneExt');
model.component('mod1').geom('ucell').run;
model.component('mod1').geom('ucell').run;
model.component('mod1').geom('ucell').run;
model.component('mod1').geom('ucell').run;

model.component('mod1').material.create('diamond');
model.component('mod1').material('diamond').label('diamond');
model.component('mod1').material('diamond').propertyGroup('def').set('density', 3500);
model.component('mod1').material('diamond').propertyGroup.create('AnisotropicVoGrp', 'Anisotropic, Voigt notation');
model.component('mod1').material('diamond').propertyGroup('AnisotropicVoGrp').set('DVo', [1.1785E12 2.25E10 1.1785E12 1.25E11 1.25E11 1.076E12 0 0 0 5.78E11 0 0 0 0 5.78E11 0 0 0 0 0 4.755E11]);
model.component('mod1').material('diamond').selection.geom('ucell', 3);
model.component('mod1').material('diamond').selection.set([1]);

model.component('mod1').physics.create('smech', 'SolidMechanics', 'ucell');
model.component('mod1').physics('smech').feature('lemm1').set('SolidModel', 'Anisotropic');
model.component('mod1').physics('smech').feature('lemm1').set('AnisotropicOption', 'AnisotropicVo');
model.component('mod1').physics('smech').feature('lemm1').set('DVo_mat', 'from_mat');
model.component('mod1').physics('smech').feature('lemm1').set('rho_mat', 'from_mat');
model.component('mod1').physics('smech').selection.set([1]);
model.component('mod1').physics('smech').create('pbcX', 'PeriodicCondition', 2);
model.component('mod1').physics('smech').feature('pbcX').label('Periodic BC, x-direction');
model.component('mod1').physics('smech').feature('pbcX').set('PeriodicType', 'Floquet');
model.component('mod1').physics('smech').feature('pbcX').selection.set([1 26]);
model.component('mod1').physics('smech').feature('pbcX').create('perBC_dest', 'DestinationDomains', 2);
model.component('mod1').physics('smech').feature('pbcX').feature('perBC_dest').selection.set([26]);
model.component('mod1').physics('smech').feature('pbcX').set('kFloquet', {'kx' 'ky' '0'});
model.component('mod1').physics('smech').create('pbcY', 'PeriodicCondition', 2);
model.component('mod1').physics('smech').feature('pbcY').label('Periodic BC, y-direction');
model.component('mod1').physics('smech').feature('pbcY').set('PeriodicType', 'Floquet');
model.component('mod1').physics('smech').feature('pbcY').selection.set([13 15]);
model.component('mod1').physics('smech').feature('pbcY').create('perBC_dest', 'DestinationDomains', 2);
model.component('mod1').physics('smech').feature('pbcY').feature('perBC_dest').selection.set([15]);
model.component('mod1').physics('smech').feature('pbcY').set('kFloquet', {'kx' 'ky' '0'});
model.component('mod1').physics('smech').create('symBCs', 'SymmetrySolid', 2);
model.component('mod1').physics('smech').feature('symBCs').label('Symmetric BC');
model.component('mod1').physics('smech').feature('symBCs').selection.set([3]);
model.component('mod1').physics('smech').create('asymBCs', 'Antisymmetry', 2);
model.component('mod1').physics('smech').feature('asymBCs').label('Anti-symmetric BC');

model.study.create('study');
model.study('study').create('std_param', 'Parametric');
model.study('study').feature('std_param').set('pname', 'k');
model.study('study').feature('std_param').set('plistarr', 'range(0,1/5,3-1/5)');
model.study('study').feature('std_param').set('punit', []);
model.study('study').create('std_eigv', 'Eigenfrequency');
model.study('study').feature('std_eigv').set('neigsactive', true);
model.study('study').feature('std_eigv').set('neigs', 12);
model.study('study').feature('std_eigv').set('shiftactive', true);
model.study('study').feature('std_eigv').set('shift', '0');

model.sol.create('solv');
model.sol('solv').study('study');
model.sol('solv').attach('study');
model.sol('solv').create('solv_stdstep', 'StudyStep');
model.sol('solv').feature('solv_stdstep').set('study', 'study');
model.sol('solv').feature('solv_stdstep').set('studystep', 'std_eigv');
model.sol('solv').create('solv_vars', 'Variables');
model.sol('solv').feature('solv_vars').set('control', 'std_eigv');
model.sol('solv').create('solv_eigv', 'Eigenvalue');
model.sol('solv').feature('solv_eigv').set('transform', 'eigenfrequency');
model.sol('solv').feature('solv_eigv').set('control', 'std_eigv');
model.sol('solv').feature('solv_eigv').set('eigref', '0');
model.sol('solv').feature('solv_eigv').feature('dDef').set('linsolver', 'spooles');
model.sol('solv').feature('solv_eigv').feature('aDef').set('complexfun', 'off');
model.sol.create('psolv');
model.sol('psolv').study('study');

model.batch.create('pbatch', 'Parametric');
model.batch('pbatch').create('pbatch_solseq', 'Solutionseq');
model.batch('pbatch').study('study');
model.batch('pbatch').attach('study');
model.batch('pbatch').set('pname', 'k');
model.batch('pbatch').set('plistarr', 'range(0,1/5,3-1/5)');
model.batch('pbatch').set('punit', []);
model.batch('pbatch').set('err', true);
model.batch('pbatch').set('control', 'std_param');
model.batch('pbatch').feature('pbatch_solseq').set('psol', 'psolv');
model.batch('pbatch').feature('pbatch_solseq').set('param', {'"k", "0"' '"k", "0.2"' '"k", "0.4"' '"k", "0.6"' '"k", "0.8"' '"k", "1"' '"k", "1.2"' '"k", "1.4"' '"k", "1.6"' '"k", "1.8"'  ...
'"k", "2"' '"k", "2.2"' '"k", "2.4"' '"k", "2.6"' '"k", "2.8"'});
model.batch('pbatch').feature('pbatch_solseq').set('seq', 'solv');

model.component('mod1').mesh.create('mesh');
model.component('mod1').mesh('mesh').autoMeshSize(5);
model.component('mod1').mesh('mesh').run;

model.sol('solv').feature('solv_stdstep').xmeshInfo;
model.sol('solv').feature('solv_stdstep').clearXmesh;
model.sol('solv').runAll;

model.batch('pbatch').run;

model.result.dataset('dset1').tag('dset');
model.result.dataset('dset').set('solution', 'solv');
model.result.dataset('dset2').tag('pdset');
model.result.dataset('pdset').set('solution', 'psolv');
model.result.dataset.create('sec1', 'Sector3D');
model.result.dataset('sec1').set('trans', 'rotrefl');
model.result.dataset('sec1').set('pddir', {'1' '0' '0'});
model.result.dataset('sec1').set('reflaxis', {'0' '1' '0'});
model.result.dataset('sec1').set('data', 'pdset');
model.result.create('dispplot', 'PlotGroup3D');
model.result('dispplot').set('data', 'sec1');
model.result('dispplot').set('solrepresentation', 'solutioninfo');
model.result('dispplot').set('titletype', 'none');
model.result('dispplot').create('dispplot_vol', 'Volume');
model.result('dispplot').feature('dispplot_vol').set('rangedataactive', 'on');
model.result('dispplot').feature('dispplot_vol').set('rangecoloractive', 'on');
model.result('dispplot').feature('dispplot_vol').create('dispplot_def', 'Deform');
model.result('dispplot').feature('dispplot_vol').set('data', 'parent');
model.result('dispplot').run;

model.sol('solv').getPVals;
model.sol('solv').getPValsImag;
model.sol('sol1').getPVals;
model.sol('sol1').getPValsImag;
model.sol('sol2').getPVals;
model.sol('sol2').getPValsImag;
model.sol('sol3').getPVals;
model.sol('sol3').getPValsImag;
model.sol('sol4').getPVals;
model.sol('sol4').getPValsImag;
model.sol('sol5').getPVals;
model.sol('sol5').getPValsImag;
model.sol('sol6').getPVals;
model.sol('sol6').getPValsImag;
model.sol('sol7').getPVals;
model.sol('sol7').getPValsImag;
model.sol('sol8').getPVals;
model.sol('sol8').getPValsImag;
model.sol('sol9').getPVals;
model.sol('sol9').getPValsImag;
model.sol('sol10').getPVals;
model.sol('sol10').getPValsImag;
model.sol('sol11').getPVals;
model.sol('sol11').getPValsImag;
model.sol('sol12').getPVals;
model.sol('sol12').getPValsImag;
model.sol('sol13').getPVals;
model.sol('sol13').getPValsImag;
model.sol('sol14').getPVals;
model.sol('sol14').getPValsImag;
model.sol('sol15').getPVals;
model.sol('sol15').getPValsImag;
model.sol('sol1').getPVals;
model.sol('sol1').getPValsImag;
model.sol('sol2').getPVals;
model.sol('sol2').getPValsImag;
model.sol('sol3').getPVals;
model.sol('sol3').getPValsImag;
model.sol('sol4').getPVals;
model.sol('sol4').getPValsImag;
model.sol('sol5').getPVals;
model.sol('sol5').getPValsImag;
model.sol('sol6').getPVals;
model.sol('sol6').getPValsImag;
model.sol('sol7').getPVals;
model.sol('sol7').getPValsImag;
model.sol('sol8').getPVals;
model.sol('sol8').getPValsImag;
model.sol('sol9').getPVals;
model.sol('sol9').getPValsImag;
model.sol('sol10').getPVals;
model.sol('sol10').getPValsImag;
model.sol('sol11').getPVals;
model.sol('sol11').getPValsImag;
model.sol('sol12').getPVals;
model.sol('sol12').getPValsImag;
model.sol('sol13').getPVals;
model.sol('sol13').getPValsImag;
model.sol('sol14').getPVals;
model.sol('sol14').getPValsImag;
model.sol('sol15').getPVals;
model.sol('sol15').getPValsImag;

model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol1').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'1' '1'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 1.7320508075695933);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 1.7320508075695933);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol1').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 2);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '2');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'2' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'2' '1'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 1.732050807569773);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 1.732050807569773);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol1').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 3);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '3');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'3' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'3' '1'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 2.474806500612668);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 2.474806500612668);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol1').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 4);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '4');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'4' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'4' '1'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 4.3430738357053285);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 4.3430738357053285);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol1').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 5);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '5');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'5' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'5' '1'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 4.378266956830871);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 4.378266956830871);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol1').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 6);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '6');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'6' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'6' '1'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 7.047569114429588);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 7.047569114429588);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol1').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 7);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '7');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'7' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'7' '1'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 7.021278628252425);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 7.021278628252425);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol1').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 8);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '8');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'8' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'8' '1'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 2.886287377819946);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 2.886287377819946);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol1').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 9);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '9');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'9' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'9' '1'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 4.641727857883282);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 4.641727857883282);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol1').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 10);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '10');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'10' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'10' '1'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 4.400028176682789);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 4.400028176682789);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol1').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 11);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '11');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'11' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'11' '1'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 4.260241994716465);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 4.260241994716465);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol1').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 12);
model.result.numerical('num1').set('outersolnum', 1);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '12');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'12' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'12' '1'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 3.3460310687503654);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 3.3460310687503654);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol6').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'1' '6'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 1.9237111573001746);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 1.9237111573001746);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol6').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 2);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '2');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'2' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'2' '6'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 2.9079059139262236);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 2.9079059139262236);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol6').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 3);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '3');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'3' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'3' '6'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 2.507051370797384);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 2.507051370797384);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol6').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 4);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '4');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'4' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'4' '6'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 4.448963705631214);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 4.448963705631214);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol6').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 5);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '5');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'5' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'5' '6'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 4.519454468794145);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 4.519454468794145);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol6').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 6);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '6');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'6' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'6' '6'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 3.8131471324625656);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 3.8131471324625656);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol6').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 7);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '7');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'7' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'7' '6'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 3.3700446692759796);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 3.3700446692759796);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol6').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 8);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '8');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'8' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'8' '6'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 6.658197910021697);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 6.658197910021697);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol6').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 9);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '9');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'9' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'9' '6'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 3.0148155398664676);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 3.0148155398664676);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol6').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 10);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '10');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'10' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'10' '6'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 7.3280434290284315);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 7.3280434290284315);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol6').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 11);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '11');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'11' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'11' '6'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 4.632501566059145);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 4.632501566059145);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol6').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 12);
model.result.numerical('num1').set('outersolnum', 6);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '12');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'12' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'12' '6'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 3.8197427728280657);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 3.8197427728280657);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol11').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'1' '11'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 2.3495599035866346);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 2.3495599035866346);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol11').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 2);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '2');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'2' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'2' '11'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 2.244307510587817);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 2.244307510587817);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol11').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 3);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '3');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'3' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'3' '11'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 2.2387871549936955);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 2.2387871549936955);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol11').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 4);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '4');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'4' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'4' '11'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 3.244956778852038);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 3.244956778852038);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol11').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 5);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '5');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'5' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'5' '11'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 2.341740181527957);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 2.341740181527957);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol11').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 6);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '6');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'6' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'6' '11'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 3.8971646669159923);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 3.8971646669159923);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol11').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 7);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '7');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'7' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'7' '11'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 3.9390395831884497);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 3.9390395831884497);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol11').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 8);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '8');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'8' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'8' '11'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 2.529441371862595);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 2.529441371862595);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol11').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 9);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '9');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'9' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'9' '11'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 5.206244633017482);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 5.206244633017482);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol11').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 10);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '10');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'10' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'10' '11'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 3.7272749969174432);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 3.7272749969174432);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol11').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 11);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '11');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'11' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'11' '11'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 3.118552510881136);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 3.118552510881136);
model.result('dispplot').run;
model.result.numerical.create('num1', 'MaxVolume');
model.result.numerical('num1').set('data', 'pdset');
model.result.numerical('num1').selection.all;

model.sol('psolv').getSizeMulti;
model.sol('psolv').getSize;
model.sol('psolv').getPVals;
model.sol('psolv').getPValsImag;

model.result.numerical('num1').set('expr', 'abs(solid.disp)');
model.result.numerical('num1').run;
model.result.numerical('num1').set('solnumindices', 1);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '1');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'1' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');

model.sol('sol11').getSizeMulti;

model.result.numerical('num1').set('solnumindices', 12);
model.result.numerical('num1').set('outersolnum', 11);
model.result.numerical('num1').set('solnum', '');
model.result.numerical('num1').set('outersolnum', '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15');
model.result.numerical('num1').set('solnumindices', '12');
model.result.numerical('num1').set('outersolnumindices', '');
model.result.numerical('num1').set('innerinput', 'manualindices');
model.result.numerical('num1').set('outerinput', 'all');
model.result.numerical('num1').set('looplevelinput', {'manualindices' 'all'});
model.result.numerical('num1').set('looplevel', {'' '' '' '' '' '' '' '' '' ''  ...
'' '' '' '' '';  ...
'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'  ...
'11' '12' '13' '14' '15'});
model.result.numerical('num1').set('looplevelindices', {'12' ''});
model.result.numerical('num1').set('interp', {'' ''});
model.result.numerical('num1').set('solrepresentation', 'solutioninfo');
model.result.numerical.remove('num1');
model.result('dispplot').set('looplevel', {'12' '11'});
model.result('dispplot').feature('dispplot_vol').set('expr', 'abs(solid.disp)');
model.result('dispplot').feature('dispplot_vol').set('rangedatamin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangedatamax', 3.8960551840560034);
model.result('dispplot').feature('dispplot_vol').set('rangecolormin', '0');
model.result('dispplot').feature('dispplot_vol').set('rangecolormax', 3.8960551840560034);
model.result('dispplot').run;

model.label('bands_evenz.mph');

model.component('mod1').geom('ucell').feature('ucellWP').geom.run('CSq');
model.component('mod1').geom('ucell').feature('ucellWP').geom.run('Llink');
model.component('mod1').geom('ucell').feature('ucellWP').geom.run('Rlink');
model.component('mod1').geom('ucell').feature('ucellWP').geom.run('Ulink');
model.component('mod1').geom('ucell').feature('ucellWP').geom.run('Blink');
model.component('mod1').geom('ucell').feature('ucellWP').geom.run('ucellComp');
model.component('mod1').geom('ucell').feature('ucellWP').geom.run('ucellFillet');
model.component('mod1').geom('ucell').run('ucsymZPlaneExt');
model.component('mod1').geom('ucell').run('ucsymZComp');

out = model;
