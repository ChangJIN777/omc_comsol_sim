function out = model
%
% trusty_boomerang_script.m
%
% Model exported on Aug 9 2026, 21:54 by COMSOL 6.3.0.290.
% Modified: geometry parameterized on a, w, r, r1, r2, th (was hard-coded
% literals -- only k, a, kx, ky were parameters). Defaults match
% omc-comsol-chang/test_Boomerang.m:9-16.
%
% HOW TO USE. This .m is an EXPORT: editing it does NOT change
% trusty_boomerang.mph. To get a parametric .mph, run this in MATLAB with
% LiveLink for COMSOL and save the result:
%
%     model = trusty_boomerang_script();
%     mphsave(model, 'trusty_boomerang.mph');
%
% Note the script solves as written (model.study('std1').runNoGen near the
% end, plus the evaluationGroup .run calls). Comment those out to build
% geometry + physics only, which is all you need to check the geometry.
%
% STILL NOT DONE -- read before varying any parameter:
%
%  1. The periodic and parity BCs are still pinned to ABSOLUTE face indices
%     (pbcX [1 22], pbcY [2 9], symBCs/asymBCs [3]). COMSOL renumbers faces
%     whenever geometry topology changes, and now that the geometry is
%     parametric it WILL change -- growing r past the cell wall splits a side
%     face, changing r1/r2 adds or removes fillet faces. The BCs then attach
%     to the WRONG faces and COMSOL solves a meaningless problem without
%     complaining. Replace them with bndindex(...)-style geometric selections
%     (see buildBoomerangUnitCell.m:106-112) before trusting any swept result.
%     Parameterizing the geometry ACTIVATES this latent bug; it was dormant
%     only because nothing ever changed.
%
%  2. Study labels are still 'Study_symmetric' and (unset ->) 'Study 1'.
%     src/acoustic_comsol_2d.py resolves studies by label and expects
%     'mech evenz' / 'mech oddz'.
%
%  3. geomRep('cadps') selects the Parasolid kernel, which needs the CAD
%     Import Module / Design Module. A headless mphserver on a
%     Structural-Mechanics-only licence will now fail at geometry REBUILD --
%     it did not before, because the geometry was never rebuilt.
%     buildBoomerangUnitCell.m never sets geomRep.
%
% VERIFY FIRST: rebuild at the defaults and confirm the eigenfrequencies are
% unchanged from the pre-parameterization model. That invariant is the only
% check that this edit preserved the geometry.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/Users/changjin/Documents/GitHub/omc_comsol_sim/python_script_2D/comsol');

model.label('trusty_boomerang.mph');

%% Parameters
% Geometry. Values match omc-comsol-chang/test_Boomerang.m:9-16, so this model
% reproduces the MATLAB reference design at its defaults.
%
% FILLET CONVENTION -- r1 is the JUNCTION fillet and r2 is the TIP fillet.
% Authority is addFillet() at buildBoomerangUnitCell.m:119-149: h_disksel1 is
% the annulus [w/(2*sqrt(2)), w] about the hole centre, catching the vertices
% at w/2 (the inner corners where the three legs meet), and h_fil1 applies r1;
% h_disksel2 is [r-selw, r+selw], catching the vertices at hypot(r, w/2) (the
% outer tip corners), and h_fil2 applies r2. Several docs had this backwards.
model.param.set('a',  '480[nm]', 'Lattice constant (rhombic primitive cell side)');
model.param.set('w',  '140[nm]', 'Boomerang leg width, transverse to the leg');
model.param.set('r',  '177[nm]', 'Boomerang leg radial length from hole centre');
model.param.set('r1', '10[nm]',  'Fillet radius at leg JUNCTIONS (inner corners)');
model.param.set('r2', '10[nm]',  'Fillet radius at leg TIPS (outer corners)');
model.param.set('th', '220[nm]', 'Slab thickness (FULL; only the top half is meshed)');

% Derived geometry helpers -- keep the geometry-node expressions readable and
% one-to-one with buildBoomerangUnitCell.m.
model.param.set('hx', 'a*(1/2+1/4)', 'Hole centre x (rhombus centroid)');
model.param.set('hy', 'a*sqrt(3)/4', 'Hole centre y (= cell in-radius)');
model.param.set('selw', '25[nm]', 'Fillet tip-annulus half-width (selection_width/2)');
model.param.set('dsel', '10[nm]', 'Box-selection tolerance');

% Brillouin-zone sweep. kx/ky MUST stay expressions in k -- the study sweeps k
% and the Python driver should set k, not kx/ky. Writing literals into kx/ky
% destroys these expressions and every sweep point then solves the same
% wavevector. Exact port of runBands_2D.m:63-64.
model.param.set('k', '0', 'BZ path parameter, 0 to 3: Gamma-M-K-Gamma');
model.param.set('kx', 'if(k<1,(-1/sqrt(3))*k*(pi/a),if(k<2,(1/sqrt(3))*pi/a*(k-2),0))', 'Floquet wavevector, x');
model.param.set('ky', 'if(k<1,(pi/a)*k,if(k<2,(k+2)*pi/(3*a),(3-k)*4*pi/(3*a)))', 'Floquet wavevector, y');

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
% Work plane at the slab's lower face, z = -th/2.
model.component('comp1').geom('geom1').feature('wp1').set('quickz', '-th/2');

% Rhombic primitive cell of the triangular lattice: side a, 60/120 deg angles.
% Its centre-to-edge distance is a*sqrt(3)/4 (= hy), NOT a*sqrt(3)/2.
model.component('comp1').geom('geom1').feature('wp1').geom.create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').label('Base plane');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pol1').set('table', {'0' '0'; 'a/2' '(a/2)*sqrt(3)'; 'a*(3/2)' '(a/2)*sqrt(3)'; 'a' '0'; '0' '0'});

% The three boomerang legs, each w wide and r long, radiating from the hole
% centre (hx,hy) at 0/120/240 deg.
%
% TAGS RENAMED r1,r2,r3 -> leg1,leg2,leg3. The original tags collided with the
% new r1/r2 PARAMETERS. COMSOL parses geometry-object names and numeric
% expressions with different parsers, so 'radius','r1' and the Compose formula
% 'pol1-r1-...' would probably still resolve correctly -- but relying on that
% is not worth it when a rename costs nothing. The Compose formula below is
% updated to match.
model.component('comp1').geom('geom1').feature('wp1').geom.create('leg1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('leg1').set('size', {'w' 'r'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('leg1').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('leg1').set('pos', {'hx' 'hy+r/2'});
model.component('comp1').geom('geom1').feature('wp1').geom.create('leg2', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('leg2').set('size', {'w' 'r'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('leg2').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('leg2').set('pos', {'hx-sqrt(3)*r/4' 'hy-r/4'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('leg2').set('rot', 120);
model.component('comp1').geom('geom1').feature('wp1').geom.create('leg3', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('leg3').set('size', {'w' 'r'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('leg3').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('leg3').set('pos', {'hx+sqrt(3)*r/4' 'hy-r/4'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('leg3').set('rot', 240);
model.component('comp1').geom('geom1').feature('wp1').geom.create('co1', 'Compose');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('co1').set('formula', 'pol1-leg1-leg2-leg3');

% Junction fillets: annulus [w/(2*sqrt(2)), w] catches the inner corners at w/2.
model.component('comp1').geom('geom1').feature('wp1').geom.create('h_disksel1', 'DiskSelection');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel1').set('entitydim', 0);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel1').set('posx', 'hx');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel1').set('posy', 'hy');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel1').set('r', 'w');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel1').set('rin', 'w/(2*sqrt(2))');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel1').set('condition', 'allvertices');
model.component('comp1').geom('geom1').feature('wp1').geom.create('h_fil1', 'Fillet');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_fil1').set('radius', 'r1');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_fil1').selection('point').named('h_disksel1');

% Tip fillets: annulus [r-selw, r+selw] must contain the tip corners, which sit
% at hypot(r, w/2) -- note the annulus is centred on r, not on hypot(r, w/2),
% so it only works while hypot(r,w/2) <= r+selw. At the reference design the
% tips are at 190.3 nm inside [152, 202] nm. src/geometry2d.py:check_feasibility
% rejects candidates outside that window. Centring on 'sqrt(r^2+(w/2)^2)'
% instead would be strictly more robust, but would depart from
% buildBoomerangUnitCell.m -- left as-is deliberately.
model.component('comp1').geom('geom1').feature('wp1').geom.create('h_disksel2', 'DiskSelection');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel2').set('entitydim', 0);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel2').set('posx', 'hx');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel2').set('posy', 'hy');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel2').set('r', 'r+selw');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel2').set('rin', 'r-selw');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_disksel2').set('condition', 'allvertices');
model.component('comp1').geom('geom1').feature('wp1').geom.create('h_fil2', 'Fillet');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_fil2').set('radius', 'r2');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('h_fil2').selection('point').named('h_disksel2');

% Extrude the full slab z = [-th/2, +th/2], then subtract the lower half so
% only z = [0, +th/2] is meshed, with the z-parity BC on the z=0 face.
model.component('comp1').geom('geom1').create('ext1', 'Extrude');
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', 'th', 0);
model.component('comp1').geom('geom1').feature('ext1').selection('input').set({'wp1'});
model.component('comp1').geom('geom1').create('symZWP', 'WorkPlane');
model.component('comp1').geom('geom1').feature('symZWP').set('quickz', '-th/2');
model.component('comp1').geom('geom1').feature('symZWP').geom.create('pol2', 'Polygon');
model.component('comp1').geom('geom1').feature('symZWP').geom.feature('pol2').set('source', 'table');
model.component('comp1').geom('geom1').feature('symZWP').geom.feature('pol2').set('table', {'0' '0'; 'a/2' '(a/2)*sqrt(3)'; 'a*(3/2)' '(a/2)*sqrt(3)'; 'a' '0'; '0' '0'});
model.component('comp1').geom('geom1').create('symZPlaneExt', 'Extrude');
model.component('comp1').geom('geom1').feature('symZPlaneExt').setIndex('distance', 'th/2', 0);
model.component('comp1').geom('geom1').feature('symZPlaneExt').selection('input').set({'symZWP'});
model.component('comp1').geom('geom1').create('symZComp', 'Compose');
model.component('comp1').geom('geom1').feature('symZComp').set('formula', 'ext1 - symZPlaneExt');

% z=0 midplane face selection.
% FIXED: the exported box spanned x,y in [-250,+250] nm with 'allvertices',
% but this cell spans x in [0, 3a/2] and y in [0, a*sqrt(3)/2] -- so it
% selected NOTHING. (Same bug in buildBoomerangUnitCell.m:95-96, unnoticed
% because runBands_2D.m uses P.zEnd from bndindex, never P.bndSel.Zsym.)
% Harmless until now since symBCs/asymBCs below use the literal index [3];
% with a correct box this selection can replace that index.
model.component('comp1').geom('geom1').create('ZsymSel', 'BoxSelection');
model.component('comp1').geom('geom1').feature('ZsymSel').set('entitydim', 2);
model.component('comp1').geom('geom1').feature('ZsymSel').set('xmin', '-dsel');
model.component('comp1').geom('geom1').feature('ZsymSel').set('xmax', 'a*(3/2)+dsel');
model.component('comp1').geom('geom1').feature('ZsymSel').set('ymin', '-dsel');
model.component('comp1').geom('geom1').feature('ZsymSel').set('ymax', 'a*sqrt(3)/2+dsel');
model.component('comp1').geom('geom1').feature('ZsymSel').set('zmin', '-dsel');
model.component('comp1').geom('geom1').feature('ZsymSel').set('zmax', 'dsel');
model.component('comp1').geom('geom1').feature('ZsymSel').set('condition', 'allvertices');
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').material.create('diamond', 'Common');
model.component('comp1').material('diamond').propertyGroup.create('AnisotropicVoGrp', 'AnisotropicVoGrp', 'Anisotropic, Voigt notation');

model.component('comp1').common.create('mpf1', 'ParticipationFactors');

model.component('comp1').physics.create('smech', 'SolidMechanics', 'geom1');

% !! FACE INDICES BELOW ARE NOT PARAMETRIC -- see item 1 in the header. These
% are absolute boundary numbers valid ONLY for the default geometry. Any change
% to a/w/r/r1/r2/th can renumber them, after which the Floquet pairs and the
% z-parity plane silently attach to the wrong faces. Replace with geometric
% selections (the ZsymSel above is now correct and can take over for [3]).
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
