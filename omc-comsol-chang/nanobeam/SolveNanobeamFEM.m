% function to solve for optical and/or mechanical modes and postprocess results
%
% Inputs
% model: COMSOL model with studies set up
% ds: data structure with FEM substructures
% 
% Outputs
% model: updated COMSOL model with solutions
% ds: data structure with FEM substructures containing postprocessing results

function [model,ds] = SolveNanobeamFEM(model,ds)

% constants
c = 299792458;
hbar = 1.05457148e-34; % J*s

% extract parameters from P
P = ds.P;
max_dof = P.max_dof;

% center of beam
xc = P.xc;

% symmetry factor
symFac = 2^(abs(P.mevenx)+abs(P.meveny)+abs(P.mevenz*strcmp(P.xsect,'rect')));

% optical
if P.solveOpt
    ofem = ds.ofem;         % optical FEM data structure
    oneigs = P.oneigs;      % no. of optical eigenmodes to solve for
    lambda = P.lambda;      % target optical wavelength
end

% mechanical
if P.solveMech
    mfem = ds.mfem;         % mechanical FEM data structure
    mneigs = P.mneigs;      % no. of mechanical eigenmodes to solve for
    freq = P.freq;          % target mechanical frequency
    bdom = mfem.dia_domind; % domain index for diamond
end

% study, solver, physics tags; geometry name
studyTags = mphmodel(model.study);
solvTags = mphmodel(model.sol);
physTags = mphmodel(model.physics);
geomnames = fieldnames(mphmodel(model.geom));
geomname = geomnames{1};
beam = model.geom(geomname);

%% === Solve for optical eigenmodes of nanobeam ===
if P.solveOpt
% optical eigenfrequency study
if ~isfield(studyTags,'ostudy')
    ostudy = model.study.create('ostudy');
    ostd_eigv = ostudy.create('ostd_eigv','Eigenvalue');
else
    ostudy = model.study('ostudy');
    ostd_eigv = ostudy.feature('ostd_eigv');
end
ostudy.label('Optical Eigenmode Study');

% set no. of eigenmodes to solve for, and target frequency to solve about
ostd_eigv.set('neigsactive',true).set('neigs',oneigs);
ostd_eigv.set('shiftactive',true).set('shift',num2str(0-1i*2*pi*c/lambda));

% disable solid mechanics in EM waves study
% and disable EM waves in solid mechanics study
if isfield(physTags,'smech')
    ostd_eigv.set('useadvanceddisable','on').set('disabledphysics',{'smech'});
end
if isfield(physTags,'emw') && isfield(studyTags,'mstudy')
    mstudy = model.study('mstudy');
    mstd_eigv = mstudy.feature('mstd_eigv');
    mstd_eigv.set('useadvanceddisable','on').set('disabledphysics',{'emw'});
end

% optical eigenfrequency solver
if ~isfield(solvTags,'osolv')
    osolv = model.sol.create('osolv');
    osolv.study('ostudy');           % connect solver sequence to study node
    osolv.attach('ostudy');          % comes from .m saved from GUI - needed?
    osolv_stdstep = osolv.create('osolv_stdstep', 'StudyStep');
    osolv_vars = osolv.create('osolv_vars', 'Variables');
    osolv_eigv = osolv.create('osolv_eigv', 'Eigenvalue');
    model.result.dataset('dset1').tag('odset'); %dset1 is default tag for new dsets
else
    osolv = model.sol('osolv');
    osolv.study('ostudy');           % connect solver sequence to study node
    osolv.attach('ostudy');          % comes from .m saved from GUI - needed?
    osolv_stdstep = osolv.feature('osolv_stdstep');
    osolv_vars = osolv.feature('osolv_vars');
    osolv_eigv = osolv.feature('osolv_eigv');
end
model.result.dataset('odset').set('solution', 'osolv');
model.result.dataset('odset').set('geom', geomname);
model.result.dataset('odset').label('Optical Solutions');
osolv.label('Optical Eigenmode Solver');
osolv_stdstep.set('study','ostudy').set('studystep','ostd_eigv');
osolv_vars.set('control','ostd_eigv');
osolv_eigv.set('control','ostd_eigv');
osolv_eigv.set('eigref',num2str(0-1i*2*pi*c/lambda));
osolv_eigv.feature('dDef').set('linsolver', 'pardiso');
osolv_eigv.feature('aDef').set('complexfun', 'on');

% disable solid mechanics in EM waves study
% and disable EM waves in solid mechanics study
if isfield(physTags,'smech')
    ostd_eigv.set('useadvanceddisable','on').set('disabledphysics',{'smech'});
end
if isfield(physTags,'emw') && isfield(studyTags,'mstudy')
    mstudy = model.study('mstudy');
    mstd_eigv = mstudy.feature('mstd_eigv');
    mstd_eigv.set('useadvanceddisable','on').set('disabledphysics',{'emw'});
end

% Create mesh
meshTags = mphmodel(model.mesh);
if ~isfield(meshTags,'mesh')
    mesh = model.mesh.create('mesh', geomname);
else
    mesh = model.mesh('mesh');
end
% The optical and the mechanical pass share this one sequence but only the
% mechanical pass calls applyMeshOverrides, so the optical pass can inherit a
% sequence that is already user-controlled. No-op unless it is user-controlled
% AND has no meshing operation - see ensureFreeTet for why that state meshes
% nothing at all.
ensureMeshOperation(mesh);

% Solve for eigenfrequencies, adjust mesh if max DOFs exceeded
ok = 0;
if P.oAdjMesh
    
    ofem.mesh = P.oMesh;
    disp(['Meshing with quality: ' num2str(ofem.mesh)]);
    mesh.feature('size').set('custom','off').set('hauto',ofem.mesh);
    mesh.feature('size').set('custom','on').set('hmax',lambda/5);
    oxmesh = mphxmeshinfo(model, 'soltag', 'osolv', ...
                             'studysteptag', 'osolv_stdstep');
    dofs = oxmesh.ndofs;
    disp(['Estimated no. of DoFs, optical: ' num2str(dofs)]);
    ofem.dofs = dofs;
    
    while (ok == 0) && (ofem.mesh < 10)
        while (dofs >= max_dof) && (ofem.mesh < 10)
            ofem.mesh = ofem.mesh + 1;
            disp(['  Meshing with quality: ' num2str(ofem.mesh)]);
            mesh.feature('size').set('custom','off').set('hauto',ofem.mesh);
            mesh.feature('size').set('custom','on').set('hmax',lambda/5);
            mesh.run;
            oxmesh = mphxmeshinfo(model, 'soltag', 'osolv', ...
                                         'studysteptag', 'osolv_stdstep');
            dofs = oxmesh.ndofs;
            disp(['Estimated no. of DoFs, optical: ' num2str(dofs)]);
            ofem.dofs = dofs;
        end 
        pause(2);
        try 
            % solve eigenvalue problem
            osolv.runAll;
            ok = 1;
        catch lasterror
            % if memory runs out while solving even though dofs < maxdof,
            % reduce mesh quality
            disp(lasterror.message)
            disp('ERROR stack:')
            for lsi=1:length(lasterror.stack)
                disp([num2str(lsi),': line ',...
                         num2str(lasterror.stack(lsi).line),...
                         ' in file: ',lasterror.stack(lsi).file])
            end
            ofem.mesh = ofem.mesh + 1;
            disp(['Meshing with quality: ' num2str(ofem.mesh)]);
            mesh.feature('size').set('custom','off').set('hauto',ofem.mesh);
            mesh.feature('size').set('custom','on').set('hmax',lambda/5);
            mesh.run;
            oxmesh = mphxmeshinfo(model, 'soltag', 'osolv', ...
                                         'studysteptag', 'osolv_stdstep');
            dofs = oxmesh.ndofs;
            disp(['Estimated no. of DoFs, optical: ' num2str(dofs)]);
            ofem.dofs = dofs;
        end 
    end %of: while (ok == 0) && (ofem.mesh < 10)
else
    osolv.runAll;
end

% rename tags of datasets
model.result.dataset('odset').set('solution', 'osolv');

%% Optical simulation results
% assemble eigenvectors and eigenfrequencies
% eigenvalue arrays given by COMSOL are column arrays 
% - transpose to get row array
oPValRe = transpose(osolv.getPVals());
oPValIm = transpose(osolv.getPValsImag());
ofem.PVals = oPValRe + 1i*oPValIm;

% postprocessing: wavelengths, frequencies, Q
ofem.lambdasC = 2*pi*c./ofem.PVals;
ofem.lambdaAll = -2*pi*c./oPValIm;
ofem.fOAll = c./ofem.lambdaAll;
ofem.QAll = (-oPValIm)./(2*oPValRe);

% extract results for high-Q optical modes
ofem.hiQSol = find(ofem.QAll>1e4);
ofem.lambda = ofem.lambdaAll(ofem.hiQSol);
ofem.Q = ofem.QAll(ofem.hiQSol);

ds.ofem = ofem;

display('Optical simulations done - Postprocessing done');
if ~isempty(ofem.hiQSol)
    display('  High Q optical mode wavelength and Q:')
    for hiQIdx = 1:length(ofem.hiQSol)
        display(['  mode ',num2str(ofem.hiQSol(hiQIdx)),': lambda = ',num2str(ofem.lambda(hiQIdx)*1e9,'%.2f'),' nm, Q = ',num2str(ofem.Q(hiQIdx),'%.2e')])
    end
else
    error('  No high-Q optical modes found')
end


end % of if P.solveOpt

%% === Solve for mechanical eigenmodes of nanobeam ===
if P.solveMech
% mechanical eigenfrequency study
if ~isfield(studyTags,'mstudy')
    mstudy = model.study.create('mstudy');
    mstd_eigv = mstudy.create('mstd_eigv','Eigenvalue');
else
    mstudy = model.study('mstudy');
    mstd_eigv = mstudy.feature('mstd_eigv');
end
mstudy.label('Mechanical Eigenmode Study');

% set no. of eigenmodes to solve for, and target frequency to solve about
mstd_eigv.set('neigsactive',true).set('neigs',mneigs);
mstd_eigv.set('shiftactive',true).set('shift',num2str(0-1i*2*pi*freq));

% disable EM waves in solid mechanics study
% and disable solid mechanics in EM waves study
if isfield(physTags,'emw')
    mstd_eigv.set('useadvanceddisable','on').set('disabledphysics',{'emw'});
end
if isfield(physTags,'smech') && isfield(studyTags,'ostudy')
    ostudy = model.study('ostudy');
    ostd_eigv = ostudy.feature('ostd_eigv');
    ostd_eigv.set('useadvanceddisable','on').set('disabledphysics',{'smech'});
end

% mechanical eigenfrequency solver
if ~isfield(solvTags,'msolv')
    msolv = model.sol.create('msolv');
    msolv.study('mstudy');           % connect solver sequence to study node
    msolv.attach('mstudy');          % comes from .m saved from GUI - needed?
    msolv_stdstep = msolv.create('msolv_stdstep', 'StudyStep');
    msolv_vars = msolv.create('msolv_vars', 'Variables');
    msolv_eigv = msolv.create('msolv_eigv', 'Eigenvalue');
    model.result.dataset('dset1').tag('mdset'); %dset1 is default tag for new dsets
else
    msolv = model.sol('msolv');
    msolv.study('mstudy');           % connect solver sequence to study node
    msolv.attach('mstudy');          % comes from .m saved from GUI - needed?
    msolv_stdstep = msolv.feature('msolv_stdstep');
    msolv_vars = msolv.feature('msolv_vars');
    msolv_eigv = msolv.feature('msolv_eigv');
end
model.result.dataset('mdset').set('solution', 'msolv');
model.result.dataset('mdset').set('geom', geomname);
model.result.dataset('mdset').label('Mechanical Solutions');
msolv.label('Mechanical Eigenmode Solver');
msolv_stdstep.set('study','mstudy').set('studystep','mstd_eigv');
msolv_vars.set('control','mstd_eigv');
msolv_eigv.set('control','mstd_eigv');
msolv_eigv.set('eigref',num2str(0-1i*2*pi*freq));
msolv_eigv.feature('dDef').set('linsolver', 'mumps');
% msolv_eigv.feature('aDef').set('complexfun', 'off');
msolv_eigv.feature('aDef').set('complexfun', 'on');


% disable EM waves in solid mechanics study
% and disable solid mechanics in EM waves study
if isfield(physTags,'emw')
    mstd_eigv.set('useadvanceddisable','on').set('disabledphysics',{'emw'});
end
if isfield(physTags,'smech') && isfield(studyTags,'ostudy')
    ostudy = model.study('ostudy');
    ostd_eigv = ostudy.feature('ostd_eigv');
    ostd_eigv.set('useadvanceddisable','on').set('disabledphysics',{'smech'});
end

% Create mesh
meshTags = mphmodel(model.mesh);
if ~isfield(meshTags,'mesh')
    mesh = model.mesh.create('mesh', geomname);
else
    mesh = model.mesh('mesh');
end

ensureMeshOperation(mesh);   % see the note in the optical pass above

% Domain-scoped mesh overrides for the PML and (if present) the phononic
% shield are applied by applyMeshOverrides below. They are Size nodes nested
% inside the sequence's FreeTet operation and override the global 'size' node
% on their own domains only; changing the global size does NOT remove the
% override, but it is re-asserted after every global change so the intent
% stays visible in the sequence.

% Solve for eigenfrequencies, adjust mesh if max DOFs exceeded
ok = 0;
if P.mAdjMesh

    mfem.mesh = P.mMesh;
    disp(['Meshing with quality: ' num2str(mfem.mesh)]);
	mesh.feature('size').set('custom','off').set('hauto',mfem.mesh);
    mesh.run;   % applied to ALL domains including PML
    % Create the domain-scoped size nodes once (reuse them if the mesh was
    % carried over from a previous run), before the first DOF estimate.
    applyMeshOverrides(mesh, geomname, P);
    mxmesh = mphxmeshinfo(model, 'soltag', 'msolv', ...
                             'studysteptag', 'msolv_stdstep');
    dofs = mxmesh.ndofs;
    disp(['Estimated no. of DoFs, mechanical: ' num2str(dofs)]);
    mfem.dofs = dofs;
    mesh.run;   % applied to ALL domains including PML
    while (ok == 0) && (mfem.mesh < 10)
        while (dofs >= max_dof) && (mfem.mesh < 10)
            mfem.mesh = mfem.mesh + 1;
            disp(['Meshing with quality: ' num2str(mfem.mesh)]);
            mesh.feature('size').set('custom','off').set('hauto',mfem.mesh);
            % Re-assert the domain overrides; the global change above leaves
            % them in place, but set them explicitly to keep intent clear.
            applyMeshOverrides(mesh, geomname, P);

            mesh.run;
            mxmesh = mphxmeshinfo(model, 'soltag', 'msolv', ...
                                         'studysteptag', 'msolv_stdstep');
            dofs = mxmesh.ndofs;
            disp(['Estimated no. of DoFs, mechanical: ' num2str(dofs)]);
            mfem.dofs = dofs;
        end 
        pause(2);
        try 
            % solve eigenvalue problem
            msolv.runAll;
            ok = 1;
            %display('Solving done')
        catch lasterror
            % if memory runs out while solving even though dofs < maxdof,
            % reduce mesh quality
            disp(lasterror.message)
            disp('ERROR stack:')
            for lsi=1:length(lasterror.stack)
                display([num2str(lsi),': line ',...
                         num2str(lasterror.stack(lsi).line),...
                         ' in file: ',lasterror.stack(lsi).file])
            end
            mfem.mesh = mfem.mesh + 1;
            disp(['Meshing with quality: ' num2str(mfem.mesh)]);
            mesh.feature('size').set('custom','off').set('hauto',mfem.mesh);
            %mesh.feature('size').set('custom','on').set('hmax',lambda/5);
            applyMeshOverrides(mesh, geomname, P);
            mesh.run;
            mxmesh = mphxmeshinfo(model, 'soltag', 'msolv', ...
                                         'studysteptag', 'msolv_stdstep');
            dofs = mxmesh.ndofs;
            disp(['Estimated no. of DoFs, mechanical: ' num2str(dofs)]);
            mfem.dofs = dofs;
        end 
    end %of: while (ok == 0) && (mfem.mesh < 10)
else
    disp(['Meshing with quality: ' num2str(P.mMesh)]);
    mesh.feature('size').set('custom','off').set('hauto',P.mMesh);
    % Domain-scoped size nodes (created once, reused if they already exist) so
    % the PML and the phononic shield can be meshed differently to the beam.
    applyMeshOverrides(mesh, geomname, P);
    mesh.run;
    msolv.runAll;
end

% rename tags of datasets
model.result.dataset('mdset').label('Mechanical Solutions');
model.result.dataset('mdset').set('solution', 'msolv');

%% Mechanical simulation results
% eigenvalue arrays given by COMSOL are column arrays 
% - transpose to get row array
mPValRe = transpose(msolv.getPVals());
mPValIm = transpose(msolv.getPValsImag());
mfem.PVals = mPValRe + 1i*mPValIm;
mfem.freqsC = mfem.PVals/(2*pi);
mfem.freqs = abs(mPValIm)/(2*pi);   % mechanical frequencies

% mechanical quality factor (radiation loss via PML)
if isfield(P,'solveMechPML') && P.solveMechPML
    mfem.QAll = abs(mPValIm) ./ (2*abs(mPValRe));
end

% extract results for localized mechanical modes
% by calculating ratio of integrated displacements in center of beam to
% that of whole beam
extractLocMechModes = P.extractLocMechModes;
if extractLocMechModes
    
    maxX = 2*P.a;
    maxY = P.w;
    stepWeight = ['if(((abs(x-',num2str(xc),')<',num2str(maxX),')',...
                    '&&(abs(y)<',num2str(maxY),')',...
                    '),1,0)'];
    totDisp = mphint2(model,'solid.disp','volume',...
                      'dataset','mdset','selection',bdom,'solnum','all');
    locDisp = mphint2(model,[stepWeight,'*solid.disp'],'volume',...
                      'dataset','mdset','selection',bdom,'solnum','all');
    mfem.locRatio = locDisp./totDisp;
	mfem.locInd = find(locDisp./totDisp>0.1);
    mfem.locFreqs = mfem.freqs(mfem.locInd);
else
    mfem.locRatio = zeros(1,length(mfem.freqs));
    mfem.locFreqs = mfem.freqs;
    mfem.locInd = 1:P.mneigs;
end

% effective mass and zero-point displacement
dispSqInt = mphint2(model,'solid.disp^2','volume','dataset','mdset','selection',bdom,'solnum','all');
dispSqMax = mphmax(model,'solid.disp^2','volume','dataset','mdset','selection',bdom,'solnum','all');
mfem.meff = symFac*P.rho*(dispSqInt./dispSqMax);
mfem.xzpf = sqrt(hbar./(2.*mfem.meff.*2.*pi.*mfem.freqs));

ds.mfem = mfem;

disp('Mechanical simulations done - Postprocessing done');
if ~isempty(mfem.locInd)
    disp('  Localized mechanical modes:')
    for locIdx = 1:length(mfem.locInd)
        mIdx = mfem.locInd(locIdx);
        if isfield(mfem,'QAll')
            disp(['  mode ',num2str(mIdx),': wM = ',num2str(mfem.locFreqs(locIdx)/1e9,'%.2f'),' GHz, Qrad = ',num2str(mfem.QAll(mIdx),'%.2e')])
        else
            disp(['  mode ',num2str(mIdx),': wM = ',num2str(mfem.locFreqs(locIdx)/1e9,'%.2f'),' GHz'])
        end
    end
else
    disp('  No localized mechanical modes found, saving most localized mode')
    [~,mfem.locInd] = max(locDisp./totDisp);
    mfem.locFreqs = mfem.freqs(mfem.locInd);
    ds.mfem = mfem;
end


    
end % of if P.solveMech

% remove extra datasets
% mphmodel(model.result.dataset)
Dnames = mphmodel(model.result.dataset);
for i=1:5
    dsetname = ['dset',num2str(i)];
    if isfield(Dnames,dsetname)
        model.result.dataset.remove(dsetname);
    end
end
% % mphmodel(model.result.dataset)

end

% -------------------------------------------------------------------------

function applyMeshOverrides(mesh, geomname, P)
%APPLYMESHOVERRIDES Domain-scoped Size nodes for the PML and the phononic shield.
%
% Every branch is field-guarded. A caller that sets only P.PMLmesh gets exactly
% the 'size_pml' hauto override this file has always applied, on exactly the
% same domains.
%
% P.PMLmeshDiv and P.shieldHmax are opt-in additions:
%
%   P.PMLmeshDiv  number of elements across the absorbing direction. hmax =
%                 P.PMLLen/P.PMLmeshDiv; 8 is the usual minimum. Takes
%                 precedence over P.PMLmesh when both are set, because a PML
%                 needs a resolution tied to its own thickness, not a global
%                 quality level.
%                 GATED to P.celltype = 'crossShield'. This field was
%                 previously read by nothing at all, and
%                 test_nanobeamRectFEM_withPML.m:125-126 already sets it
%                 (PMLmeshDiv = 20) with P.PMLmesh commented out. Honouring
%                 it for every caller would silently move that script's PML
%                 from the global mesh quality to hmax = 500 nm and shift its
%                 results. Legacy callers keep the P.PMLmesh path untouched.
%
%   P.shieldHmax  max element size in P.domSel.shield [m]. The cross-shield
%                 ligaments are the narrowest solid feature in the model and
%                 carry the whole shield response, so they need their own
%                 resolution independent of the global quality level.
%
%   P.beamHmax    max element size in the beam (P.domSel.beam minus the
%                 shield) [m]. Default min(P.th, min(hx))/3. GATED to
%                 P.celltype = 'crossShield'. Without it the beam is the only
%                 domain the mAdjMesh coarsening loop can touch - see the long
%                 comment at the node itself for why that loop cannot converge
%                 once the shield carries its own hmax.
%
% WHY THE DECISIONS ARE ALL MADE BEFORE ANY NODE IS CREATED
% A mesh sequence with nothing in it but the global 'size' node is
% PHYSICS-CONTROLLED and meshes itself. Creating even one extra feature flips
% it to user-controlled, and from then on it runs only what it contains - see
% ensureFreeTet. So this function must not touch the sequence at all unless it
% is actually going to override something, and when it does it must add the
% FreeTet operation as well as the Size nodes. Hence: work out every condition
% first, bail out if none of them fired, then create.

isCrossShield = isfield(P,'celltype') && strcmp(P.celltype,'crossShield');
hasDomSel     = isfield(P,'domSel');

% --- PML -----------------------------------------------------------------
% See the P.PMLmeshDiv note in the header: gated to the cross-shield path so
% that pre-existing callers which set the field are unaffected.
doPML     = false;
usePMLdiv = false;
if isfield(P,'solveMechPML') && P.solveMechPML && ...
        hasDomSel && isfield(P.domSel,'PML') && ~isempty(P.domSel.PML)
    usePMLdiv = isCrossShield && ...
                isfield(P,'PMLmeshDiv') && ~isempty(P.PMLmeshDiv) && ...
                P.PMLmeshDiv > 0 && isfield(P,'PMLLen');
    doPML = usePMLdiv || isfield(P,'PMLmesh');
end

% --- Beam ----------------------------------------------------------------
% P.domSel.beam is beam + pad + shield; the shield has its own node, so take
% the difference rather than relying on Size-node ordering.
doBeam   = false;
beamOnly = [];
if isCrossShield && hasDomSel && ...
        isfield(P.domSel,'beam') && ~isempty(P.domSel.beam)
    beamOnly = P.domSel.beam;
    if isfield(P.domSel,'shield') && ~isempty(P.domSel.shield)
        beamOnly = setdiff(beamOnly, P.domSel.shield);
    end
    if isempty(beamOnly)
        warning('SolveNanobeamFEM:noBeamDomain', ...
            ['P.domSel.beam minus P.domSel.shield is empty, so the beam got ' ...
             'no mesh override and the mAdjMesh loop can coarsen it away. ' ...
             'Check the builder''s domain selections.']);
    else
        doBeam = true;
    end
end

% --- Shield --------------------------------------------------------------
doShield = isfield(P,'shieldHmax') && ~isempty(P.shieldHmax) && ...
           hasDomSel && isfield(P.domSel,'shield') && ~isempty(P.domSel.shield);

% Nothing to override: leave the sequence physics-controlled. This is the path
% every legacy caller takes (test_nanobeamRectFEM.m has P.solveMechPML = 0;
% test_nanobeamRectFEM_withPML.m sets only P.PMLmeshDiv, which is gated to
% crossShield), and it must stay byte-for-byte the mesh they have always had.
if ~(doPML || doBeam || doShield)
    return;
end

% From here on the sequence is user-controlled, so it needs a real meshing
% operation and the Size nodes must be nested inside it.
ftet = ensureFreeTet(mesh);

if doPML
    szPML = getOrCreateSizeNode(ftet, 'size_pml');
    szPML.selection.geom(geomname, 3).set(P.domSel.PML);
    if usePMLdiv
        szPML.set('custom','on').set('hmaxactive',true);
        szPML.set('hmax', P.PMLLen/P.PMLmeshDiv);
        disp(['PML mesh: hmax = PMLLen/',num2str(P.PMLmeshDiv),' = ', ...
              num2str(P.PMLLen/P.PMLmeshDiv*1e9,'%.0f'),' nm']);
    else
        szPML.set('custom','off').set('hauto', P.PMLmesh);
        disp(['PML mesh quality: ' num2str(P.PMLmesh)]);
    end
end

% Beam. WHY THIS EXISTS: the mAdjMesh loop above raises the GLOBAL hauto while
% dofs >= max_dof, capped only by mfem.mesh < 10, i.e. up to 9 - "extremely
% coarse". The shield and the PML each carry their own hmax and are immune to
% that, so the shield's contribution - which dominates the DOF count - never
% shrinks and the loop cannot converge. It therefore runs hauto all the way up,
% and the beam, the one part of the model with no override, is the only thing
% that actually gets coarsened. At hauto = 9 the global element size is set by
% an ~18 um model against a 250 nm slab, so the beam ends up unmeshed or
% degenerate while the shield and PML still look perfectly meshed.
%
% Pinning the beam makes it independent of that loop: an undersized max_dof can
% still stop the solve, but it can no longer silently destroy the cavity mesh,
% which is the part the physics of interest lives in.
%
% Gated to crossShield for the same reason P.PMLmeshDiv is (see the header):
% legacy callers reach this function with no beam override and must keep the
% mesh they have always had.
if doBeam
    if isfield(P,'beamHmax') && ~isempty(P.beamHmax)
        beamHmax = P.beamHmax;
    else
        % Three elements across the narrowest beam feature, matching the
        % convention P.shieldHmax uses for the ligaments. The slab
        % thickness and the smallest hole dimension are the candidates.
        narrow = P.th;
        if isfield(P,'geomHalf') && ~isempty(P.geomHalf)
            narrow = min(narrow, min(P.geomHalf(:,1)));
        end
        beamHmax = narrow/3;
    end
    szB = getOrCreateSizeNode(ftet, 'size_beam');
    szB.selection.geom(geomname, 3).set(beamOnly);
    szB.set('custom','on').set('hmaxactive',true);
    szB.set('hmax', beamHmax);
    disp(['Beam mesh: hmax = ',num2str(beamHmax*1e9,'%.1f'),' nm on ', ...
          num2str(numel(beamOnly)),' domain(s)']);
end

if doShield
    szS = getOrCreateSizeNode(ftet, 'size_shield');
    szS.selection.geom(geomname, 3).set(P.domSel.shield);
    szS.set('custom','on').set('hmaxactive',true);
    szS.set('hmax', P.shieldHmax);
    disp(['Shield mesh: hmax = ',num2str(P.shieldHmax*1e9,'%.1f'),' nm']);
end
end

% -------------------------------------------------------------------------

function ensureMeshOperation(mesh)
%ENSUREMESHOPERATION Repair a user-controlled sequence that meshes nothing.
%
% No-op in both of the states this code normally produces:
%   {size}              physics-controlled, COMSOL generates its own operations
%   {size, ftet_all, …} already carries an explicit meshing operation
% It only acts on {size, size_pml, …} with no operation - the state a model
% built by an earlier version of this file is left in. Detected by inspecting
% tags rather than by calling mesh.isAutomatic(), so it does not depend on that
% accessor being present in the installed COMSOL release.
tags = cellstr(char(mesh.feature.tags()));
if any(~strcmp(tags, 'size')) && ~any(strcmp(tags, 'ftet_all'))
    ensureFreeTet(mesh);
end
end

% -------------------------------------------------------------------------

function ftet = ensureFreeTet(mesh)
%ENSUREFREETET Guarantee the mesh sequence contains a meshing operation.
%
% WHY THIS IS NEEDED. model.mesh.create(tag, geom) returns a sequence holding a
% single global 'size' node, with isAutomatic = 1: it is PHYSICS-CONTROLLED and
% COMSOL generates the meshing operations itself at run time. Writing to 'size'
% keeps it that way - which is why setting hauto has always worked. Creating any
% OTHER feature, including a domain-scoped Size node, flips the sequence to
% user-controlled, and a user-controlled sequence runs only the features it
% actually holds. A Size node is not a meshing operation, so the run then
% produces no elements anywhere and COMSOL reports
%   "No mesh on domains 1-N in the meshing sequence with tag mesh"
% for every domain in the geometry - not for a subset, which is the tell that
% this is a sequence fault and not a selection fault.
%
% Measured on COMSOL 6.3, LiveLink, two-domain test block, hauto = 5:
%   'size' only                   isAutomatic = 1    4877 elements
%   'size' + a top-level Size     isAutomatic = 0       0 elements
%   + an explicit FreeTet         isAutomatic = 0    3446 elements, hmax honoured
%
% The FreeTet's selection is deliberately left at its default ("Remaining" -
% every domain no earlier operation has meshed). With this as the sole
% operation that is all of them, and it stays correct if a swept PML mesh is
% ever inserted ahead of it.
%
% Same construction the band-structure pipeline already uses:
% bands/optmechbands1D/RunOpticalBands.m:370-374, Size nested inside FreeTet.
tag  = 'ftet_all';
tags = cellstr(char(mesh.feature.tags()));
if any(strcmp(tag, tags))
    ftet = mesh.feature(tag);
    return;
end
% Discard top-level Size overrides left by an earlier version of this file.
% They are re-created nested inside the FreeTet; left where they are they would
% sit AFTER the operation in the sequence and so would not apply to it.
for stale = {'size_pml','size_shield','size_beam'}
    if any(strcmp(stale{1}, tags))
        mesh.feature.remove(stale{1});
    end
end
ftet = mesh.feature.create(tag, 'FreeTet');
ftet.label('Free Tetrahedral - all domains');
disp('Mesh: added explicit FreeTet (sequence is user-controlled once a Size override exists).');
end

% -------------------------------------------------------------------------

function sz = getOrCreateSizeNode(parent, tag)
%GETORCREATESIZENODE Fetch a Size feature under PARENT, creating it on first use.
% PARENT is the FreeTet operation, NOT the mesh sequence: a Size node nested
% inside the operation always runs before it, whereas a top-level Size created
% after the operation would never apply to it.
if any(strcmp(tag, cellstr(char(parent.feature.tags()))))
    sz = parent.feature(tag);
else
    sz = parent.create(tag, 'Size');
end
end