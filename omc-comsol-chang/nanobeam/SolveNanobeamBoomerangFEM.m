% function to solve for optical and/or mechanical modes and postprocess results
%
% Inputs
% model: COMSOL model with studies set up
% ds: data structure with FEM substructures
% 
% Outputs
% model: updated COMSOL model with solutions
% ds: data structure with FEM substructures containing postprocessing results

function [model,ds] = SolveNanobeamBoomerangFEM(model,ds)

% constants
c = 299792458;
hbar = 1.05457148e-34; % J*s

% extract parameters from P
P = ds.P;
max_dof = P.max_dof;

% center of beam
% xc = P.xc;
xc = 0;

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
    % mstd_eigv.set('shiftactive',true).set('shift',num2str(2*pi*freq));

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
    msolv_eigv.feature('aDef').set('complexfun', 'off');

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

    % Solve for eigenfrequencies, adjust mesh if max DOFs exceeded
    ok = 0;
    if P.mAdjMesh

        mfem.mesh = P.mMesh;
        disp(['Meshing with quality: ' num2str(mfem.mesh)]);
        mesh.feature('size').set('custom','off').set('hauto',mfem.mesh);
        mxmesh = mphxmeshinfo(model, 'soltag', 'msolv', ...
                                 'studysteptag', 'msolv_stdstep');
        dofs = mxmesh.ndofs;
        disp(['Estimated no. of DoFs, mechanical: ' num2str(dofs)]);
        mfem.dofs = dofs;
        while (ok == 0) && (mfem.mesh < 10)
            while (dofs >= max_dof) && (mfem.mesh < 10)
                mfem.mesh = mfem.mesh + 1;
                disp(['Meshing with quality: ' num2str(mfem.mesh)]);
                mesh.feature('size').set('custom','off').set('hauto',mfem.mesh);

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
                mesh.run;
                mxmesh = mphxmeshinfo(model, 'soltag', 'msolv', ...
                                             'studysteptag', 'msolv_stdstep');
                dofs = mxmesh.ndofs;
                disp(['Estimated no. of DoFs, mechanical: ' num2str(dofs)]);
                mfem.dofs = dofs;
            end 
        end %of: while (ok == 0) && (mfem.mesh < 10)
    else
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

    % % mechanical quality factor
    % mfem.post.mechQ = imag(mfem.sol.lambda)./(2*real(mfem.sol.lambda));

    % extract results for localized mechanical modes
    % by calculating ratio of integrated displacements in center of beam to
    % that of whole beam
    extractLocMechModes = 1;
    if extractLocMechModes
        
        % defines the boundaries of confined modes 
        maxX = (2*P.TN-1)*P.a;
        maxY = P.b;

        stepWeight = ['if(((abs(x-',num2str(xc),')<',num2str(maxX),')',...
                        '&&(abs(y)<',num2str(maxY),')',...
                        '),1,0)'];
        totDisp = mphint2(model,'solid.disp','volume',...
                          'dataset','mdset','selection',bdom,'solnum','all');
        locDisp = mphint2(model,[stepWeight,'*solid.disp'],'volume',...
                          'dataset','mdset','selection',bdom,'solnum','all');
        mfem.locRatio = locDisp./totDisp;
        mfem.locInd = find(locDisp./totDisp>0.35);
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
            disp(['  mode ',num2str(mfem.locInd(locIdx)),': wM = ',num2str(mfem.locFreqs(locIdx)/1e9,'%.2f'),' GHz'])
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