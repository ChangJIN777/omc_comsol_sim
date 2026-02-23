% function to calculate optomechanical coupling from optical and mechanical FEM simulations of nanobeam
% Coupling calculated between optical/mechanical modes of interest 
% (as specified in oModes/mModes) by adding contributions from 
% moving boundaries (MB/Bnd) and photoelasticity (PE/Str).
% 
% Inputs
% ds: data structure with optical and mechanical FEM substructures (ds.ofem and ds.mfem)
% model: COMSOL model with optical and mechanical FEM solutions
% oModes: optical solution nos. to calculate OM coupling for (as indexed for all solutions, including high- and low-Q modes)
% mModes: mechanical solution nos. to calculate OM coupling for (as indexed for all solutions, including localized and unlocalized modes)
% 
% Outputs
% ds: updated data structure with additional coupling substructure ds.cpl
% model: updated COMSOL model with solutions

function [ds,model] = CalcGOM(ds,model,oModes,mModes)

if ~isfield(ds,'cpl')
    cpl = [];
else
    cpl = ds.cpl;
end

% check both optical and mechanical FEM data structs exist
if ~isfield(ds,'ofem') || ~isfield(ds,'mfem')
    error('ERROR: missing ofem or mfem in data struct')
end
ofem = ds.ofem;
mfem = ds.mfem;
P = ds.P;

% constants
c = 299792458; % m/s
hbar = 1.05457148e-34; % J*s

% init gMax to search for max g over for loops
cpl.gMax = 0;

rxtal = P.rxtal;

% material constants
n = ofem.n{end};
rho = mfem.rho;
Dn2 = n^2 - 1;
Dn2_1 = 1/n^2 - 1;

% symmetry factor
symFac = 2^(abs(P.mevenx)+abs(P.meveny)+abs(P.mevenz*strcmp(P.xsect,'rect')));

% generate signs for each octant of structure
% for summing components of coupling while taking symmetry into account
% 2nd element of each vector corresponds to symmetry about that axis
XsymVec = [1 P.mevenx];
YsymVec = [1 P.meveny];
ZsymVec = [1 P.mevenz];
% create 3x8 array, where each column contains x/y/z symmetry of structure
[XV,YV,ZV] = meshgrid(XsymVec,YsymVec,ZsymVec);
symAll = transpose([XV(:),YV(:),ZV(:)]); 
% use prod to get sign for each octant
sgnAll = prod(symAll);
sgnCpl = sum(sgnAll);

disp(' ')
disp('Calculating optomechanical coupling...');

%% Domains and boundaries
% choose domain corresponding to beam
bdomM = mfem.dia_domind;
adomO = ofem.air_domind;
bdomO = ofem.dia_domind;

% geometry name
geomnames = fieldnames(mphmodel(model.geom));
geomname = geomnames{1};
beam = model.geom(geomname);

% get boundaries on beam, not overlapping with symmetry planes
model.geom(P.geomname).create('beamBndsAll','AdjacentSelection');
model.geom(P.geomname).feature('beamBndsAll').set('entitydim',3);
model.geom(P.geomname).feature('beamBndsAll').set('input','beamSel');
model.geom(P.geomname).feature('beamBndsAll').set('outputdim',2);
model.geom(P.geomname).feature('beamBndsAll').set('selkeep','on');
beam.runCurrent;
bndsM = model.selection([P.geomname,'_beamBndsAll']).inputEntities(); % get output entities

% exclude symmetry planes
if (abs(P.mevenx) && abs(P.oevenx))
    % bndsM = setdiff(bndsM,P.bndSel.Xsym_l);
    bndsM = setdiff(bndsM,P.bndSel.cylXsym);
end

if (abs(P.meveny) && abs(P.oeveny))
    bndsM = setdiff(bndsM,P.bndSel.cylYsym);
end

if (abs(P.mevenz) && abs(P.oevenz))
    bndsM = setdiff(bndsM,P.bndSel.cylZsym);
end

%% Form datasets for volume and boundary integrals
% if required datasets not present, create datasets
% first: duplicate datasets, set geometry dimension
% then join datasets: first optical, second mechanical
% expressions called using data(...) or data2(...) will apply to
% optical or mechanical solutions
dsetTags = mphmodel(model.result.dataset);

% datasets for boundary integrals
if ~isfield(dsetTags,'jdset_bnd')
    if ~isfield(dsetTags,'mdset_bnd')
        mdset_bnd = model.result.dataset.duplicate('mdset_bnd', 'mdset');
    else
        mdset_bnd = model.result.dataset('mdset_bnd');
    end
    mdset_bnd.label('Mechanical Solutions (Boundary)');
    mdset_bnd.selection.geom(P.geomname, 2).set(bndsM);
    if ~isfield(dsetTags,'odset_bnd')
        odset_bnd = model.result.dataset.duplicate('odset_bnd', 'odset');
    else
        odset_bnd = model.result.dataset('odset_bnd');
    end
    odset_bnd.label('Optical Solutions (Boundary)');
    odset_bnd.selection.geom(P.geomname, 2).set(bndsM);
    jdset_bnd = model.result.dataset.create('jdset_bnd', 'Join');
    jdset_bnd.label('OM Solutions (Boundary)');
    jdset_bnd.set('data', 'odset_bnd').set('solutions', 'one');
    jdset_bnd.set('data2', 'mdset_bnd').set('solutions2', 'one');
    jdset_bnd.set('method', 'explicit');
else
    jdset_bnd = model.result.dataset('jdset_bnd');
end

% datasets for volume integrals
if ~isfield(dsetTags,'jdset_vol')
    if ~isfield(dsetTags,'mdset_vol')
        mdset_vol = model.result.dataset.duplicate('mdset_vol', 'mdset');
    else
        mdset_vol = model.result.dataset('mdset_vol');
    end
    mdset_vol.label('Mechanical Solutions (Volume)');
    mdset_vol.selection.geom(P.geomname, 3).set(bdomM);
    if ~isfield(dsetTags,'odset_vol')
        odset_vol = model.result.dataset.duplicate('odset_vol', 'odset');
    else
        odset_vol = model.result.dataset('odset_vol');
    end
    odset_vol.label('Optical Solutions (Volume)');
    odset_vol.selection.geom(P.geomname, 3).set(bdomO);
    jdset_vol = model.result.dataset.create('jdset_vol', 'Join');
    jdset_vol.label('OM Solutions (Volume)');
    jdset_vol.set('data', 'odset_vol').set('solutions', 'one');
    jdset_vol.set('data2', 'mdset_vol').set('solutions2', 'one');
    jdset_vol.set('method', 'explicit');
else
    jdset_vol = model.result.dataset('jdset_vol');
end

%% Optical normalization and effective volume
lambdaAll = ofem.lambdaAll;
wO = c./lambdaAll;
epsE2Str = 'abs(emw.normE)^2*emw.epsrAv*epsilon0_const';
LV = symFac*mphint2(model,epsE2Str,'volume','dataset','odset','selection','all','solnum',oModes);
LVmax = mphmax(model,epsE2Str,'volume','dataset','odset','selection','all','solnum',oModes);
cpl.Veff = (LV./LVmax)./((lambdaAll(oModes)./(2*n)).^3);   %norm to cubic eff wavelength

%% Displacement normalization and zero-point displacement
cpl.xzpf = mfem.xzpf;
maxDisp = mphmax(model,'abs(solid.disp)','volume','dataset','mdset','selection','all','solnum',mModes);

%% Calculate g
% initialize expressions for moving boundary integrand
mDispExpr = 'data2(u*nX + v*nY + w*nZ)';
oEtExpr = ['(',num2str(Dn2),')*epsilon0_const*data1(',...
            'abs(emw.normE)^2-abs(nx*emw.Ex+ny*emw.Ey+nz*emw.Ez)^2)'];
oDnExpr = ['(',num2str(Dn2_1),')/epsilon0_const*data1(',...
            'abs(nx*emw.Dx)^2+abs(ny*emw.Dy)^2+abs(nz*emw.Dz)^2)'];
MB = [mDispExpr,'*(',oEtExpr,'-',oDnExpr,')'];

% strain components, Voigt notation
S = cell(6,1);
S{1,1} = 'data2(solid.eXX)';
S{2,1} = 'data2(solid.eYY)';
S{3,1} = 'data2(solid.eZZ)';
S{4,1} = 'data2(2*solid.eYZ)';
S{5,1} = 'data2(2*solid.eXZ)';
S{6,1} = 'data2(2*solid.eXY)';

% electric field components
Ex2 = 'data1(epsilon0_const*abs(emw.Ex)^2)';
Ey2 = 'data1(epsilon0_const*abs(emw.Ey)^2)';
Ez2 = 'data1(epsilon0_const*abs(emw.Ez)^2)';
Eyz = 'data1(epsilon0_const*real(emw.Ey*conj(emw.Ez)))';
Exz = 'data1(epsilon0_const*real(emw.Ex*conj(emw.Ez)))';
Exy = 'data1(epsilon0_const*real(emw.Ex*conj(emw.Ey)))';

pcompts = [P.p11,    0,    0;...
               0,P.p12,    0;...
               0,    0,P.p44];%...
%            P.p11,P.p12,P.p44];
[nR,~] = size(pcompts);

% loop over optical modes
oIdx = 1;
for oi = oModes
    disp(['  Processing optical mode ',num2str(oi),'...']);
    jdset_bnd.set('solnum', oi);
    jdset_vol.set('solnum', oi);
    
    % loop over mechanical modes
    mIdx = 1;
    for mi = mModes
        disp(['  Processing mechanical mode ',num2str(mi),'...']);
        jdset_bnd.set('solnum2', mi);   % select solution 
        jdset_vol.set('solnum2', mi);
        wM = mfem.freqs(mi);    %mech freq
        
        %% moving boundary integral
        LMB = mphint2(model,MB,'surface','dataset','jdset_bnd','selection','all');
        cpl.LMB(oi,mi) = LMB;
        cpl.gMB(oi,mi) = -sgnCpl*cpl.xzpf(mi)*0.5*wO(oi)*(LMB/LV(oIdx)/maxDisp(mIdx));
        
        %% photoelastic integral
        % loop over photoelastic components
        for pci = 1:nR
            p11 = pcompts(pci,1);
            p12 = pcompts(pci,2);
            p44 = pcompts(pci,3);
            
            % Compile and rotate photoelastic tensor for cubic crystal
            p = [p11,p12,p12,  0,  0,  0;
                 p12,p11,p12,  0,  0,  0;
                 p12,p12,p11,  0,  0,  0;
                   0,  0,  0,p44,  0,  0;
                   0,  0,  0,  0,p44,  0;
                   0,  0,  0,  0,  0,p44];
            [pR,~] = RotateXtalTensor(p,rxtal);
            
            % generate text expressions for pS elements for evaluation in
            % COMSOL
            % corresponding to Jasper Chan's thesis
            % - only write non-zero terms; special formatting for first term;
            % additional brackets in string only if expression has more than one term
            pS = cell(6,1);
            for i = 1:6
                isFirstTerm = 1;
                nTermCnt = 0;
                for j = 1:6
                    if pR(i,j) ~= 0
                        if isFirstTerm
                            pS{i,1} = ['(',num2str(pR(i,j)),'*',S{j,1},')'];
                            isFirstTerm = 0;
                        else
                            pS{i,1} = [pS{i,1},' + (',num2str(pR(i,j)),'*',S{j,1},')'];
                        end
                        nTermCnt = nTermCnt + 1;
                    end
                end

                %bracket term
                if nTermCnt > 1
                    pS{i,1} = ['( ',pS{i,1},' )'];
                elseif nTermCnt == 0
                    pS{i,1} = '0';
                end
            end

            % extract pS elements
            pS1 = pS{1}; pS2 = pS{2}; pS3 = pS{3};
            pS4 = pS{4}; pS5 = pS{5}; pS6 = pS{6};
            
            % generate expression for phodattoelastic integrand E.pS.E
            EpSEDiv =   [pS1,'*',Ex2,'+',pS2,'*',Ey2,'+',pS3,'*',Ez2];
            EpSEShear = [ '2*',pS4,'*',Eyz,'+2*',pS5,'*',Exz,...
                         '+2*',pS6,'*',Exy];
            % photoelastic integrals
            LPEDc(oi,mi,pci) = mphint2(model,EpSEDiv,  'volume','dataset','jdset_vol','selection','all');
            LPESc(oi,mi,pci) = mphint2(model,EpSEShear,'volume','dataset','jdset_vol','selection','all');
            cpl.gPEc(oi,mi,pci) = sgnCpl*cpl.xzpf(mi)*0.5*wO(oi)*(n^4*(LPEDc(oi,mi,pci)+LPESc(oi,mi,pci))/LV(oIdx)/maxDisp(mIdx));
            
        end
        cpl.LPED(oi,mi) = sum(LPEDc(oi,mi,1:3));
        cpl.LPES(oi,mi) = sum(LPESc(oi,mi,1:3));
        cpl.gPE(oi,mi) = sum(cpl.gPEc(oi,mi,1:3));
        cpl.gOM(oi,mi) = cpl.gMB(oi,mi) + cpl.gPE(oi,mi);
        
        % display results
        disp(['  wM = ',num2str(wM*1e-9,'%.2f'),' GHz, g0 = ',...
                 num2str(real(cpl.gMB(oi,mi))*1e-3),' + ',...
                 num2str(real(cpl.gPE(oi,mi))*1e-3),' = ',...
                 num2str(real(cpl.gOM(oi,mi))*1e-3),' kHz']);
        
        %% save params at max OM coupling     
        if (sgnCpl~=0 && abs(real(cpl.gOM(oi,mi))) > cpl.gMax)
            
            cpl.gMax = abs(real(cpl.gOM(oi,mi)));
            cpl.gOMmax = real(cpl.gOM(oi,mi));
            cpl.gMBmax = real(cpl.gMB(oi,mi));
            cpl.gPEmax = real(cpl.gPE(oi,mi));
            cpl.oSol = oi; 
            cpl.mSol = mi;
            cpl.optWvl = lambdaAll(oi);
            cpl.mechFreq = wM;
            cpl.Q = ofem.QAll(oi);
            
        end
        
        mIdx = mIdx + 1;
    end
    
    if sgnCpl==0
        cpl.gMax = abs(real(cpl.gOM(oi,mfem.locInd)));
        cpl.gOMmax = real(cpl.gOM(oi,mfem.locInd));
        cpl.gMBmax = real(cpl.gMB(oi,mfem.locInd));
        cpl.gPEmax = real(cpl.gPE(oi,mfem.locInd));
        cpl.oSol = oi; 
        cpl.mSol = mfem.locInd;
        cpl.optWvl = lambdaAll(oi);
        cpl.mechFreq = mfem.freqs(mfem.locInd);
        cpl.Q = ofem.QAll(oi);
    end
    
    oIdx = oIdx + 1;
end

ds.cpl = cpl;

end