% readme:
% zSiV = cell containing combinations of [1 1 1],[-1 1 1],[-1 -1 1],[1 -1 1] specifying
% orientation of SiV axis

% add varargout for saving of strain tensor component - helpful for
% plotting YZ profiles

function [ds,model,dsSiV] = CalcStrCplSiV(ds,model,mModes,xpos,ypos,zpos,zSiVAll,runStats,dispOutput,varargin)
% timer
if dispOutput
    tCalcStart = tic;
end

% preliminaries
if ~isfield(ds,'cpl')
    cpl = [];
else
    cpl = ds.cpl;
end

% extract results
mfem = ds.mfem;     % mechanical FEM solution data struct
P = ds.P;           % param data struct
diadom = mfem.dia_domind;   % domain index for diamond in mechanical FEM simulations
cpl.xzpf = mfem.xzpf;       % zero-point displacements

% strain susceptibilities
dg = P.dg;    
fg = P.fg;

% conversion factor from orbital to spin strain susceptibility
dspin_dorb = sqrt(2)*mfem.freqs/46e9;

% crystal orientation of beam
rxtal = P.rxtal; 

% generate coords to extract strain tensor over
[xm,ym,zm] = meshgrid(xpos,ypos,zpos); % array of size length(yCoords) x length(xCoords(xi)) x length(zCoords)
coord = transpose([xm(:),ym(:),zm(:)]);   % array of size 3x(length(xCoords(xi))*length(yCoords)*length(zCoords))
ncoord = size(coord,2);

% get max disp of each mode for normalization
maxDisp = mphmax(model,'abs(solid.disp)','volume','dataset','mdset',...
                 'selection','all','solnum',1:P.mneigs);
cpl.maxDisp = maxDisp;

% create 3x8 array, where each column contains x/y/z symmetry of sector
% [x+ x+ x- x- x+ x+ x- x-;...
%  y+ y- y+ y- y+ y- y+ y-;...
%  z+ z+ z+ z+ z- z- z- z-]

% generate geometry signs
[XG,YG,ZG] = meshgrid([1 -1],[1 -1],[1 -1]);
secSgnAll = transpose([XG(:),YG(:),ZG(:)]);

% generate symmetries for all sectors
[XV,YV,ZV] = meshgrid([1 P.mevenx],[1 P.meveny],[1 P.mevenz]);
secSymAll = transpose([XV(:),YV(:),ZV(:)]); 
secIdxs = find(prod(secSymAll)~=0); % find sectors where symmetries apply

% if calc mean coupling in aperture, find coordinates to sample coupling
% and run stats over
if ~isempty(runStats)
    xIdxs = find(coord(1,:)<=runStats.xmax & coord(1,:)>=runStats.xmin);
    yIdxs = find(coord(2,:)<=runStats.ymax & coord(2,:)>=runStats.ymin);
    zIdxs = find(coord(3,:)<=runStats.zmax & coord(3,:)>=runStats.zmin);
    aperIdxs = intersect(xIdxs, yIdxs);
    aperIdxs = intersect(aperIdxs, zIdxs);
end

% get string representation of symmetry sector
secStr = cell(1,8);
for si = secIdxs
    secStr{si} = num2str(secSgnAll(:,si)');
    secStr{si} = strrep(secStr{si},'-1','m');
    secStr{si} = strrep(secStr{si},'1','p');
    secStr{si} = strrep(secStr{si},' ','');
end

% get string representation of orientation
zSiVstr = cell(1,length(zSiVAll));
for oi = 1:length(zSiVAll)
    zSiV = zSiVAll{oi};

    % for non-standard orientations, convert to row vector and invert signs
    [nr,nc] = size(zSiV);
    if nr==3 && nc==1
        zSiV = zSiV';
    end
    if (prod(zSiV==[-1 -1 -1]) || prod(zSiV==[ 1 -1 -1]) || ...
        prod(zSiV==[ 1  1 -1]) || prod(zSiV==[-1  1 -1]))
        zSiV = -1*zSiV;
    end
    zSiVAll{oi} = zSiV;
    
    zSiVstr{oi} = num2str(zSiV);
    zSiVstr{oi} = strrep(zSiVstr{oi},'-','m');
    zSiVstr{oi} = strrep(zSiVstr{oi},' ','');
end

% 1st varargin argument - coordinates
% specify as 1x3 array of coords
% or as 3 x 1 x (no. of mech modes) x (no. of SiV orientation) array
if ~isempty(varargin)
    posT = varargin{1};
    posTlen = length(posT);
    
end

%% get rotated strain tensors in SiV basis
if dispOutput
    disp(' ')
    disp('Calculating spin-phonon coupling ...')
end

%% loop over mechanical modes

% get crystal rotation tensor from simulation-->geom basis
RZ = [cos(rxtal) -sin(rxtal) 0;...
      sin(rxtal)  cos(rxtal) 0;...
      0           0          1];

for mIdx = 1:length(mModes) % mi = mModes
    mi = mModes(mIdx);
    if dispOutput
        disp(' ')
        disp(['  mode ',num2str(mi),', wM=',num2str(mfem.freqs(mi)/1e9,'%.2f'),'GHz: '])
    end
    
    % for each mode, get strain tensor elements in simulation basis 
    % at all coordinates (Voigt ordering used)
    % eIJ's are row vectors of length = no. of coords (size(coord,2))
%     tInterpStart = tic;
    [eXXsim,eYYsim,eZZsim,eYZsim,eXZsim,eXYsim] = mphinterp(model,...
        {'solid.eXX','solid.eYY','solid.eZZ', ...
         'solid.eYZ','solid.eXZ','solid.eXY'},'coord',coord,...
         'dataset','mdset','edim',3,'selection',diadom,'solnum',mi);
     
    % convert to real array of size 6 x no. of coords to save memory
    % row ordering follows Voigt notation
    esimAll = real([eXXsim;eYYsim;eZZsim;eYZsim;eXZsim;eXYsim]);
    clear eXXsim eYYsim eZZsim eYZsim eXZsim eXYsim
%     tInterpEnd = toc(tInterpStart);
%     disp(['Interp time = ',num2str(tInterpEnd,'%.5f'),' s']) 
    
    %% loop over specified orientations of SiV
    for oi = 1:length(zSiVAll)
        zSiV = zSiVAll{oi};
        
        %% get rotation matrices for each SiV orientation
        % Goal is to transform strain tensor from simulation basis to SiV basis
        % From simulations: [xyz]_geom corresponds to x'y'z' = R1(rxtal)
        % and rotation matrix R2 = [xSiV ySiV zSiV] corresponds to 
        % rotation from [xyz]_SiV = I to [xyz]_geom
        % Therefore, perform this sequence to transform strain tensor:
        % 1. simulation --> geom basis: eSim --> eGeom = (RZ)*eSim*(RZ^-1)
        % 2. geom --> SiV basis: eGeom --> eSiV = (RSiV)*eGeom*(RSiV^-1)
        % i.e. eSiV = RSiV*RZ*eSim*RZ^-1*RSiV^-1
        % where RZ = R1^-1, RSiV = R2^-1
        % for orthonormal basis M, M^-1 = M^T
        
        % get normalized SiV basis vectors
        if prod(zSiV==[1 1 1]) 
            ySiV = [1 -1 0];
        elseif prod(zSiV==[-1 1 1])
            ySiV = [1 1 0];
        elseif prod(zSiV==[-1 -1 1])
            ySiV = [-1 1 0];
        elseif prod(zSiV==[1 -1 1])
            ySiV = [-1 -1 0];
        end
        xSiV = cross(ySiV,zSiV);

        normalize = @(v) v / sqrt(sum(v .* v)); % normalize vectors
        zSiV = normalize(zSiV);
        xSiV = normalize(xSiV);
        ySiV = normalize(ySiV);
        
        % rotation from SiV to cartesian basis
        % columns correspond to normalized basis vectors
        % Z aligned along SiV axis; Y aligned along in-plane SiV bond
        USiV = [xSiV' ySiV' zSiV'];
        RSiV = transpose(USiV);
        
        % rotation from simulation-->geom-->SiV basis
        RTot = RSiV*RZ;

        %% loop over sectors where symmetries apply
        % multiply components by factors depending on sector and symmetries
        % pj = parity about j-axis
        for si = secIdxs
            % default sector: x>0, y>0, z>0
            eSgn = ones(1,6);

            % x --> -x: multiply -px for shear terms with x, +px otherwise
            if secSgnAll(1,si) == -1
                eSgn = eSgn.*P.mevenx.*[ 1  1  1  1 -1 -1];
            end

            % y --> -y: multiply -py for shear terms with y, +py otherwise
            if secSgnAll(2,si) == -1
                eSgn = eSgn.*P.meveny.*[ 1  1  1 -1  1 -1];
            end

            % z --> -z: multiply -pz for shear terms with z, +pz otherwise
            if secSgnAll(3,si) == -1
                eSgn = eSgn.*P.mevenz.*[ 1  1  1 -1 -1  1];
            end
            
            % apply symmetry operations on strain tensor components
%             tSecStart = tic;
            eSgn = eSgn';
            esimSec = eSgn.*esimAll;
%             eXX = eSgn(1)*eXXsim;
%             eYY = eSgn(2)*eYYsim;
%             eZZ = eSgn(3)*eZZsim;
%             eYZ = eSgn(4)*eYZsim;
%             eXZ = eSgn(5)*eXZsim;
%             eXY = eSgn(6)*eXYsim;
%             esimAll = real([eXX;eYY;eZZ;eYZ;eXZ;eXY]);
%             tSecEnd = toc(tSecStart);
%             disp(['Sec time = ',num2str(tSecEnd,'%.5f'),' s'])
            
            %% loop over coordinates
            eSiVAll = zeros(6,ncoord);
%             tRotStart = tic;
            for ci = 1:ncoord
                % assemble strain tensor in simulation basis
                eSim = [esimSec(1,ci), esimSec(6,ci), esimSec(5,ci);
                        esimSec(6,ci), esimSec(2,ci), esimSec(4,ci);
                        esimSec(5,ci), esimSec(4,ci), esimSec(3,ci)];
                    
                % rotate strain tensor to SiV basis
                eSiVTensor = RTot*eSim*transpose(RTot);
                eSiVVoigt = eSiVTensor([1,5,9,8,7,4]);
                eSiVAll(1:6,ci) = eSiVVoigt(1:6)';
            end % of loop over coordinates
            clear esimSec
%             tRotEnd = toc(tRotStart);
%             disp(['Rot time = ',num2str(tRotEnd,'%.5f'),' s'])
            
            %% Calculate strain coupling for each orientation and sector
%             tCalcStart = tic;
            e_Egx = dg*(eSiVAll(1,:)-eSiVAll(2,:)) + fg*eSiVAll(4,:);
            e_Egy = -2*dg*eSiVAll(6,:) + fg*eSiVAll(5,:);
            LSiV = real(dspin_dorb(mi).*sqrt(abs(e_Egx).^2+abs(e_Egy).^2)...
                        .*cpl.xzpf(mi)./maxDisp(mi));
            
            % indexing: compts, coords, mech mode, orientation, sector
            dsSiV.e_Egx(1,1:ncoord,mIdx,oi,si) = e_Egx(1:ncoord);
            dsSiV.e_Egy(1,1:ncoord,mIdx,oi,si) = e_Egy(1:ncoord);
            dsSiV.LSiV(1,1:ncoord,mIdx,oi,si) = LSiV(1:ncoord);
            
            % find max coupling in sector and extract coord and strain
            % tensor
            [LSiVMax,LSiVMaxIdx] = max(LSiV);
            dsSiV.LSiVMax(1,1,mIdx,oi,si) = LSiVMax;
            dsSiV.LSiVMaxPos(1:3,1,mIdx,oi,si) = secSgnAll(:,si).*coord(:,LSiVMaxIdx);
            dsSiV.eSiVMax(1:6,1,mIdx,oi,si) = eSiVAll(1:6,LSiVMaxIdx);
            
            clear e_Egx e_Egy
            
            % find min/max/mean/sdev of coupling in aperture in each sector
            if ~isempty(runStats)
                LSiVRun = LSiV(aperIdxs);
                LSiVAperMax = max(LSiVRun,[],'omitnan');
                LSiVAperMin = min(LSiVRun,[],'omitnan');
                LSiVAperMean = mean(LSiVRun,'omitnan');
                LSiVAperSDev = std(LSiVRun,0,'omitnan');
                
                dsSiV.LSiVAperMax(1,1,mIdx,oi,si) = LSiVAperMax;
                dsSiV.LSiVAperMin(1,1,mIdx,oi,si) = LSiVAperMin;
                dsSiV.LSiVAperMean(1,1,mIdx,oi,si) = LSiVAperMean;
                dsSiV.LSiVAperSDev(1,1,mIdx,oi,si) = LSiVAperSDev;
                
            end
            
            clear LSiV
            
            % save strain tensor at specified location
            if ~isempty(varargin)
                if posTlen == 3
                    posSec = posT;
                elseif posTlen ~=3 && length(size(posT)) == 4
                    posSec = posT(1:3,1,mIdx,oi);
                end
%                 disp(posSec)
                
                % invert signs if symmetry plane
%                 posSec(1) = posSec(1)*(sign(posSec(1)))^(abs(P.mevenx));
%                 posSec(2) = posSec(2)*(sign(posSec(2)))^(abs(P.meveny));
%                 posSec(3) = posSec(3)*(sign(posSec(3)))^(abs(P.mevenz));
%                 disp(posSec)
                
                cxIdxs = find(abs(posSec(1)-coord(1,:))<1e-9);
                cyIdxs = find(abs(posSec(2)-coord(2,:))<1e-9);
                czIdxs = find(abs(posSec(3)-coord(3,:))<1e-9);
                cIdx = intersect(cxIdxs,cyIdxs);
                cIdx = intersect(cIdx,czIdxs);
%                 disp(coord(1,1))
%                 disp(['cx: ',num2str(length(cxIdxs))])
%                 disp(['cy: ',num2str(length(cyIdxs))])
%                 disp(['cz: ',num2str(length(czIdxs))])
                
                
                dsSiV.eSiVpos(1:6,1,mIdx,oi,si) = eSiVAll(1:6,cIdx(1)); 
                    % index 1 for cIdx is for cases where cIdx has length 
                    % greater than 1, i.e. it finds two coords close to 
                    % specified location
                
            end
%             tCalcEnd = toc(tCalcStart);
%             disp(['Calc time = ',num2str(tCalcEnd,'%.5f'),' s'])
            
        end % of loop over symmetry sectors
        
        %% find max coupling in full beam, across all sectors
        [maxFull,maxSec] = max(dsSiV.LSiVMax(1,1,mIdx,oi,:));
        dsSiV.LSiVMaxSec(1,1,mIdx,oi) = maxSec;
%         dsSiV.LSiVMaxFull(1,1,mIdx,oi) = maxFull;
%         dsSiV.LSiVMaxFSec(1,1,mIdx,oi) = maxSec;
%         dsSiV.eSiVMaxFSec(1:6,1,mIdx,oi) = dsSiV.eSiVMax(1:6,1,mIdx,oi,maxSec);
%         dsSiV.LSiVMaxFPos(1:3,1,mIdx,oi) = dsSiV.LSiVMaxPos(1:3,1,mIdx,oi,maxSec);
        
        if dispOutput
            disp(['  SiV // ',zSiVstr{oi},': ',...
                  'max g_sp=',num2str(maxFull/1e6,'%.2f'),'MHz ',...
                  'at (',num2str(dsSiV.LSiVMaxPos(1,1,mIdx,oi,maxSec)*1e9,'%.1f'),',',...
                         num2str(dsSiV.LSiVMaxPos(2,1,mIdx,oi,maxSec)*1e9,'%.1f'),',',...
                         num2str(dsSiV.LSiVMaxPos(3,1,mIdx,oi,maxSec)*1e9,'%.1f'),')nm'])
        end
        
        %% find sector with highest mean coupling in aperture
        if ~isempty(runStats)
            [aperMeanFull,aperMeanSec] = max(dsSiV.LSiVAperMean(1,1,mIdx,oi,:));
            dsSiV.LSiVAperMeanSec(1,1,mIdx,oi) = aperMeanSec;
%             dsSiV.LSiVAperMeanFull(1,1,mIdx,oi) = aperMeanFull;
%             dsSiV.LSiVAperMeanSec(1,1,mIdx,oi) = aperMeanSec;
%             dsSiV.LSiVAperMaxFull(1,1,mIdx,oi) = dsSiV.LSiVAperMax(1,1,mIdx,oi,aperMeanSec);
%             dsSiV.LSiVAperMinFull(1,1,mIdx,oi) = dsSiV.LSiVAperMin(1,1,mIdx,oi,aperMeanSec);
%             dsSiV.LSiVAperSDevFull(1,1,mIdx,oi) = dsSiV.LSiVAperSDev(1,1,mIdx,oi,aperMeanSec);
            if dispOutput
                disp(['        mean g_sp=',num2str(aperMeanFull/1e6,'%.2f'),'MHz ',...
                      'in sector ',secStr{aperMeanSec},' ',...
                      'in aperture ',...
                      'x=(',num2str(runStats.xmin*1e9,'%.1f'),',',...
                            num2str(runStats.xmax*1e9,'%.1f'),'), ',...
                      'y=(',num2str(runStats.ymin*1e9,'%.1f'),',',...
                            num2str(runStats.ymax*1e9,'%.1f'),'), ',...
                      'z=(',num2str(runStats.zmin*1e9,'%.1f'),',',...
                            num2str(runStats.zmax*1e9,'%.1f'),')'])
            end
        end
    end % of loop over SiV orientations
    clear esimAll eSiVAll
    
    % find max for each mode
    LSiVMode = dsSiV.LSiVMax(1,1,mIdx,:,:);
    dsSiV.LSiVMaxMode(mi) = max(LSiVMode(:));
end % of loop over mModes

%% find max coupling in full beam, across all sectors
% sector is 5th dimension - take max over 5th dimension
% LSiVMaxFull, LSiVMaxFullSec = arrays of size 1 x 1 x length(mModes) x
% length(zSiVAll)
% [dsSiV.LSiVMaxFull,dsSiV.LSiVMaxFullSec] = max(dsSiV.LSiVMax,[],5,'omitnan');
% dsSiV.eSiVMaxFull = dsSiV.eSiVMax(1:6,1,:,:,si)

%% 
% find global max
[globalMax,globalInd] = max(dsSiV.LSiVMax(:));
[~,~,mmodeIdx,oriIdx,secIdx] = ind2sub(size(dsSiV.LSiVMax),globalInd);
dsSiV.LSiVMaxGlobal = globalMax;
dsSiV.LSiVMaxGlobalMode = mModes(mmodeIdx);
dsSiV.LSiVMaxGlobalOri = zSiVAll{oriIdx};
dsSiV.LSiVMaxGlobalSec = secSgnAll(:,secIdx);

% copy fields over to cpl
flds = {'LSiVMax','LSiVMaxPos','eSiVMax','LSiVMaxSec',...
        'LSiVMaxMode','LSiVMaxGlobal','LSiVMaxGlobalMode',...
        'LSiVMaxGlobalOri','LSiVMaxGlobalSec'};
if ~isempty(runStats)
    flds = [flds,'LSiVAperMax','LSiVAperMin',...
        'LSiVAperMean','LSiVAperSDev','LSiVAperMeanSec'];
end
for fi = 1:length(flds)
    cpl.SiV.(flds{fi}) = dsSiV.(flds{fi});
end

% extract max cpl and put into ds
if ~P.calcG
    mfreqIdx = mModes(mmodeIdx);
    cpl.mechFreq = mfem.freqs(mfreqIdx);
    cpl.lambdaGMax = globalMax;
end
ds.cpl = cpl;

if dispOutput
    tCalcEnd = toc(tCalcStart);
    disp(['  Calc time = ',num2str(tCalcEnd,'%.5f'),' s'])
end

% disp('done')


end

























