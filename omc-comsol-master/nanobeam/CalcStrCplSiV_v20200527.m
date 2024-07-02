% workflow:
% port scripts to new version folder

% loop over orientations
% at each orientation, create new subfield cpl.lSiV_(orientation)
% loop over coordinates (and high symmetry planes)
% interpolate e_IJ (in simulation basis) at coordinates
% assemble strain tensor e = [e_IJ]
% rotate for specified orientations
% populate cpl.lSiV_(ori) with subfields coords, e_ij (in SiV basis)



% reshape SiV strain into 3D array?


% in plotting script,


% readme:
% zSiV = cell containing combinations of [1 1 1],[-1 1 1],[-1 -1 1],[1 -1 1] specifying
% orientation of SiV axis


function [ds,model,dsSiV] = CalcStrCplSiV_v20200527(ds,model,mModes,xpos,ypos,zpos,zSiVAll,runStats,dispOutput)
tCalcStart = tic;

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

% generate coords to extract strain tensor over
[xm,ym,zm] = meshgrid(xpos,ypos,zpos); % array of size length(yCoords) x length(xCoords(xi)) x length(zCoords)
coord = transpose([xm(:),ym(:),zm(:)]);   % array of size 3x(length(xCoords(xi))*length(yCoords)*length(zCoords))

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

%% get rotated strain tensors in SiV basis
if dispOutput
    disp(' ')
    disp('Calculating strain coupling: rotating strain into SiV basis...')
end
% tAllStart = tic;
% tRotStart = tic;
[~,~,dsSiV] = rotateStrainIntoSiV(ds,model,mModes,xpos,ypos,zpos,zSiVAll);
% tRotEnd = toc(tRotStart);
% disp(['Rotation time = ',num2str(tRotEnd,'%.1f'),' s'])

%% Calculate strain coupling
% rotateStrainIntoSiV outputs dsSiV with subfields (zSiV).(sector)
% each with fields:
% eSiVAll - size 3 x 3 x nCoords x max(mModes) (nCoords =
% length(xpos)*length(ypos)*length(zpos))
% e_Egx, e_Egy, LSiV - size nCoords x max(mModes)
% nested for loops: orientation of SiV --> symmetry
% sectors in beam

%% loop over orientations
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
    
    if dispOutput
        disp(' ')
        disp(['  for SiV axis along ',num2str(zSiV)])
    end
    
    % get string representation of orientation
    % for reference to subfield in data struct
    zSiVstr = num2str(zSiV);
    zSiVstr = strrep(zSiVstr,'-','m');
    zSiVstr = strrep(zSiVstr,' ','');
    %% loop over symmetry sectors
    % counters for max strain coupling over sectors
    eSiVFullMax = zeros(3,3,max(mModes));
    LSiVFullMax = zeros(1,max(mModes));
    LSiVFullMaxSec = zeros(3,max(mModes));
    LSiVFullMaxIdx = zeros(1,max(mModes));
    LSiVFullMaxPos = zeros(3,max(mModes));
    
    % cumulative array for strain coupling in apertures
    LSiVFullApers = [];
    %%
    for si = secIdxs
        % get string representation of symmetry sector
        % for reference to subfield in data struct
        secStr = num2str(secSgnAll(:,si)');
        secStr = strrep(secStr,'-1','m');
        secStr = strrep(secStr,'1','p');
        secStr = strrep(secStr,' ','');
        
        %% find max SiV coupling in sector for all modes
        LSiV = dsSiV.(['z',zSiVstr]).(secStr).LSiV; % size = max(mModes)*size(coords,2)
        [LSiVSecMax, LSiVSecMaxIdx] = max(LSiV,[],2);
        LSiVSecMaxPos = coord(:,LSiVSecMaxIdx);
        % apply coordinate sign flips based on sector
        for i = 1:3
            LSiVSecMaxPos(i,:) = secSgnAll(i,si)*LSiVSecMaxPos(i,:);
        end
        
        dsSiV.(['z',zSiVstr]).(secStr).LSiVMax = LSiVSecMax;
        dsSiV.(['z',zSiVstr]).(secStr).LSiVMaxIdx = LSiVSecMaxIdx;
        dsSiV.(['z',zSiVstr]).(secStr).LSiVMaxPos = LSiVSecMaxPos;
        
        cpl.(['SiV_',zSiVstr]).(secStr).LSiVMax = LSiVSecMax;
        cpl.(['SiV_',zSiVstr]).(secStr).LSiVMaxIdx = LSiVSecMaxIdx;
        cpl.(['SiV_',zSiVstr]).(secStr).LSiVMaxPos = LSiVSecMaxPos;
        
        % update counter for max strain coupling over sectors
%         mIdx = 1;
        for mi = mModes
            eSiVSecMax = dsSiV.(['z',zSiVstr]).(secStr).eSiV(:,:,LSiVSecMaxIdx(mi),mi);
            dsSiV.(['z',zSiVstr]).(secStr).eSiVMax(:,:,mi) = eSiVSecMax;
            cpl.(['SiV_',zSiVstr]).(secStr).eSiVMax = eSiVSecMax;
            
            
            if LSiVSecMax(mi) > LSiVFullMax(mi)
                eSiVFullMax(:,:,mi) = eSiVSecMax;
                LSiVFullMax(mi) = LSiVSecMax(mi);
                LSiVFullMaxSec(:,mi) = secSgnAll(:,si);
                LSiVFullMaxIdx(mi) = LSiVSecMaxIdx(mi);
                posTmp = LSiVSecMaxPos(:,mi);
                LSiVFullMaxPos(:,mi) = posTmp;
%                 mIdx = mIdx + 1;
%                 disp(LSiVFullMaxPos)
            end
        end
        
        %% find min/max/mean/std dev of coupling in aperture
        if ~isempty(runStats)
            xIdxs = find(coord(1,:)<=runStats.xmax & coord(1,:)>=runStats.xmin);
            yIdxs = find(coord(2,:)<=runStats.ymax & coord(2,:)>=runStats.ymin);
            zIdxs = find(coord(3,:)<=runStats.zmax & coord(3,:)>=runStats.zmin);
            idxs = intersect(xIdxs, yIdxs);
            idxs = intersect(idxs, zIdxs);
            LSiVRun = LSiV(:,idxs);
            LSiVFullApers = [LSiVFullApers, LSiVRun];
            LSiVAperMax = max(LSiVRun,[],2,'omitnan');
            LSiVAperMin = min(LSiVRun,[],2,'omitnan');
            LSiVAperMean = mean(LSiVRun,2,'omitnan');
            LSiVAperSDev = std(LSiVRun,0,2,'omitnan');
            
%             if dispOutput
%                 disp(['  sec ',secStr,' - mean cpl (+/-1 sd) in aper for'])
%                 for mi = mModes
%                     disp(['  mode ',num2str(mi),', wM=',num2str(mfem.freqs(mi)/1e9,'%.2f'),'GHz: (',...
%                           num2str(mean(LSiVAperMean(mi))*1e-6,'%.2f'),'+/-',...
%                           num2str(mean(LSiVAperSDev(mi))*1e-6,'%.2f'),') MHz'])
%                 end
%             end
            
            dsSiV.(['z',zSiVstr]).(secStr).LSiVAperMax = LSiVAperMax;
            dsSiV.(['z',zSiVstr]).(secStr).LSiVAperMin = LSiVAperMin;
            dsSiV.(['z',zSiVstr]).(secStr).LSiVAperMean = LSiVAperMean;
            dsSiV.(['z',zSiVstr]).(secStr).LSiVAperSDev = LSiVAperSDev;
            
            cpl.(['SiV_',zSiVstr]).(secStr).LSiVAperMax = LSiVAperMax;
            cpl.(['SiV_',zSiVstr]).(secStr).LSiVAperMin = LSiVAperMin;
            cpl.(['SiV_',zSiVstr]).(secStr).LSiVAperMean = LSiVAperMean;
            cpl.(['SiV_',zSiVstr]).(secStr).LSiVAperSDev = LSiVAperSDev;
        end
        
    end
    
    %% find max SiV coupling in full beam for all modes
    dsSiV.(['z',zSiVstr]).full.eSiVMax = eSiVFullMax;
    dsSiV.(['z',zSiVstr]).full.LSiVMax = LSiVFullMax;
    dsSiV.(['z',zSiVstr]).full.LSiVMaxSec = LSiVFullMaxSec;
    dsSiV.(['z',zSiVstr]).full.LSiVMaxIdx = LSiVFullMaxIdx;
    dsSiV.(['z',zSiVstr]).full.LSiVMaxPos = LSiVFullMaxPos;

    cpl.(['SiV_',zSiVstr]).full.eSiVMax = eSiVFullMax;
    cpl.(['SiV_',zSiVstr]).full.LSiVMax = LSiVFullMax;
    cpl.(['SiV_',zSiVstr]).full.LSiVMaxSec = LSiVFullMaxSec;
    cpl.(['SiV_',zSiVstr]).full.LSiVMaxIdx = LSiVFullMaxIdx;
    cpl.(['SiV_',zSiVstr]).full.LSiVMaxPos = LSiVFullMaxPos;
    
    if dispOutput
        disp('  Max coupling in full beam for')
        for mi = mModes
            secsgn = LSiVFullMaxSec(:,mi);
            secStr = num2str(LSiVFullMaxSec(:,mi)');
            secStr = strrep(secStr,'-1','m');
            secStr = strrep(secStr,'1','p');
            secStr = strrep(secStr,' ','');

            disp(['  mode ',num2str(mi),', wM=',num2str(mfem.freqs(mi)/1e9,'%.2f'),'GHz: ',...
                  num2str(LSiVFullMax(mi)*1e-6,'%.2f'),' MHz ',...
                  'in sec ',secStr,' ',...
                  'at pos (',num2str(secsgn(1)*LSiVFullMaxPos(1,mi)*1e9,'%.1f'),',',...
                             num2str(secsgn(2)*LSiVFullMaxPos(2,mi)*1e9,'%.1f'),',',...
                             num2str(secsgn(3)*LSiVFullMaxPos(3,mi)*1e9,'%.1f'),')nm'])
        end
    end
    
    %% find min/max/mean/std dev of coupling in apertures for full beam
    if ~isempty(runStats)
        
        LSiVAperMax = max(LSiVFullApers,[],2,'omitnan');
        LSiVAperMin = min(LSiVFullApers,[],2,'omitnan');
        LSiVAperMean = mean(LSiVFullApers,2,'omitnan');
        LSiVAperSDev = std(LSiVFullApers,0,2,'omitnan');
        
        dsSiV.(['z',zSiVstr]).full.LSiVAperMax = LSiVAperMax;
        dsSiV.(['z',zSiVstr]).full.LSiVAperMin = LSiVAperMin;
        dsSiV.(['z',zSiVstr]).full.LSiVAperMean = LSiVAperMean;
        dsSiV.(['z',zSiVstr]).full.LSiVAperSDev = LSiVAperSDev;

        cpl.(['SiV_',zSiVstr]).full.LSiVAperMax = LSiVAperMax;
        cpl.(['SiV_',zSiVstr]).full.LSiVAperMin = LSiVAperMin;
        cpl.(['SiV_',zSiVstr]).full.LSiVAperMean = LSiVAperMean;
        cpl.(['SiV_',zSiVstr]).full.LSiVAperSDev = LSiVAperSDev;
        
        if dispOutput
            disp(['  mean cpl for full beam in aperture: ',...
                  'x=(',num2str(runStats.xmin*1e9,'%.1f'),',',...
                        num2str(runStats.xmax*1e9,'%.1f'),'), ',...
                  'y=(',num2str(runStats.ymin*1e9,'%.1f'),',',...
                        num2str(runStats.ymax*1e9,'%.1f'),'), ',...
                  'z=(',num2str(runStats.zmin*1e9,'%.1f'),',',...
                        num2str(runStats.zmax*1e9,'%.1f'),') for'])
            for mi = mModes
                disp(['  mode ',num2str(mi),', wM=',num2str(mfem.freqs(mi)/1e9,'%.2f'),'GHz: (',...
                      num2str(mean(LSiVAperMean(mi))*1e-6,'%.2f'),'+/-',...
                      num2str(mean(LSiVAperSDev(mi))*1e-6,'%.2f'),') MHz'])
            end
        end
    end
end
% tAllEnd = toc(tAllStart);
% disp(['Total time = ',num2str(tAllEnd,'%.1f'),' s'])


% ds.dsSiV = dsSiV;
if ~P.calcG
    % find max str cpl across all modes
    if ~isempty(runStats)
        [LSiVMaxAll,LSiVMaxMode] = max(dsSiV.(['z',zSiVstr]).full.LSiVAperMax);
        cpl.mechFreq = mfem.freqs(LSiVMaxMode);
        cpl.lambdaGMax = dsSiV.(['z',zSiVstr]).full.LSiVMax(LSiVMaxMode);
    else
        [LSiVMaxAll,LSiVMaxMode] = max(dsSiV.(['z',zSiVstr]).full.LSiVMax);
        cpl.mechFreq = mfem.freqs(LSiVMaxMode);
        cpl.lambdaGMax = LSiVMaxAll;
    end
%     [~,maxIdx] = max(real(cpl.lambdaGMax));
%     cpl.mechFreq = mfem.freqs(maxIdx);
end
% cpl.(['SiV_',zSiVstr]) = dsSiV.(['z',zSiVstr]);
ds.cpl = cpl;


tCalcEnd = toc(tCalcStart);
disp(['Calc time = ',num2str(tCalcEnd,'%.5f'),' s'])

end