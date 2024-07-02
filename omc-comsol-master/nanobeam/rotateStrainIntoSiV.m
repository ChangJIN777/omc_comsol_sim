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
% zSiV = row vector [1 1 1],[-1 1 1],[-1 -1 1],[1 -1 1] specifying
% orientation of SiV axis


function [ds,model,dsSiV] = rotateStrainIntoSiV(ds,model,mModes,xpos,ypos,zpos,zSiVAll)
% preliminaries
if ~isfield(ds,'cpl')
    cpl = [];
else
    cpl = ds.cpl;
end

mfem = ds.mfem;     % mechanical FEM solution data struct
P = ds.P;           % param data struct
diadom = mfem.dia_domind;   % domain index for diamond in mechanical FEM simulations
cpl.xzpf = mfem.xzpf;       % zero-point displacements

% strain susceptibilities
dg = P.dg;    
fg = P.fg;

% conversion factor from orbital to spin strain susceptibility
ds_do = sqrt(2)*mfem.freqs/46e9;

% crystal orientation of beam
rxtal = P.rxtal; 

% get max disp of each mode for normalizatin
maxDisp = mphmax(model,'abs(solid.disp)','volume','dataset','mdset',...
                 'selection','all','solnum',1:P.mneigs);
cpl.maxDisp = maxDisp;

% generate coords to extract strain tensor over
[xm,ym,zm] = meshgrid(xpos,ypos,zpos); % array of size length(yCoords) x length(xCoords(xi)) x length(zCoords)
coord = transpose([xm(:),ym(:),zm(:)]);   % array of size 3x(length(xCoords(xi))*length(yCoords)*length(zCoords))

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

% initialize strain coupling array
LSiV = zeros(max(mModes),size(coord,2));

%% get rotated strain tensor in SiV basis
% nested for loops: mechanical modes --> orientation of SiV --> symmetry
% sectors in beam

%% loop over modes
for mi = mModes
    % for each mode, get elements of strain tensor 
    % in simulation basis at all coordinates
    % (Voigt ordering used)
    [eXXsim,eYYsim,eZZsim,eYZsim,eXZsim,eXYsim] = mphinterp(model,...
        {'solid.eXX','solid.eYY','solid.eZZ', ...
         'solid.eYZ','solid.eXZ','solid.eXY'},'coord',coord,...
         'dataset','mdset','edim',3,'selection',diadom,'solnum',mi);
    
    %% loop over specified orientations of SiV
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

        % get string representation of orientation
        % for reference to subfield in data struct
        zSiVstr = num2str(zSiV);
        zSiVstr = strrep(zSiVstr,'-','m');
        zSiVstr = strrep(zSiVstr,' ','');
        % disp(['calculating strain coupling for ',zSiVstr,' oriented SiV ...'])
        
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

        % get crystal rotation tensor from simulation-->geom basis
        RZ = [cos(rxtal) -sin(rxtal) 0;...
              sin(rxtal)  cos(rxtal) 0;...
              0           0          1];

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

        % save parameters to new data struct
        dsSiV.(['z',zSiVstr]).xSiV = xSiV;
        dsSiV.(['z',zSiVstr]).ySiV = ySiV;
        dsSiV.(['z',zSiVstr]).zSiV = zSiV;
        dsSiV.(['z',zSiVstr]).coord = coord;
        dsSiV.(['z',zSiVstr]).RZ = RZ;
        dsSiV.(['z',zSiVstr]).RSiV = RSiV;
        dsSiV.(['z',zSiVstr]).RTot = RTot;
        
        %% loop over sectors where symmetries apply
        % multiply components by factors depending on sector and symmetries
        % pj = parity about j-axis
        for si = secIdxs
            % get string representation of symmetry sector
            % for reference to subfield in data struct
            secStr = num2str(secSgnAll(:,si)');
            secStr = strrep(secStr,'-1','m');
            secStr = strrep(secStr,'1','p');
            secStr = strrep(secStr,' ','');

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
            eXX = eSgn(1)*eXXsim;
            eYY = eSgn(2)*eYYsim;
            eZZ = eSgn(3)*eZZsim;
            eYZ = eSgn(4)*eYZsim;
            eXZ = eSgn(5)*eXZsim;
            eXY = eSgn(6)*eXYsim;

            %% loop over coordinates
            eSiVAll = zeros(3,3,size(coord,2));
            for ci = 1:size(coord,2)
                % assemble strain tensor in simulation basis
                eSim = [eXX(ci), eXY(ci), eXZ(ci);
                        eXY(ci), eYY(ci), eYZ(ci);
                        eXZ(ci), eYZ(ci), eZZ(ci)];
                    
                % rotate strain tensor to SiV basis
                eSiV = RTot*eSim*transpose(RTot);
                eSiVAll(:,:,ci) = eSiV(:,:);
            end
            dsSiV.(['z',zSiVstr]).(secStr).eSiV(:,:,:,mi) = eSiVAll(:,:,:);
            
            %% Calculate strain coupling for each orientation and sector
            % indexing of eSiVAll: i-component, j-component, coordinate, mode
%             eSiVAll = dsSiV.(['z',zSiVstr]).(secStr).eSiV;
            
            e_Egx = dg*(eSiVAll(1,1,:)-eSiVAll(2,2,:)) + fg*eSiVAll(2,3,:);
            e_Egy = -2*dg*eSiVAll(1,2,:) + fg*eSiVAll(1,3,:);
            
            % here, eSiVAll has size (3,3,nCoords)
            % and e_Egx and e_Egy have size (1,1,nCoords)
            % --> remove singleton dimensions
            % --> now have size (nCoords)
            e_Egx = squeeze(e_Egx);
            e_Egy = squeeze(e_Egy);
            
            LSiVm = real(ds_do(mi).*sqrt(abs(e_Egx).^2+abs(e_Egy).^2)...
                        .*cpl.xzpf(mi)./maxDisp(mi));
            LSiV(mi,:) = LSiVm(:);
            
            % save data
            dsSiV.(['z',zSiVstr]).(secStr).e_Egx = e_Egx;
            dsSiV.(['z',zSiVstr]).(secStr).e_Egy = e_Egy;
            dsSiV.(['z',zSiVstr]).(secStr).LSiV = LSiV;

        end
    end
    
end

%% Calculate strain coupling
% loop over orientations --> sectors --> modes for efficiency
% as MATLAB handles operations on large arrays better than looping over
% smaller arrays

% for oi = 1:length(zSiVAll)
%     zSiV = zSiVAll{oi};
%     disp(['for SiV axis along ',num2str(zSiV)])
%     
%     % get string representation of orientation
%     % for reference to subfield in data struct
%     zSiVstr = num2str(zSiV);
%     zSiVstr = strrep(zSiVstr,'-','m');
%     zSiVstr = strrep(zSiVstr,' ','');
%     for si = secIdxs
%         % get string representation of symmetry sector
%         % for reference to subfield in data struct
%         secStr = num2str(secSgnAll(:,si)');
%         secStr = strrep(secStr,'-1','m');
%         secStr = strrep(secStr,'1','p');
%         secStr = strrep(secStr,' ','');
%         
%         %% Calculate strain coupling for each orientation and sector
%         % indexing of eSiVAll: i-component, j-component, coordinate, mode
%         eSiVAll = dsSiV.(['z',zSiVstr]).(secStr).eSiV;
%         
%         e_Egx = dg*(eSiVAll(1,1,:,:)-eSiVAll(2,2,:,:)) + fg*eSiVAll(1,3,:,:);
%         e_Egy = -2*dg*eSiVAll(1,2,:,:) + fg*eSiVAll(2,3,:,:);
%         
%         % here, eSiVAll has size (3,3,nCoords,max(mModes))
%         % and e_Egx and e_Egy have size (1,1,nCoords,max(mModes))
%         % --> remove singleton dimensions 
%         % --> now have size (nCoords,max(mModes))
%         e_Egx = squeeze(e_Egx);
%         e_Egy = squeeze(e_Egy);
%         
%         LSiV = zeros(size(e_Egx));
%         for mi = mModes
%             LSiVm = real(sqrt(abs(e_Egx(:,mi)).^2+abs(e_Egy(:,mi)).^2)...
%                     .*cpl.xzpf(mi)./maxDisp(mi));
%             LSiV(:,mi) = LSiVm;
%         end
%         
%         % save data
%         dsSiV.(['z',zSiVstr]).(secStr).e_Egx = e_Egx;
%         dsSiV.(['z',zSiVstr]).(secStr).e_Egy = e_Egy;
%         dsSiV.(['z',zSiVstr]).(secStr).LSiV = LSiV;
%     end
% end

end