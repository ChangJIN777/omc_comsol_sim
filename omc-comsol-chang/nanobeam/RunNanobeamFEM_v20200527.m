% function to run optical and/or mechanical FEM simulations for nanobeam
% input: data struct (P) with geometry and simulation parameters
% output: data struct (ds) with simulation results, COMSOL model object
% (model); saves model file (.mph) and plots (.png, .fig)

function [ds,model] = RunNanobeamFEM_v20200527(P,datLoc)
tStart = tic;
%% Init
%append \ to end of datLoc if not present
if ~strcmp(datLoc(end),'\')
    datLoc = [datLoc,'\'];
end
% create directory to save files
if ~exist(datLoc,'dir')
    mkdir(datLoc)
end

% import COMSOL class
import com.comsol.model.*
import com.comsol.model.util.*

ModelUtil.showProgress(true);

ModelUtil.clear();
clear model

ds = [];
model = [];

% create COMSOL model named 'model' from which COMSOL methods can be called, 
% e.g. model.save
model = ModelUtil.create('model');

if isfield(P,'scriptLoc')
    addpath(P.scriptLoc); % directory containing scripts used together with FDTD simulations
end

% create base filename for saving of files and plots
if ~isfield(P,'fileBase')
    P = CreateFileBase(P);
    disp(P.fileBase)
end

%% main script
try
    %% Check if simulations have been run
    matFile = dir([datLoc,P.fileBase,'*.mat']);
    mphFile = dir([datLoc,P.fileBase,'*.mph']);
    
    if isempty(matFile) || isempty(mphFile)
        %% Create geometry: generate array of geometry dimensions for all unit cells
        if isfield(P,'celltype') && strcmp(P.celltype,'blockTet')
            P = CreateNanobeamBlockTetGeom(P);              % block-tether geometry
        else
            % create beam-hole geometry by default
            if isfield(P,'asymCav') && P.asymCav && isfield(P,'PL')
                P = CreateNanobeamGeom_asym(P,P.PL);        % creates asymmetric cavity
            else
                P = CreateNanobeamGeom(P);                  % creates symmetric cavity
            end

            % shave geometry - to simulate effect of O2 anneal
            % shave > 0 corresponds to removal of material --> smaller beam
            % width, larger holes
            if isfield(P,'shave')
                % shave holes and width
                P.w = P.w - P.shave;
                P.hx = P.hx + P.shave;
                P.hy = P.hy + P.shave;

                hxG = P.geom(:,1);
                hyG = P.geom(:,2);
                hxGH = P.geomHalf(:,1);
                hyGH = P.geomHalf(:,2);

                hxG = hxG + P.shave;
                hyG = hyG + P.shave;
                hxGH = hxGH + P.shave;
                hyGH = hyGH + P.shave;

                P.geom(:,1) = hxG;
                P.geom(:,2) = hyG;
                P.geomHalf(:,1) = hxGH;
                P.geomHalf(:,2) = hyGH;
            end
            PlotDefectCells(P);                             % plots beam geometry
            saveas(gcf,[datLoc,P.fileBase,'_geom.png']);
        end

        %% Load material params
        P = LoadMaterialParams(P);

        %% Run simulations
        [model,P] = BuildNanobeamFEM(model,P);              % generates nanobeam in COMSOL

        % optional - plot geometry
        if P.plotgeom
            figure; 
            mphgeom(model)
    %         title(P.fileBase)
            saveas(gcf,[datLoc,P.fileBase,'_mphgeom.png']);
        end

        [model,ds] = SetupNanobeamFEM(model,P);         % set up nanobeam FEM simulations
        [model,ds] = SolveNanobeamFEM(model,ds);        % solve and postprocess FEM simulations
        
    elseif ~isempty(matFile) && ~isempty(mphFile)
        disp('loading existing simulation results...')
        load([datLoc,matFile(1).name])
        model = mphload([datLoc,mphFile(1).name]);
        ds.P.calcS = P.calcS;
        ds.P.plotStrCpl = P.plotStrCpl;
        ds.P.plotMech = P.plotMech;
        ds.P.zSiV = P.zSiV;
        ds.P.xSlc = P.xSlc;
        ds.P.ySlc = P.ySlc;
        ds.P.zSlc = P.zSlc;
        ds.P.LStats = P.LStats;
        P = ds.P;
    end
    
    
    %% Calculate coupling(s) - optomechanical, strain
    % optomechanical coupling
    if P.calcG
        [ds,model] = CalcGOM(ds,model,1,ds.mfem.locInd);
    end
    
    %% select modes
    mModes = ds.mfem.locInd;
    bndfreq = find(abs(ds.mfem.freqs-ds.P.freq)<2e9);
    mModes=intersect(mModes,bndfreq);
%     disp(mModes)
    if P.calcG && isfield(P,'g0min')
        highGmodes = find(abs(ds.cpl.gOM)>ds.P.g0min);
        mModes = intersect(mModes,highGmodes);
    end
    %%
    % (orbital) strain coupling
    % find in (half) unit cell at center of beam, and top z-half of beam
    if P.calcS
        % define coordinates to interpolate strain over
        interpResX = 10e-9;
        interpResY = 10e-9;
        interpResZ = 10e-9;
        
%         yposHalf = linspace(0,(max(P.w_hole))/2,ceil((max(P.w_hole))/interpResY)); % span half of width of beam
%         zposHalf = linspace(0,P.th/2,ceil((P.th/2)/interpResZ)); % span half of thickness of beam
        yposHalf = [0:interpResY:(max(P.w_hole))/2,(max(P.w_hole))/2];
        zposHalf = [0:interpResZ:P.th/2,P.th/2];
        
        % for beam with holes, interpolate between center holes
        if ~isfield(P,'celltype')
            ctrHoleRX = ds.P.xc;% + ds.P.holeatctr*ds.P.hx_hole(1)/2;   % central hole right-x-coordinate
            a1 = ds.P.a_hole(1);
            a1(isnan(a1)) = 0;
            adjHoleLX = ds.P.xc ...
                        + (1-ds.P.holeatctr)*a1/2 ...
                        + (ds.P.holeatctr)*ds.P.a_hole(2) ...
                        + ds.P.hx_hole(1+ds.P.holeatctr)/2; % adjacent hole left-x-coordinate
%             xpos = ds.P.xc + [ctrHoleRX:interpResX:adjHoleLX,adjHoleLX];
            xpos = ds.P.xc + [0:interpResX:ds.P.a,ds.P.a];
        else
            xpos = ds.P.xc + linspace(0,P.a,ceil(P.a/interpResX));
        end

        if abs(P.meveny)
            ypos = yposHalf;
        else
            ypos = [-1*fliplr(yposHalf(2:end)) yposHalf];
        end
        
%         zpos = zposHalf;

        if abs(P.mevenz)
            zpos = zposHalf;
        else
            zpos = [-1*fliplr(zposHalf(2:end)) zposHalf];
        end
        
        % if specified, run stats for strain coupling in custom-defined
        % aperture
        if isfield(P,'LStats')
            % define extents of aperture for hole-based geometry
            if ~isfield(P,'celltype') && (isempty(P.LStats.xmin) || isempty(P.LStats.xmax))
%                 ctrHoleRX = P.xc + ds.P.hx_hole(1)/2;   % central hole right-x-coordinate
%                 adjHoleLX = P.xc + ds.P.a_hole(2) - ds.P.hx_hole(2)/2; % adjacent hole left-x-coordinate
%                 
                ctrHoleRX = ds.P.xc + ds.P.holeatctr*ds.P.hx_hole(1)/2;   % central hole right-x-coordinate
                a1 = ds.P.a_hole(1);
                a1(isnan(a1)) = 0;
                adjHoleLX = ds.P.xc ...
                            + (1-ds.P.holeatctr)*a1/2 ...
                            + (ds.P.holeatctr)*ds.P.a_hole(2) ...
                            - ds.P.hx_hole(1+ds.P.holeatctr)/2; % adjacent hole left-x-coordinate

                P.LStats.xmin = ctrHoleRX;
                P.LStats.xmax = adjHoleLX;
            end
            LStats = P.LStats;
        else
            LStats = [];
        end
        
        [ds,model,~] = CalcStrCplSiV_v20200527(ds,model,mModes,xpos,ypos,zpos,P.zSiV,LStats,1);
        
    end
    %
    %% if only calculating strain coupling: save plot for max strain coupling
    if ~P.calcG && P.calcS && P.plotStrCpl
        maxLIdxs = find(ds.mfem.freqs==ds.cpl.mechFreq);
        mModes = intersect(mModes,maxLIdxs);
%         disp(mModes)
    end
    %% Save file
    % filenames with coupling results
    fileName = P.fileBase;
    if P.calcG
        CfBase = ['wl_',num2str(ds.cpl.optWvl*1e9,'%.0f'),'nm_',...
                  'wM_',num2str(ds.cpl.mechFreq*1e-9,'%.2f'),'GHz_',...
                  'gO_',num2str(ds.cpl.gMax*1e-3,'%.0f'),'kHz'];
        fileName = [fileName,'_',CfBase];
    end

    if P.calcS
        LSiVMax = ds.cpl.SiV_111.full.LSiVMax;
        if P.calcG
            CfBase = ['lS_',num2str(real(LSiVMax(ds.cpl.mSol)*1e-6),'%.2f'),'MHz'];
        else
            CfBase = ['wM_',num2str(ds.cpl.mechFreq*1e-9,'%.2f'),'GHz_',...
                      'lS_',num2str(real(max(LSiVMax*1e-6)),'%.2f'),'MHz'];
        end
        fileName = [fileName,'_',CfBase];
    end

    mphsave(model,[datLoc,fileName,'.mph']);
    save([datLoc,fileName,'.mat'],'ds');
    
    %% Plot - electric field, displacement, strain, strain coupling
    % plot (y-component of) electric field
    if P.plotOpt
        PlotEy(model,ds,1,datLoc,1);
    end
    
    % plot mechanical displacement and strain profiles
    if P.plotMech
        PlotDispStr(model,ds,mModes,datLoc,1);
    end
    %%
    if P.plotStrCpl
        % pSlc = data struct containing fields z111, zm111, zm1m11, and/or
        % z1m11, each with subfield pos (array of size 3 x length(modes))
        
        for i = 1:length(P.zSiV)
            zSiVstr = num2str(P.zSiV{i});
            zSiVstr = strrep(zSiVstr,'-','m');
            zSiVstr = strrep(zSiVstr,' ','');
            pSlcSiV = ds.cpl.(['SiV_',zSiVstr]).full.LSiVMaxPos(:,mModes);
            pSlc.(['z',zSiVstr]).pos = abs(pSlcSiV);
        end
        
        PlotStrCplSiVXY_v20200527(ds,model,mModes,pSlc,P.zSiV,datLoc)
        PlotStrCplSiVYZ_v20200527(ds,model,mModes,pSlc,P.zSiV,datLoc)
        
        % plot at custom coordinates
        if isfield(P,'ySlc') && isfield(P,'zSlc') 
            % plot at center between two holes for hole-based geometry
            if ~isfield(P,'celltype') 
                ctrHoleRX = P.xc + ds.P.hx_hole(1)/2;   % central hole right-x-coordinate
                adjHoleLX = P.xc + ds.P.a_hole(2) - ds.P.hx_hole(2)/2; % adjacent hole left-x-coordinate
                P.xSlc = 0.5*(ctrHoleRX+adjHoleLX)*P.holeatctr;
            else
%                 P.xSlc = 0;
            end
            
            slc = [P.xSlc;P.ySlc;P.zSlc];
            for i = 1:length(P.zSiV)
                zSiVstr = num2str(P.zSiV{i});
                zSiVstr = strrep(zSiVstr,'-','m');
                zSiVstr = strrep(zSiVstr,' ','');
                pSlcSiV = repmat(slc,size(ds.mfem.locInd));
                pSlc.(['z',zSiVstr]).pos = abs(pSlcSiV);
            end
            
            PlotStrCplSiVXY_v20200527(ds,model,mModes,pSlc,P.zSiV,datLoc)
            PlotStrCplSiVYZ_v20200527(ds,model,mModes,pSlc,P.zSiV,datLoc)
        end
    end
    
    %% Save file
    mphsave(model,[datLoc,fileName,'.mph']);
    save([datLoc,fileName,'.mat'],'ds');
catch lasterror % if geometry cannot be fabbed
    disp(lasterror.message)
    disp('ERROR stack:')
    for lsi=1:length(lasterror.stack)
        display([num2str(lsi),': line ',...
                 num2str(lasterror.stack(lsi).line),...
                 ' in file: ',lasterror.stack(lsi).file])
    end
    
    
    mphsave(model,[datLoc,'error_',P.fileBase,'.mph']);
    save([datLoc,'error_',P.fileBase,'.mat'],'ds');
end

tEnd = toc(tStart);
disp(['Simulation time = ',num2str(tEnd/60,'%.2f'),' mins'])
disp(['Files saved in ',datLoc])
close all
