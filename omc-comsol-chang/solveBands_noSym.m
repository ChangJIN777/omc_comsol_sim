function ds = solveBands_noSym(P)
%SOLVEBANDS Summary of this function goes here
%   Detailed explanation goes here
if ~strcmp(P.datLoc(end),'\')
    datLoc = [P.datLoc,'\'];
else
    datLoc = P.datLoc;
end
% create directory to save files
if ~exist(datLoc,'dir')
    mkdir(datLoc)
end

% create base folder name (TODO: need to update)
if ~isfield(P,'fileBase')
    if strcmp(P.celltype,'rib')
        fBase = ['hole_unitCell','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            's_',num2str(P.s*1e9,'%.0f'),'nm',...
            'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
            't_',num2str(P.t*1e9,'%.0f'),'nm_', ...
            'th_',num2str(P.th*1e9,'%.0f'),'nm_', ...
            'd_',num2str(P.d*180/pi,'%.0f'),'deg'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    end
end
fBase = P.fileBase;
if isempty(dir([datLoc,fBase,'_bds.mat']))
    tStart = tic;
    disp('Solving with no symmetric boundary condition');
    full = runBands_noSym(P);

    %find gaps
    [full.midGap,full.gapSize] = findGaps(full);
    
    %% test
    % write to data structure
    ds.full = full;
    
    %% find complete bandgaps
    
    [full.midGap,full.gapSize] = findGaps(full);
    ds.full = full;
    
    % save data structures
    if P.savedat
        fname_mat = [fBase,'_bds.mat'];
        pathMat = [P.datLoc,fname_mat];
        save(pathMat,'ds');
    end
    
    %% plot bandstructure
    if P.savebndplot
       figure; hold on
       maxFreqs = [0 0 0 0];

       p1 = plot(full.k_norm,full.F*1e-9,'-k','linewidth',2,'DisplayName','full');

        if P.bandStruct_2D

            % plot full bandgaps
            if ~isempty(full.midGap)
                for k = 1:length(full.gapSize)
                    bgp = patch([0 3 3 0],1e-9.*(full.midGap(k) + 0.5*[full.gapSize(k) ...
                        full.gapSize(k) -full.gapSize(k) -full.gapSize(k)]),180/255*[1 0 0],'EdgeColor','none');
                    alpha(bgp,0.2);
                end
            end

            % plot complete midgap frequencies
            for k = 1:length(full.midGap)
                midfreqs = full.midGap(k)*ones(length(full.k_norm),1);
                plot(full.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
            end

            xlabel('k','FontSize',12);
            ylabel('Frequency (GHz)','FontSize',12);
            amax = max([full.F(:)])*1e-9;
            axis([0 3 0 amax]);
            %         axis tight
            set(gca,'XTick',[0; 1; 2; 3]);
            %         set(gca,'XTickLabel',{'G','C'},'fontname','symbol','fontsize',16)
            set(gca,'XTickLabel',{'\Gamma','X','M','\Gamma'},'fontsize',12)
        else 

            % plot full bandgaps
            if ~isempty(full.midGap)
                for k = 1:length(full.gapSize)
                    bgp = patch([0 1 1 0],1e-9.*(full.midGap(k) + 0.5*[full.gapSize(k) ...
                        full.gapSize(k) -full.gapSize(k) -full.gapSize(k)]),180/255*[1 0 0],'EdgeColor','none');
                    alpha(bgp,0.2);
                end
            end

            % plot complete midgap frequencies
            for k = 1:length(full.midGap)
                midfreqs = full.midGap(k)*ones(length(full.k_norm),1);
                plot(full.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
            end

            xlabel('k','FontSize',12);
            ylabel('Frequency (GHz)','FontSize',12);
            amax = max([full.F(:)])*1e-9;
            axis([0 1 0 amax]);
            %         axis tight
            set(gca,'XTick',[0; 1]);
            %         set(gca,'XTickLabel',{'G','C'},'fontname','symbol','fontsize',16)
            set(gca,'XTickLabel',{'\Gamma','X'},'fontsize',12)
            
        if strcmp(P.celltype,'rib')
            bandtitle = {['hole_unitCell_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            's_',num2str(P.s*1e9,'%.0f'),'nm',...
            'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
            't_',num2str(P.t*1e9,'%.0f'),'nm_', ...
            'th_',num2str(P.th*1e9,'%.0f'),'nm_', ...
            'd_',num2str(P.d*180/pi,'%.0f'),'deg']};
        else 
            bandtitle = 'Default plot';
        end
        title(bandtitle);
        box on
        hold off
       
        end
        % save band diagram as .png and .fig
        pathFig = [P.datLoc,fBase,'_fullBands'];
        saveas(gcf,[pathFig,'.png']);
        saveas(gcf,[pathFig,'.fig']);
    end
    tEnd = toc(tStart);
    disp(['Simulation time = ',num2str(tEnd/60,'%.2f'),'mins'])
else
    disp('Data folder exists in working directory')
    ds.full = [];
    ds.afull = [];
    
end

end
