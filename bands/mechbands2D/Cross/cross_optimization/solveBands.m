function ds = solveBands(vars,P)
%SOLVEBANDS Summary of this function goes here
%   Detailed explanation goes here


P.beamMat = 'diamond';                  % beam material name
P.anisoMat = 1;clc
P.rxtal = 45;  

P.a = vars.a*1e-9;%580e-9;                           % nominal lattice constant
P.hc = vars.hc*1e-9;%929e-9;                           % beam width
P.wc = vars.wc*1e-9;%250e-9;                          % nominal hole height (along x-axis)

P.th = 160e-9; %This should be modified
P.r1 = 10e-9;
P.r2 = 10e-9;

P.freq = 0;                           % target frequency - set to 0 for bandstructure simulations
P.kpts = 1;                             % no. of k-points, EXCLUDING gamma point
P.nbands = 6;                           % no. of bands to solve for
P.max_dof = 3e6;                        % max # of degrees of freedom
P.meshSize = 5;

P.saveplots = 0;
P.saveMPH = 0;
P.savebndplot = 1;
P.savedat = 0;

if ~strcmp(P.datLoc(end),'\')
    datLoc = [P.datLoc,'\'];
else
    datLoc = P.datLoc;
end
% create directory to save files
if ~exist(datLoc,'dir')
    mkdir(datLoc)
end

% create base folder name
if ~isfield(P,'fileBase')
    fBase = ['a_',num2str(P.a*1e9,'%.0f'),'nm_',...
        'hc_',num2str(P.hc*1e9,'%.0f'),'nm_', ...
        'wc_',num2str(P.wc*1e9,'%.0f'),'nm_', ...
        'th_',num2str(P.th*1e9,'%.0f'),'nm_', ...
        'r1_',num2str(P.r1*1e9,'%.0f'),'nm_',...
        'r2_',num2str(P.r2*1e9,'%.0f'),'nm_'];
    if isfield(P,'prefname')
        fBase = [P.prefname,'_',fBase];
    end
    P.fileBase = fBase;
end
fBase = P.fileBase;
if isempty(dir([datLoc,fBase,'_bds.mat']))
    tStart = tic;
    disp('Solving with symmetric boundary condition');
    P.mbevenz = 1;
    sym = runBands(P);
    
    %find gaps
    [sym.midGap,sym.gapSize] = findGaps(sym);
    
    disp('Solving with anti-symmetric boundary condition');
    P.mbevenz = -1;
    asym = runBands(P);
    
    %find gaps
    [asym.midGap,asym.gapSize] = findGaps(asym);
    %% test
    % write to data structure
    ds.sym = sym;
    ds.asym = asym;
    
    %% find complete bandgaps

    full.F = [sym.F,asym.F];
    
    [full.midGap,full.gapSize] = findGaps(full);
    ds.full = full;
    
    % save data structures
    if P.savedat
        fname_mat = [fBase,'_bds.mat'];
        pathMat = [P.datLoc,fname_mat];
        save(pathMat,'ds');
    end
    
    %%
    % plot bandstructure
    if P.savebndplot
        figure; hold on
        plot(sym.k_norm,sym.F*1e-9,'-k','linewidth',3);
        plot(asym.k_norm,asym.F*1e-9,'--b','linewidth',1);
        
        % plot asymmetric bandgaps
        for k = 1:length(asym.gapSize)
            bgp = patch([0 3 3 0],1e-9.*(asym.midGap(k) + 0.5*[asym.gapSize(k) ...
                asym.gapSize(k) -asym.gapSize(k) -asym.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
            alpha(bgp,0.4);
        end
        
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
            midfreqs = full.midGap(k)*ones(length(sym.k_norm),1);
            plot(sym.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
        end
        
        xlabel('k','FontSize',12);
        ylabel('Frequency (GHz)','FontSize',12);
        amax = max([sym.F(:);asym.F(:)])*1e-9;
        axis([0 3 0 amax]);
        %         axis tight
        set(gca,'XTick',[0; 1; 2; 3]);
        %         set(gca,'XTickLabel',{'G','C'},'fontname','symbol','fontsize',16)
        set(gca,'XTickLabel',{'\Gamma','X','M','\Gamma'},'fontsize',12)
        
        
        bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
            'hc = ',num2str(P.hc*1e9,'%.0f'),'nm, ',...
            'wc = ',num2str(P.wc*1e9,'%.0f'),'nm'];
            ['th = ',num2str(P.th*1e9,'%.0f'),'nm, ', ...
            'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm, ',...
            'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm, ']};
        title(bandtitle);
        box on
        hold off
        
        % save band diagram as .png and .fig
        pathFig = [P.datLoc,fBase,'_fullBands'];
        saveas(gcf,[pathFig,'.png']);
        saveas(gcf,[pathFig,'.fig']);
    end
    tEnd = toc(tStart);
    disp(['Simulation time = ',num2str(tEnd/60,'%.2f'),'mins'])
else
    disp('Data folder exists in working directory')
    ds.sym = [];
    ds.asym = [];
    
end

end

