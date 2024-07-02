% solves bandstructure along Gamma-X-M wavevectors for 2D unit cell


function ds = solveBands2D(P,datLoc)
%append \ to end of datLoc if not present
if ~strcmp(datLoc(end),'\')
    datLoc = [datLoc,'\'];
end
% create directory to save files
if ~exist(datLoc,'dir')
    mkdir(datLoc)
end

fBase = [P.xsect,'_',P.celltype,'X_',...
    'a_',num2str(P.a*1e9,'%.0f'),'nm_',...
    'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
    't_',num2str(P.th*1e9,'%.0f'),'nm_', ...
    'cx_',num2str(P.cx*1e9,'%.0f'),'nm_', ...
    'cy_',num2str(P.cy*1e9,'%.0f'),'nm_',...
    'tx_',num2str(P.tx*1e9,'%.0f'),'nm_', ...
    'ty_',num2str(P.ty*1e9,'%.0f'),'nm'];
if isfield(P,'prefname')
    fBase = [P.prefname,'_',fBase];
end
if P.filRad > 0
    fBase = [fBase,'_Rfil_',num2str(P.filRad*1e9,'%.1f'),'nm'];
end
display(fBase)

if isempty(dir([datLoc,fBase,'*.mat']))
    display('Solving with symmetric boundary condition');
    P.mbevenz = 1;
    sym = RunCrossBands(P,datLoc);
    
    %find gaps
    [sym.midGap,sym.gapSize] = findGaps(sym);
    
    if P.solveasym == 1
        display('Solving with anti-symmetric boundary condition');
        P.mbevenz = -1;
        asym = RunCrossBands(P,datLoc);
        
        %find gaps
        [asym.midGap,asym.gapSize] = findGaps(asym);
    end
    
    %% test
	% write to data structure
    ds.sym = sym;
    if P.solveasym == 1
        ds.asym = asym;
    end
    
    % save data structures
    if P.savedat
        fname_mat = [fBase,'_bds.mat'];
        pathMat = [datLoc,fname_mat];
        save(pathMat,'ds');
    end

    %% 
    % plot bandstructure
    if P.savebndplot
        figure; hold on
        plot(sym.k_norm,sym.F*1e-9,'-k','linewidth',3);
        if P.solveasym == 1
            plot(asym.k_norm,asym.F*1e-9,'--b','linewidth',1);
        end
        
        % plot symmetric bandgaps
        for k = 1:length(sym.gapSize)
            bgp = patch([0 3 3 0],1e-9.*(sym.midGap(k) + 0.5*[sym.gapSize(k) ...
                sym.gapSize(k) -sym.gapSize(k) -sym.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
            alpha(bgp,0.5);
        end
        
        if P.solveasym == 1
            % plot asymmetric bandgaps
            for k = 1:length(asym.gapSize)
                bgp = patch([0 3 3 0],1e-9.*(asym.midGap(k) + 0.5*[asym.gapSize(k) ...
                    asym.gapSize(k) -asym.gapSize(k) -asym.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
                alpha(bgp,0.4);
            end
        end
        
        % plot symmetric midgap frequencies
        for k = 1:length(sym.midGap)
            midfreqs = sym.midGap(k)*ones(length(sym.k_norm),1);
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
        
        if P.filRad > 0
            filRadTtl = [', Rf = ',num2str(P.filRad*1e9,'%.1f'),'nm'];
        else
            filRadTtl = '';
        end
        bandtitle = {[P.celltype,' cross, ',...
                      'a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                      'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
                      'th = ',num2str(P.th*1e9,'%.0f'),'nm'];
                     ['cx = ',num2str(P.cx*1e9,'%.0f'),'nm, ', ...
                      'cy = ',num2str(P.cy*1e9,'%.0f'),'nm, ',...
                      'tx = ',num2str(P.tx*1e9,'%.0f'),'nm, ',...
                      'ty = ',num2str(P.ty*1e9,'%.0f'),'nm',...
                      filRadTtl]};
        title(bandtitle);
        box on
        hold off
        
        % save band diagram as .png and .fig
        pathFig = [datLoc,fBase,'_fullBands'];
        saveas(gcf,[pathFig,'.png']);
        saveas(gcf,[pathFig,'.fig']);
    end

else
    display('Data folder exists in working directory')
    ds.sym = [];
    ds.asym = [];
    
end

