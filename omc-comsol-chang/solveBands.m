function ds = solveBands(P)
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
    if strcmp(P.celltype,'2D_ribs')
        fBase = ['2D_ribs_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'w_',num2str(P.w*1e9,'%.0f'),'nm',...
            'hi_',num2str(P.hi*1e9,'%.0f'),'nm_', ...
            'wi_',num2str(P.wi*1e9,'%.0f'),'nm_', ...
            'th_',num2str(P.th*1e9,'%.0f'),'nm_', ...
            'ho_',num2str(P.ho*1e9,'%.0f'),'nm_',...
            'wo_',num2str(P.wo*1e9,'%.0f'),'nm_',...
            'ai_',num2str(P.ai*1e9,'%.0f'),'nm_'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    elseif strcmp(P.celltype,'hole')
        fBase = ['hole_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'hx_',num2str(P.hx*1e9,'%.0f'),'nm',...
            'hy_',num2str(P.hy*1e9,'%.0f'),'nm_', ...
            'beamWidth_',num2str(P.beam_width*1e9,'%.0f'),'nm_', ...
            'th_',num2str(P.th*1e9,'%.0f'),'nm_'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    elseif strcmp(P.celltype,'cross')
        fBase = ['cross_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'h_',num2str(P.h*1e9,'%.0f'),'nm',...
            'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
            'th_',num2str(P.th*1e9,'%.0f'),'nm_', ...
            'r1_',num2str(P.r1*1e9,'%.0f'),'nm_', ...s
            'r2_',num2str(P.r2*1e9,'%.0f'),'nm_'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    elseif strcmp(P.celltype,'boomerang')
        fBase = ['boomerang_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'r_',num2str(P.r*1e9,'%.0f'),'nm',...
            'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
            'th_',num2str(P.th*1e9,'%.0f'),'nm_', ...
            'r1_',num2str(P.r1*1e9,'%.0f'),'nm_', ...
            'r2_',num2str(P.r2*1e9,'%.0f'),'nm_'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    elseif strcmp(P.celltype,'snowflake')
        fBase = ['snowflake_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'r_',num2str(P.r*1e9,'%.0f'),'nm',...
            'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
            'th_',num2str(P.th*1e9,'%.0f'),'nm_', ...
            'r1_',num2str(P.r1*1e9,'%.0f'),'nm_', ...
            'r2_',num2str(P.r2*1e9,'%.0f'),'nm_'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    elseif strcmp(P.celltype,'boomerang_lower')
        fBase = ['boomerang_lower_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'r_',num2str(P.r*1e9,'%.0f'),'nm',...
            'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
            'd_',num2str(P.d*1e9,'%.0f'),'nm_', ...
            'h_',num2str(P.h*1e9,'%.0f'),'nm_', ...
            'th_',num2str(P.th*1e9,'%.0f'),'nm_', ...
            'r1_',num2str(P.r1*1e9,'%.0f'),'nm_', ...
            'r2_',num2str(P.r2*1e9,'%.0f'),'nm_'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    elseif strcmp(P.celltype,'boomerang_strip')
        fBase = ['boomerang_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'r_',num2str(P.r*1e9,'%.0f'),'nm',...
            'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
            'th_',num2str(P.th*1e9,'%.0f'),'nm_', ...
            'wo_',num2str(P.wo*1e9,'%.0f'),'nm_', ...
            'wi_',num2str(P.wi*1e9,'%.0f'),'nm_', ...
            'ho_',num2str(P.ho*1e9,'%.0f'),'nm_', ...
            'hi_',num2str(P.hi*1e9,'%.0f'),'nm_', ...
            'r1_',num2str(P.r1*1e9,'%.0f'),'nm_', ...
            'r2_',num2str(P.r2*1e9,'%.0f'),'nm_'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    elseif strcmp(P.celltype,'Snowflake_strip')
        fBase = ['boomerang_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'r_',num2str(P.r*1e9,'%.0f'),'nm',...
            'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
            'wo_',num2str(P.wo*1e9,'%.0f'),'nm_', ...
            'wi_',num2str(P.wi*1e9,'%.0f'),'nm_', ...
            'ho_',num2str(P.ho*1e9,'%.0f'),'nm_', ...
            'hi_',num2str(P.hi*1e9,'%.0f'),'nm_', ...
            'd_',num2str(P.d*1e9,'%.0f'),'nm_', ...
            'th_',num2str(P.th*1e9,'%.0f'),'nm_', ...
            'r1_',num2str(P.r1*1e9,'%.0f'),'nm_', ...
            'r2_',num2str(P.r2*1e9,'%.0f'),'nm_'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    elseif strcmp(P.celltype,'boomerang_strip_v2')
        fBase = ['boomerang_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'r_',num2str(P.r*1e9,'%.0f'),'nm',...
            'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
            'wo_',num2str(P.wo*1e9,'%.0f'),'nm_', ...
            'wi_',num2str(P.wi*1e9,'%.0f'),'nm_', ...
            'ho_',num2str(P.ho*1e9,'%.0f'),'nm_', ...
            'hi_',num2str(P.hi*1e9,'%.0f'),'nm_', ...
            'd_',num2str(P.d*1e9,'%.0f'),'nm_', ...
            'th_',num2str(P.th*1e9,'%.0f'),'nm_', ...
            'r1_',num2str(P.r1*1e9,'%.0f'),'nm_', ...
            'r2_',num2str(P.r2*1e9,'%.0f'),'nm_'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    end
end
fBase = P.fileBase;
%% set up the symmetry conditions 
if isempty(dir([datLoc,fBase,'_bds.mat']))
    tStart = tic;
    if P.TwoSymPlanes
        disp('Solving symmetric y and symmetric z');
        P.mbevenz = 1;
        P.mbeveny = 1;
        if P.bandStruct_2D
            symy_symz = runBands_2D(P);
        else 
            symy_symz = runBands(P);
        end
        %find gaps
        [symy_symz.midGap,symy_symz.gapSize] = findGaps(symy_symz);

        disp('Solving symmetric y and asymmetric z');
        P.mbevenz = -1;
        P.mbeveny = 1;
        if P.bandStruct_2D
            symy_asymz = runBands_2D(P);
        else 
            symy_asymz = runBands(P);
        end
        %find gaps
        [symy_asymz.midGap,symy_asymz.gapSize] = findGaps(symy_asymz);
        
        disp('Solving asymmetric y and symmetric z');
        P.mbevenz = 1;
        P.mbeveny = -1;
        if P.bandStruct_2D
            asymy_symz = runBands_2D(P);
        else 
            asymy_symz = runBands(P);
        end
        
        %find gaps
        [asymy_symz.midGap,asymy_symz.gapSize] = findGaps(asymy_symz);

        disp('Solving asymmetric y and asymmetric z');
        P.mbevenz = -1;
        P.mbeveny = -1;
        if P.bandStruct_2D
            asymy_asymz = runBands_2D(P);
        else 
            asymy_asymz = runBands(P);
        end
        
        %find gaps
        [asymy_asymz.midGap,asymy_asymz.gapSize] = findGaps(asymy_asymz);
    else
        disp('Solving with symmetric boundary condition');
        if P.zSymCondition
            disp('Solving with z symmetric boundary condition');
            P.mbevenz = 1;
        else 
            disp('Solving with y symmetric boundary condition');
            P.mbeveny = 1;
        end
        if P.bandStruct_2D
            sym = runBands_2D(P);
        else 
            sym = runBands(P);
        end
        
        %find gaps
        [sym.midGap,sym.gapSize] = findGaps(sym);
        
        if P.solveasym
            disp('Solving with anti-symmetric boundary condition');
            if P.zSymCondition
                P.mbevenz = -1;
            else
                P.mbeveny = -1;
            end
            if P.bandStruct_2D
                asym = runBands_2D(P);
            else
                asym = runBands(P);
            end

            %find gaps
            [asym.midGap,asym.gapSize] = findGaps(asym);
        end
    end
    %% test
    % write to data structure
    if P.TwoSymPlanes
        ds.symy_symz = symy_symz;
        ds.asymy_symz = asymy_symz;
        ds.symy_asymz = symy_asymz;
        ds.asymy_asymz = asymy_asymz;
    else
        ds.sym = sym;
        if P.solveasym
            ds.asym = asym;
        end
    end

    %% find complete bandgaps
    if P.TwoSymPlanes
        full.F = [symy_symz.F,asymy_symz.F,symy_asymz.F,asymy_asymz.F];
    else
        if P.solveasym
            full.F = [sym.F,asym.F];
        else
            full.F = sym.F;
        end
    end

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
            
            if P.TwoSymPlanes
                p1 = plot(symy_symz.k_norm,symy_symz.F*1e-9,'-k','linewidth',2,'DisplayName','symy_symz');
                p2 = plot(symy_asymz.k_norm,symy_asymz.F*1e-9,'--b','linewidth',2,'DisplayName','symy_asymz');
                p3 = plot(asymy_symz.k_norm,asymy_symz.F*1e-9,'--r','linewidth',2,'DisplayName','asymy_symz');
                p4 = plot(asymy_asymz.k_norm,asymy_asymz.F*1e-9,'--m','linewidth',2,'DisplayName','asymy_asymz');
            else
                p1 = plot(sym.k_norm,sym.F*1e-9,'-k','linewidth',2,'DisplayName','sym');
                if P.solveasym
                    p2 = plot(asym.k_norm,asym.F*1e-9,'--b','linewidth',2,'DisplayName','asym');
                end
            end

        if P.bandStruct_2D

            % plot symmetric bandgaps
            for k = 1:length(symy_symz.gapSize)
                bgp = patch([0 3 3 0],1e-9.*(symy_symz.midGap(k) + 0.5*[symy_symz.gapSize(k) ...
                    symy_symz.gapSize(k) -symy_symz.gapSize(k) -symy_symz.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
                alpha(bgp,0.5);
            end

            % plot asymmetric bandgaps
            for k = 1:length(symy_asymz.gapSize)
                bgp = patch([0 3 3 0],1e-9.*(symy_asymz.midGap(k) + 0.5*[symy_asymz.gapSize(k) ...
                    symy_asymz.gapSize(k) -symy_asymz.gapSize(k) -symy_asymz.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
                alpha(bgp,0.2);
            end

            % plot asymmetric bandgaps
            for k = 1:length(asymy_symz.gapSize)
                bgp = patch([0 3 3 0],1e-9.*(asymy_asymz.midGap(k) + 0.5*[asymy_asymz.gapSize(k) ...
                    asymy_asymz.gapSize(k) -asymy_asymz.gapSize(k) -asymy_asymz.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
                alpha(bgp,0.2);
            end

            % plot asymmetric bandgaps
            for k = 1:length(asymy_asymz.gapSize)
                bgp = patch([0 3 3 0],1e-9.*(asymy_asymz.midGap(k) + 0.5*[asymy_asymz.gapSize(k) ...
                    asymy_asymz.gapSize(k) -asymy_asymz.gapSize(k) -asymy_asymz.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
                alpha(bgp,0.2);
            end

            % plot full bandgaps
            if ~isempty(full.midGap)
                for k = 1:length(full.gapSize)
                    bgp = patch([0 3 3 0],1e-9.*(full.midGap(k) + 0.5*[full.gapSize(k) ...
                        full.gapSize(k) -full.gapSize(k) -full.gapSize(k)]),180/255*[1 0 0],'EdgeColor','none');
                    alpha(bgp,0.2);
                end
            end

            % plot symmetric midgap frequencies
            for k = 1:length(sym.midGap)
                midfreqs = sym.midGap(k)*ones(length(sym.k_norm),1);
                plot(sym.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
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
        else 
            if P.TwoSymPlanes
                % --- two symmetry planes (four structs) ---
                % plot symmetric bandgaps
                for k = 1:length(symy_symz.gapSize)
                    bgp = patch([0 1 1 0],1e-9.*(symy_symz.midGap(k) + 0.5*[symy_symz.gapSize(k) ...
                        symy_symz.gapSize(k) -symy_symz.gapSize(k) -symy_symz.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
                    alpha(bgp,0.5);
                end

                % plot asymmetric bandgaps
                for k = 1:length(asymy_asymz.gapSize)
                    bgp = patch([0 1 1 0],1e-9.*(asymy_asymz.midGap(k) + 0.5*[asymy_asymz.gapSize(k) ...
                        asymy_asymz.gapSize(k) -asymy_asymz.gapSize(k) -asymy_asymz.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
                    alpha(bgp,0.2);
                end

                % plot full bandgaps
                if ~isempty(full.midGap)
                    for k = 1:length(full.gapSize)
                        bgp = patch([0 1 1 0],1e-9.*(full.midGap(k) + 0.5*[full.gapSize(k) ...
                            full.gapSize(k) -full.gapSize(k) -full.gapSize(k)]),180/255*[1 0 0],'EdgeColor','none');
                        alpha(bgp,0.2);
                    end
                end

                % plot symmetric midgap frequencies
                for k = 1:length(symy_symz.midGap)
                    midfreqs = symy_symz.midGap(k)*ones(length(symy_symz.k_norm),1);
                    plot(symy_symz.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
                end

                % plot asymmetric midgap frequencies
                for k = 1:length(asymy_asymz.midGap)
                    midfreqs = asymy_asymz.midGap(k)*ones(length(asymy_asymz.k_norm),1);
                    plot(asymy_asymz.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
                end

                % plot complete midgap frequencies
                for k = 1:length(full.midGap)
                    midfreqs = full.midGap(k)*ones(length(symy_symz.k_norm),1);
                    plot(symy_symz.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
                end

                amax = max([symy_symz.F(:);asymy_asymz.F(:);symy_asymz.F(:);asymy_symz.F(:)])*1e-9;
            else
                % --- single symmetry plane (sym / asym) ---
                % plot symmetric bandgaps
                for k = 1:length(sym.gapSize)
                    bgp = patch([0 1 1 0],1e-9.*(sym.midGap(k) + 0.5*[sym.gapSize(k) ...
                        sym.gapSize(k) -sym.gapSize(k) -sym.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
                    alpha(bgp,0.5);
                end

                % plot asymmetric bandgaps
                if P.solveasym
                    for k = 1:length(asym.gapSize)
                        bgp = patch([0 1 1 0],1e-9.*(asym.midGap(k) + 0.5*[asym.gapSize(k) ...
                            asym.gapSize(k) -asym.gapSize(k) -asym.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
                        alpha(bgp,0.2);
                    end
                end

                % plot full bandgaps
                if ~isempty(full.midGap)
                    for k = 1:length(full.gapSize)
                        bgp = patch([0 1 1 0],1e-9.*(full.midGap(k) + 0.5*[full.gapSize(k) ...
                            full.gapSize(k) -full.gapSize(k) -full.gapSize(k)]),180/255*[1 0 0],'EdgeColor','none');
                        alpha(bgp,0.2);
                    end
                end

                % plot symmetric midgap frequencies
                for k = 1:length(sym.midGap)
                    midfreqs = sym.midGap(k)*ones(length(sym.k_norm),1);
                    plot(sym.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
                end

                % plot asymmetric midgap frequencies
                if P.solveasym
                    for k = 1:length(asym.midGap)
                        midfreqs = asym.midGap(k)*ones(length(asym.k_norm),1);
                        plot(asym.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
                    end
                end

                % plot complete midgap frequencies
                for k = 1:length(full.midGap)
                    midfreqs = full.midGap(k)*ones(length(sym.k_norm),1);
                    plot(sym.k_norm,midfreqs*1e-9,'.--r','linewidth',0.5);
                end

                if P.solveasym
                    amax = max([sym.F(:);asym.F(:)])*1e-9;
                else
                    amax = max(sym.F(:))*1e-9;
                end
            end

            xlabel('k','FontSize',12);
            ylabel('Frequency (GHz)','FontSize',12);
            axis([0 1 0 amax]);
            %         axis tight
            set(gca,'XTick',[0; 1]);
            %         set(gca,'XTickLabel',{'G','C'},'fontname','symbol','fontsize',16)
            set(gca,'XTickLabel',{'\Gamma','X'},'fontsize',12)
            
        if strcmp(P.celltype,'2D_ribs')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.w*1e9,'%.0f'),'nm',...
                'hi = ',num2str(P.hi*1e9,'%.0f'),'nm, ',...
                'wi = ',num2str(P.wi*1e9,'%.0f'),'nm',...
                'ho = ',num2str(P.ho*1e9,'%.0f'),'nm, ',...
                'wo = ',num2str(P.wo*1e9,'%.0f'),'nm',...
                'ai = ',num2str(P.ai*1e9,'%.0f'),'nm'];
                ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'hole')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'hx = ',num2str(P.hx*1e9,'%.0f'),'nm, ',...
                'hy = ',num2str(P.hy*1e9,'%.0f'),'nm',...
                'beam_width = ',num2str(P.beam_width*1e9,'%.0f'),'nm'];
                ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'cross')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'h = ',num2str(P.h*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.wc*1e9,'%.0f'),'nm',...
                'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
                'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm'];
                ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'boomerang')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
                'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
                'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm'];
                ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'boomerang_lower')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
                'd = ',num2str(P.d*1e9,'%.0f'),'nm',...
                'h = ',num2str(P.h*1e9,'%.0f'),'nm',...
                'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
                'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm'];
                ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'snowflake')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
                'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
                'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm'];
                ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'Snowflake_strip')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
                'wo = ',num2str(P.wo*1e9,'%.0f'),'nm',...
                'wi = ',num2str(P.wi*1e9,'%.0f'),'nm',...
                'ho = ',num2str(P.ho*1e9,'%.0f'),'nm',...
                'hi = ',num2str(P.hi*1e9,'%.0f'),'nm',...
                'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
                'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm'];
                ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'boomerang_strip')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
                'wo = ',num2str(P.wo*1e9,'%.0f'),'nm',...
                'wi = ',num2str(P.wi*1e9,'%.0f'),'nm',...
                'ho = ',num2str(P.ho*1e9,'%.0f'),'nm',...
                'hi = ',num2str(P.hi*1e9,'%.0f'),'nm',...
                'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
                'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm'];
                ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'boomerang_strip_v2')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
                'wo = ',num2str(P.wo*1e9,'%.0f'),'nm',...
                'wi = ',num2str(P.wi*1e9,'%.0f'),'nm',...
                'ho = ',num2str(P.ho*1e9,'%.0f'),'nm',...
                'hi = ',num2str(P.hi*1e9,'%.0f'),'nm',...
                'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
                'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm'];
                ['th = ',num2str(P.th*1e9,'%.0f'),'nm']};
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
    ds.sym = [];
    ds.asym = [];
    
end

end

