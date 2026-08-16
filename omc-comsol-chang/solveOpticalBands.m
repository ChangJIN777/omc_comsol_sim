function ds = solveOpticalBands(P)

%% make the file directory
% Normalise the separators before testing or creating anything, then write the
% result back into P so every downstream consumer sees the same path. The old
% version tested strcmp(P.datLoc(end),'\') against a hard-coded backslash: a
% datLoc ending in '/' failed that test, so a stray '\' was appended and the
% directory that actually got created was not the one the figures were later
% written into. That is the "Unable to create output file ... No such file or
% directory" you get from the plotgeom block in runOpticalBand_*, which fires
% before this function's own saveplots mkdir further down.
%
% Same normalise-and-create pattern runBands_2D already uses. On Windows filesep
% is '\', so this is a no-op there and the behaviour is unchanged.
if isfield(P,'datLoc') && ~isempty(P.datLoc)
    P.datLoc = strrep(strrep(P.datLoc,'\',filesep),'/',filesep);
    if ~strcmp(P.datLoc(end),filesep)
        P.datLoc = [P.datLoc,filesep];
    end
end
datLoc = P.datLoc;
% create directory to save files
if ~exist(datLoc,'dir')
    mkdir(datLoc)
end

% create base folder name (TODO: need to update)
if ~isfield(P,'fileBase')
    if strcmp(P.celltype,'2D_ribs')
        fBase = ['optical_2D_ribs_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'w_',num2str(P.w*1e9,'%.0f'),'nm',...
            'hi_',num2str(P.hi*1e9,'%.0f'),'nm_', ...
            'wi_',num2str(P.wi*1e9,'%.0f'),'nm_', ...
            'ho_',num2str(P.ho*1e9,'%.0f'),'nm_',...
            'wo_',num2str(P.wo*1e9,'%.0f'),'nm_',...
            'ai_',num2str(P.ai*1e9,'%.0f'),'nm_'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    elseif strcmp(P.celltype,'cross')
        fBase = ['optical_cross_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'hc_',num2str(P.hc*1e9,'%.0f'),'nm',...
            'wc_',num2str(P.wc*1e9,'%.0f'),'nm_', ...
            'r1_',num2str(P.r1*1e9,'%.0f'),'nm_', ...
            'r2_',num2str(P.r2*1e9,'%.0f'),'nm_'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    elseif strcmp(P.celltype,'boomerang')
        fBase = ['optical_boomerang_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'r_',num2str(P.r*1e9,'%.0f'),'nm',...
            'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
            'r1_',num2str(P.r1*1e9,'%.0f'),'nm_', ...
            'r2_',num2str(P.r2*1e9,'%.0f'),'nm_'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    elseif strcmp(P.celltype,'boomerang_strip_v2')
        fBase = ['boomerang_strip_optical_',...
            'a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'r_',num2str(P.r*1e9,'%.0f'),'nm',...
            'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
            'hi_',num2str(P.hi*1e9,'%.0f'),'nm_', ...
            'ho_',num2str(P.ho*1e9,'%.0f'),'nm_', ...
            'wi_',num2str(P.wi*1e9,'%.0f'),'nm_', ...
            'wo_',num2str(P.wo*1e9,'%.0f'),'nm_', ...
            'th_',num2str(P.th*1e9,'%.0f'),'nm_', ...
            'r1_',num2str(P.r1*1e9,'%.0f'),'nm_', ...
            'r2_',num2str(P.r2*1e9,'%.0f'),'nm_'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
     elseif strcmp(P.celltype,'boomerang_strip')
        fBase = ['boomerang_strip_optical_',...
            'a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'r_',num2str(P.r*1e9,'%.0f'),'nm',...
            'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
            'h_',num2str(P.h*1e9,'%.0f'),'nm_', ...
            'd_',num2str(P.d1*1e9,'%.0f'),'nm_', ...
            'th_',num2str(P.th*1e9,'%.0f'),'nm_', ...
            'r1_',num2str(P.r1*1e9,'%.0f'),'nm_', ...
            'r2_',num2str(P.r2*1e9,'%.0f'),'nm_'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    elseif strcmp(P.celltype,'snowflake')
        fBase = ['optical_snowflake_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'r_',num2str(P.r*1e9,'%.0f'),'nm',...
            'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
            'r1_',num2str(P.r1*1e9,'%.0f'),'nm_', ...
            'r2_',num2str(P.r2*1e9,'%.0f'),'nm_'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
        
    elseif strcmp(P.celltype,'hole')
        fBase = ['optical_hole_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'hx_',num2str(P.hx*1e9,'%.0f'),'nm',...
            'hy_',num2str(P.hy*1e9,'%.0f'),'nm',...
            'beam_width_',num2str(P.beam_width*1e9,'%.0f'),'nm',...
            'th_',num2str(P.th*1e9,'%.0f'),'nm'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    elseif strcmp(P.celltype,'hole_strip')
        fBase = ['optical_holeStrip_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'r_',num2str(P.r*1e9,'%.0f'),'nm',...
            'b_',num2str(P.b*1e9,'%.0f'),'nm'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    elseif strcmp(P.celltype,'hole_strip_wvg')
        fBase = ['optical_holeStrip_wvg_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'r_',num2str(P.r*1e9,'%.0f'),'nm',...
            'b_',num2str(P.b*1e9,'%.0f'),'nm'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    elseif strcmp(P.celltype,'Snowflake_strip_1d')
        fBase = ['optical_snowflake_strip_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
            'b_',num2str(P.b*1e9,'%.0f'),'nm_',...
            'b_wvg_',num2str(P.b_wvg*1e9,'%.0f'),'nm_',...
            'r_',num2str(P.r*1e9,'%.0f'),'nm_',...
            'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
            'r1_',num2str(P.r1*1e9,'%.0f'),'nm_', ...
            'r2_',num2str(P.r2*1e9,'%.0f'),'nm_', ...
            'th_',num2str(P.th*1e9,'%.0f'),'nm_'];
        if isfield(P,'prefname')
            fBase = [P.prefname,'_',fBase];
        end
        P.fileBase = fBase;
    elseif strcmp(P.celltype,'rib')
        fBase = ['optical_rib_','a_',num2str(P.a*1e9,'%.0f'),'nm_',...
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

%% solve the optical band structure
if isempty(dir([datLoc,fBase,'_bds.mat']))
    tStart = tic;
    if P.bandStructureDim==3
        OpticalBand = runOpticalBand_3D(P);
    elseif P.bandStructureDim==2
        OpticalBand = runOpticalBand_2D(P);
    else 
        OpticalBand = runOpticalBand_1D(P);
    end
    
    % Decide the k-path shape ONCE, here, because both the light line and the
    % plot branch further down depend on it and they must not disagree.
    % bandStruct_2D is the flag the runners themselves use to decide whether
    % ds.k_norm is the 0->3 circuit or is collapsed to kx_norm.
    % runOpticalBand_1D collapses unconditionally and never reads it, hence the
    % dimension test as well.
    use2Dpath = ismember(P.bandStructureDim,[2 3]) && ...
                isfield(P,'bandStruct_2D') && P.bandStruct_2D;

    % filter data below light line
    %
    % The light line used for the FILTER is the same one the plot DRAWS, so that
    % the modes feeding findGaps_optical are exactly the modes a reader sees
    % below the blue curves. On the 2D circuit that means reproducing the three
    % per-segment expressions from the plot branch below, evaluated at each
    % k_norm rather than on a 100-point display grid:
    %
    %   Gamma-X  (0<=q<1):  t = q/2,      f = t
    %   X-M      (1<=q<2):  t = (q-1)/2,  f = sqrt(1/4 + (t/sqrt(3))^2)
    %   M-Gamma  (2<=q<=3): t = (q-2)/2,  f = sqrt((1/2-t)^2 + (1/(2*sqrt(3))-t/sqrt(3))^2)
    %
    % all times c/a. The three agree at the joins (0.5 at X, 0.5774 at M, 0 at
    % the closing Gamma), so the assembled curve is continuous.
    %
    % NOTE: this is NOT the same as c*hypot(kx_norm,ky_norm)/(2a), which is what
    % this line computed before. The two agree on Gamma-X and at X, but diverge
    % along X-M - at M they differ by 2/sqrt(3) = 1.1547. They encode different
    % conventions for where the high-symmetry points sit, and only one can be
    % right. Matching the drawn curve is what was asked for and makes the figure
    % self-consistent; see the note reported alongside this change before
    % trusting either as absolute physics.
    % The formula itself lives in opticalLightLine so test_plotOpticalBands can
    % reproduce this filter exactly against saved data, instead of keeping a
    % second copy that quietly diverges.
    lightline = opticalLightLine(OpticalBand,P);
    TE.F0 = OpticalBand.F;
    TEbelow = OpticalBand.F < lightline;    % check which bands are below lightline
    TE.F = OpticalBand.F.*TEbelow;        % filter out data below lightline
    TE.F(TE.F==0) = NaN;      % replace zeros with NaN so they don't get plotted

    % find gaps
    %
    % findGaps_belowLightLine, NOT findGaps_optical. The question being asked
    % here is "which frequency regions contain no guided data point", and that
    % is a statement about the SET of surviving frequencies - it does not depend
    % on which band a point was assigned to. That independence matters: the
    % runners assign a point to a band purely by its position in COMSOL's
    % eigenvalue list at that k-point, with no mode tracking and no sort, so
    % column identity is not a safe thing to build on. findGaps_optical does
    % build on it (it collapses each COLUMN to an interval). It is left
    % untouched on disk but now has no callers anywhere - the mechanical
    % pipeline uses findGaps, which is a different function.
    %
    % Two independent criteria, both applied inside the search so a region that
    % fails is never constructed in the first place:
    %
    %   P.minOpticalGap    minimum WIDTH,  default 20 THz. A gap is measured
    %                      between sampled k-points only - with kpts = 9 the
    %                      circuit has 27 samples - so a narrow region may just
    %                      be the spacing between two samples of one band rather
    %                      than a real void in the spectrum.
    %   P.minOpticalMidGap minimum CENTRE, default 100 THz. Width alone does not
    %                      separate a usable gap from the wide, sparsely
    %                      populated stretches low in the spectrum, where few
    %                      modes survive the light-line mask and the voids
    %                      between them are large but nowhere near the design
    %                      target.
    if isfield(P,'minOpticalGap') && ~isempty(P.minOpticalGap)
        minOpticalGap = P.minOpticalGap;
    else
        minOpticalGap = 20e12;   % Hz
    end
    if isfield(P,'minOpticalMidGap') && ~isempty(P.minOpticalMidGap)
        minOpticalMidGap = P.minOpticalMidGap;
    else
        minOpticalMidGap = 100e12;   % Hz
    end

    [OpticalBand.midGap,OpticalBand.gapSize,OpticalBand.gapEdges,maxRejected] = ...
        findGaps_belowLightLine(TE.F,minOpticalGap,minOpticalMidGap);

    if isempty(OpticalBand.gapSize)
        if isempty(maxRejected)
            fprintf('  no optical gap found below the light line\n');
        else
            fprintf(['  no optical gap meeting width >= %.2f THz AND centre ' ...
                     '>= %.2f THz;\n  widest region found was %.2f THz ' ...
                     '(it may have failed on centre, not width)\n'], ...
                minOpticalGap*1e-12,minOpticalMidGap*1e-12,maxRejected*1e-12);
        end
    else
        fprintf('  %d optical gap(s): width >= %.2f THz, centre >= %.2f THz\n', ...
            numel(OpticalBand.gapSize),minOpticalGap*1e-12,minOpticalMidGap*1e-12);
    end
    % write to data structure
    ds.opticalBand = OpticalBand;

    % Everything needed to re-derive the gaps WITHOUT re-solving travels with
    % the result. solveBands gets this for free - its sub-structs come straight
    % out of runBands and its gap criterion has no parameters - but the optical
    % gaps are the product of a light line and a minimum-gap threshold, and a
    % file that records the answer without recording the criterion cannot be
    % compared against a later run scored under different rules.
    ds.P             = P;              % geometry, k-path flags, backend choice
    ds.lightline     = lightline;      % Hz at each k-point, as filtered AND drawn
    ds.bandsBelow    = TE.F;           % NaN-masked matrix findGaps_optical scored
    ds.minOpticalGap    = minOpticalGap;    % Hz minimum width applied
    ds.minOpticalMidGap = minOpticalMidGap; % Hz minimum centre applied
    ds.use2Dpath     = use2Dpath;      % which k-path shape the above assumed
    ds.solveTimeMin  = toc(tStart)/60; % elapsed at save time, not at return

    % Persist, mirroring solveBands.m:248-252 on the mechanical side. Without
    % this the cache test at the top of this block - which looks for exactly
    % this file - could never hit, because nothing in the optical chain ever
    % wrote one: a run left behind figures and a returned struct and nothing
    % else, so every repeat re-solved from scratch and no band data survived to
    % be re-analysed later.
    %
    % Deliberately BEFORE the plotting block: plotting is the part most likely
    % to throw (an unlisted celltype, a headless session, a full disk), and a
    % crash there must not cost the hours of solve time that produced the data.
    % isfield-guarded rather than solveBands' bare P.savedat so a config that
    % predates the flag warns instead of erroring.
    if ~isfield(P,'savedat')
        warning('solveOpticalBands:noSavedat', ...
            ['P.savedat is not set, so no _bds.mat will be written and this ' ...
             'solve will have to be repeated to get the data back. Set ' ...
             'P.savedat = 1 in the test script.']);
    elseif P.savedat
        pathMat = [P.datLoc,fBase,'_bds.mat'];
        save(pathMat,'ds');
        fprintf('  saved band data -> %s\n',pathMat);
    end
%% plot bandstructure
if P.savebndplot
    % Which PLOT to draw is set by the shape of the k-path, not by the model
    % dimension. Those are different questions and this used to conflate them:
    % the condition was bandStructureDim==2, so a dim = 3 run - a 3D slab model
    % sweeping the full hexagonal Gamma-X-M-Gamma circuit - fell into the 1D
    % branch below. That branch fixes the axis at [0 1] while k_norm runs 0 to 3,
    % so only the Gamma-X third was ever visible, under two tick labels taken
    % from a four-label list.
    %
    % use2Dpath was decided up at the light-line filter and is reused verbatim
    % here. Computing it once is deliberate: the filter and the plot must agree
    % about the shape of the k-path, or the figure shows one light line while the
    % gaps were scored against another.
    if use2Dpath
        figure; hold on
        maxFreqs = [0 0 0 0];

        % Markers only, no connecting line. Connecting consecutive k-points
        % would assert a band continuity the data does not carry: a point is
        % assigned to a column by its position in COMSOL's eigenvalue list at
        % that k-point, with no mode tracking, so a column is a frequency LEVEL
        % rather than a single physical mode and the line would swap modes at
        % every crossing.
        p1 = plot(OpticalBand.k_norm,OpticalBand.F*1e-12,'ko','linewidth',2,'DisplayName','sym','MarkerSize',5);
        % plot the light lines 
        lightx = linspace(0,0.5,100);
        lighty1 = lightx*(3e8)/(P.a*(1e12));
        lighty2 = sqrt(1/4+(lightx./sqrt(3)).^2).*3e8./(P.a*(1e12));
        lighty3 = sqrt((0.5-lightx).^2 + (0.5/sqrt(3)-(lightx./sqrt(3))).^2)*(3e8)/(P.a*(1e12));
        hold on;
        light1 = plot(2*lightx,lighty1,'b-','linewidth',1);
        light2 = plot(1+2*lightx,lighty2,'b-','linewidth',1);
        light3 = plot(2+2*lightx,lighty3,'b-','linewidth',1);
        
        % plot optical bandgaps
        for k = 1:length(OpticalBand.gapSize)
            bgp = patch([0 3 3 0],(1e-12)*(OpticalBand.midGap(k) + 0.5*[OpticalBand.gapSize(k) ...
                OpticalBand.gapSize(k) -OpticalBand.gapSize(k) -OpticalBand.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
            alpha(bgp,0.5);
        end
        
        % plot symmetric midgap frequencies
        for k = 1:length(OpticalBand.midGap)
            midfreqs = OpticalBand.midGap(k)*ones(length(OpticalBand.k_norm),1);
            plot(OpticalBand.k_norm,midfreqs*1e-12,'.--r','linewidth',0.5);
        end
      
        xlabel('k','FontSize',12);
        ylabel('Frequency (THz)','FontSize',12);
        amax = max([OpticalBand.F(:)])*1e-12;
        axis([0 3 0 amax]);
        %         axis tight
        % One tick per high-symmetry point, matching the four XTickLabels below.
        % This was [0; 1], so MATLAB consumed only the first two labels and the
        % circuit was annotated as if it ran Gamma to X.
        set(gca,'XTick',[0; 1; 2; 3]);
        %         set(gca,'XTickLabel',{'G','C'},'fontname','symbol','fontsize',16)
        set(gca,'XTickLabel',{'\Gamma','X','M','\Gamma'},'fontsize',12);
        
        if strcmp(P.celltype,'2D_ribs')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.w*1e9,'%.0f'),'nm',...
                'hi = ',num2str(P.hi*1e9,'%.0f'),'nm, ',...
                'wi = ',num2str(P.wi*1e9,'%.0f'),'nm',...
                'ho = ',num2str(P.ho*1e9,'%.0f'),'nm, ',...
                'wo = ',num2str(P.wo*1e9,'%.0f'),'nm',...
                'ai = ',num2str(P.ai*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'cross')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'hc = ',num2str(P.hc*1e9,'%.0f'),'nm, ',...
                'wc = ',num2str(P.wc*1e9,'%.0f'),'nm',...
                'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
                'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'boomerang_strip_v2')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
                'hi = ',num2str(P.hi*1e9,'%.0f'),'nm, ',...
                'wi = ',num2str(P.wi*1e9,'%.0f'),'nm',...
                'ho = ',num2str(P.ho*1e9,'%.0f'),'nm, ',...
                'wo = ',num2str(P.wo*1e9,'%.0f'),'nm',...
                'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
                'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'boomerang_lower')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
                'd = ',num2str(P.d*1e9,'%.0f'),'nm',...
                'h = ',num2str(P.h*1e9,'%.0f'),'nm',...
                'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
                'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'snowflake')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
                'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
                'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'hole')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm, ',...
                'b = ',num2str(P.b*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'hole_strip')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm, ',...
                'b = ',num2str(P.b*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'rib')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                's_',num2str(P.s*1e9,'%.0f'),'nm',...
                'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
                't_',num2str(P.t*1e9,'%.0f'),'nm_', ...
                'th_',num2str(P.th*1e9,'%.0f'),'nm_', ...
                'd_',num2str(P.d*180/pi,'%.0f'),'deg']};
        else
            % Fallback so an unlisted celltype cannot leave bandtitle undefined
            % and take the title() call down with it. 'boomerang' reaches this
            % branch now that dim = 3 plots the 2D circuit, and previously there
            % was neither a case for it nor an else. Only a, r, w are named
            % because they are the fields every celltype in this chain shares.
            bandtitle = {[P.celltype,': a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'r_',num2str(P.r*1e9,'%.0f'),'nm_',...
                'w_',num2str(P.w*1e9,'%.0f'),'nm']};
        end
        title(bandtitle);
        box on
        hold off
        
        % save band diagram as .png and .fig
        pathFig = [P.datLoc,fBase,'_fullBands'];
        saveas(gcf,[pathFig,'.png']);
        saveas(gcf,[pathFig,'.fig']);
    else 
        figure; hold on
        maxFreqs = [0 0 0 0];

        % Markers only, no connecting line. Connecting consecutive k-points
        % would assert a band continuity the data does not carry: a point is
        % assigned to a column by its position in COMSOL's eigenvalue list at
        % that k-point, with no mode tracking, so a column is a frequency LEVEL
        % rather than a single physical mode and the line would swap modes at
        % every crossing.
        p1 = plot(OpticalBand.k_norm,OpticalBand.F*1e-12,'ko','linewidth',2,'DisplayName','sym','MarkerSize',5);
        
        % plot the light line 
        hold on;
        lightx = linspace(0,1,100);
        lighty1 = lightx*(3e8)/(P.a*(1e12))/2;
        light1 = plot(OpticalBand.kx_norm,lightline*1e-12,'b-','linewidth',1);
        
        % plot optical bandgaps
        for k = 1:length(OpticalBand.gapSize)
            bgp = patch([0 1 1 0],(1e-12)*(OpticalBand.midGap(k) + 0.5*[OpticalBand.gapSize(k) ...
                OpticalBand.gapSize(k) -OpticalBand.gapSize(k) -OpticalBand.gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
            alpha(bgp,0.5);
        end
        
        % plot symmetric midgap frequencies
        for k = 1:length(OpticalBand.midGap)
            midfreqs = OpticalBand.midGap(k)*ones(length(OpticalBand.k_norm),1);
            plot(OpticalBand.k_norm,midfreqs*1e-12,'.--r','linewidth',0.5);
        end
      
        xlabel('k','FontSize',12);
        ylabel('Frequency (THz)','FontSize',12);
        amax = max([OpticalBand.F(:)])*1e-12;
        axis([0 1 0 amax]);
        %         axis tight
        set(gca,'XTick',[0; 1]);
        %         set(gca,'XTickLabel',{'G','C'},'fontname','symbol','fontsize',16)
        % Two labels for two ticks. This branch draws a 1D path that ends at X,
        % so the M and Gamma labels it used to carry were both surplus and
        % misleading - MATLAB silently dropped them, labelling k = 1 as X, which
        % happens to be right, but only by accident.
        set(gca,'XTickLabel',{'\Gamma','X'},'fontsize',12)
        
        if strcmp(P.celltype,'2D_ribs')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.w*1e9,'%.0f'),'nm',...
                'hi = ',num2str(P.hi*1e9,'%.0f'),'nm, ',...
                'wi = ',num2str(P.wi*1e9,'%.0f'),'nm',...
                'ho = ',num2str(P.ho*1e9,'%.0f'),'nm, ',...
                'wo = ',num2str(P.wo*1e9,'%.0f'),'nm',...
                'ai = ',num2str(P.ai*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'cross')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'hc = ',num2str(P.hc*1e9,'%.0f'),'nm, ',...
                'wc = ',num2str(P.wc*1e9,'%.0f'),'nm',...
                'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
                'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'boomerang_strip_v2')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
                'hi = ',num2str(P.hi*1e9,'%.0f'),'nm, ',...
                'wi = ',num2str(P.wi*1e9,'%.0f'),'nm',...
                'ho = ',num2str(P.ho*1e9,'%.0f'),'nm, ',...
                'wo = ',num2str(P.wo*1e9,'%.0f'),'nm',...
                'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
                'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'snowflake')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
                'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
                'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'boomerang_lower')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
                'd = ',num2str(P.d*1e9,'%.0f'),'nm',...
                'h = ',num2str(P.h*1e9,'%.0f'),'nm',...
                'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
                'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'hole')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'hx = ',num2str(P.hx*1e9,'%.0f'),'nm, ',...
                'hy = ',num2str(P.hy*1e9,'%.0f'),'nm, ',...
                'th = ',num2str(P.th*1e9,'%.0f'),'nm, ',...
                'beamWidth = ',num2str(P.beam_width*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'hole_strip')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm, ',...
                'b = ',num2str(P.b*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'hole_strip_wvg')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm, ',...
                'b = ',num2str(P.b*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'Snowflake_strip_2d')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                'b = ',num2str(P.b*1e9,'%.0f'),'nm, ',...
                'bbase = ',num2str(P.b_base*1e9,'%.0f'),'nm, ',...
                'w = ',num2str(P.w*1e9,'%.0f'),'nm, ',...
                'r = ',num2str(P.r*1e9,'%.0f'),'nm',...
                'r1 = ',num2str(P.r1*1e9,'%.0f'),'nm',...
                'r2 = ',num2str(P.r2*1e9,'%.0f'),'nm',...
                'th = ',num2str(P.th*1e9,'%.0f'),'nm']};
        elseif strcmp(P.celltype,'rib')
            bandtitle = {['a = ',num2str(P.a*1e9,'%.0f'),'nm, ',...
                's_',num2str(P.s*1e9,'%.0f'),'nm',...
                'w_',num2str(P.w*1e9,'%.0f'),'nm_', ...
                't_',num2str(P.t*1e9,'%.0f'),'nm_', ...
                'th_',num2str(P.th*1e9,'%.0f'),'nm_', ...
                'd_',num2str(P.d*180/pi,'%.0f'),'deg']}; 
        else 
            bandtitle = 'blank title';
        end
        title(bandtitle);
        box on
        hold off
        
        % save band diagram as .png and .fig
        pathFig = [P.datLoc,fBase,'_fullBands'];
        saveas(gcf,[pathFig,'.png']);
        saveas(gcf,[pathFig,'.fig']);
end
end
    tEnd = toc(tStart);
    disp(['Simulation time = ',num2str(tEnd/60,'%.2f'),'mins'])
else
    % Cache hit. Two things were wrong here and both mattered once the solve
    % branch above started actually writing the file it tests for:
    %
    %   1. The field was spelled ds.OpticalBand - capital O - while the solve
    %      branch returns ds.opticalBand. A caller testing
    %      isfield(ds,'opticalBand') therefore got true on a fresh run and
    %      false on a cached one, for the same design.
    %   2. It returned empty rather than the data sitting in the file whose
    %      existence is the entire reason this branch was taken. solveBands
    %      does the same, but it had the excuse that nothing had been saved to
    %      load; that is no longer true here.
    pathMat = [datLoc,fBase,'_bds.mat'];
    fprintf('Cached optical band data found, loading instead of re-solving:\n  %s\n', ...
        pathMat);
    loaded = load(pathMat,'ds');
    if isfield(loaded,'ds') && isfield(loaded.ds,'opticalBand')
        ds = loaded.ds;
    else
        % A file written before the wrapper existed, or by something else
        % entirely. Report it rather than returning a struct that silently
        % lacks the field every caller expects.
        warning('solveOpticalBands:staleCache', ...
            ['%s exists but carries no ds.opticalBand, so it cannot be ' ...
             'reused. Delete it to force a fresh solve.'],pathMat);
        ds.opticalBand = [];
    end
end