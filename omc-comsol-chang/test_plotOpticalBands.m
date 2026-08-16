function results = test_plotOpticalBands(datLocs,varargin)
%TEST_PLOTOPTICALBANDS Re-analyse and plot SAVED optical band structures.
%
% Loads *_bds.mat files written by solveOpticalBands, re-applies the
% below-light-line filter and the gap search, draws the band diagram with the
% light line and gap shading, and prints every gap to the terminal.
%
% Runs NO simulation. COMSOL and LiveLink are not needed, so this is the cheap
% way to re-score a finished sweep - try a different minimum gap, or check what
% a run actually produced - without spending solve time.
%
% Gaps are recomputed from the stored band matrix rather than read from the
% file, so data saved before the light-line and minimum-gap changes is scored
% by today's rules, and so MinGap can be swept without re-running anything.
% Whatever midGap/gapSize the file already carries is reported alongside, and a
% disagreement is flagged rather than hidden.
%
% INPUTS
%   datLocs   Directory, or cell array of directories, to search. A directory
%             is searched one level deep, so passing the parent of a set of
%             dated run folders works. Defaults to the location
%             test_Boomerang.m writes to: fullfile('.','test','boomerang').
%
% NAME-VALUE OPTIONS
%   'MinGap'   Minimum gap to report, in Hz. Default 1e13 (10 THz), matching
%              solveOpticalBands. Set 0 to see every gap found.
%   'Pattern'  File pattern to match. Default '*_bds.mat'.
%   'Plot'     true (default) to draw a figure per file, false for text only.
%   'Verbose'  true (default) to print the per-file header and gap table.
%
% OUTPUT
%   results   Struct array, one entry per file loaded, with fields:
%             file, datLoc, a, midGap, gapSize, midGapTHz, gapSizeTHz,
%             gapRatio (gapSize/midGap), nBandsBelowLightLine, storedGapSize.
%             Empty struct if nothing matched.
%
% EXAMPLES
%   % everything under the default results root
%   test_plotOpticalBands;
%
%   % one dated run, no minimum, text only
%   test_plotOpticalBands(fullfile('.','test','boomerang','08152026'), ...
%                         'MinGap',0,'Plot',false);
%
%   % collect for further analysis
%   R = test_plotOpticalBands({'./test/boomerang','./test/boomerang_bayesopt'});
%   [~,best] = max([R.gapRatio]);
%
% Requires only base MATLAB. Depends on opticalLightLine.m and
% findGaps_belowLightLine.m from this directory.

%% ---------------------------------------------------------------- arguments
if nargin < 1 || isempty(datLocs)
    datLocs = fullfile('.','test','boomerang');
end
if ~iscell(datLocs)
    datLocs = {char(datLocs)};
end

opt = struct('MinGap',1e13,'Pattern','*_bds.mat','Plot',true,'Verbose',true);
if mod(numel(varargin),2) ~= 0
    error('test_plotOpticalBands:badOptions', ...
        'Name-value options must come in pairs; got %d trailing argument(s).', ...
        numel(varargin));
end
for ii = 1:2:numel(varargin)
    name = validatestring(varargin{ii},fieldnames(opt), ...
        'test_plotOpticalBands','option name');
    opt.(name) = varargin{ii+1};
end
validateattributes(opt.MinGap,{'numeric'},{'scalar','nonnegative','finite'}, ...
    'test_plotOpticalBands','MinGap');

%% ------------------------------------------------------------ collect files
files = {};
for ii = 1:numel(datLocs)
    root = datLocs{ii};
    if ~exist(root,'dir')
        warning('test_plotOpticalBands:noDir','Not a directory, skipped: %s',root);
        continue
    end
    % the directory itself, then one level of subdirectories - dated run
    % folders sit one level under the geometry folder
    found = dir(fullfile(root,opt.Pattern));
    files = [files; fullfileCell(found)];                          %#ok<AGROW>
    subs  = dir(root);
    subs  = subs([subs.isdir] & ~ismember({subs.name},{'.','..'}));
    for jj = 1:numel(subs)
        found = dir(fullfile(root,subs(jj).name,opt.Pattern));
        files = [files; fullfileCell(found)];                      %#ok<AGROW>
    end
end

results = struct([]);
if isempty(files)
    fprintf(['No files matching %s under:\n'],opt.Pattern);
    fprintf('    %s\n',datLocs{:});
    fprintf(['\nIf the optical runs are recent, note that solveOpticalBands only\n' ...
             'started writing _bds.mat once P.savedat is set - runs from before\n' ...
             'that left figures but no band data. Check P.savedat = 1.\n']);
    return
end

fprintf('Found %d file(s) matching %s\n\n',numel(files),opt.Pattern);

%% ------------------------------------------------------------ analyse each
keep = 0;
for ii = 1:numel(files)
    fpath = files{ii};
    [thisDir,thisName] = fileparts(fpath);

    S = load(fpath);
    [OB,P,why] = extractOpticalBand(S);
    if isempty(OB)
        warning('test_plotOpticalBands:notOptical', ...
            'Skipped %s: %s',thisName,why);
        continue
    end

    % --- light line and below-light-line filter, identical to solveOpticalBands
    lightline = opticalLightLine(OB,P);
    below     = OB.F < lightline;          % implicit expansion over bands
    TE.F      = OB.F.*below;
    TE.F(TE.F==0) = NaN;

    % findGaps_belowLightLine, matching solveOpticalBands. It reports frequency
    % regions containing no guided data point, which is a property of the SET of
    % surviving frequencies and so is immune to how points were assigned to
    % bands - the runners do that by eigenvalue index with no mode tracking.
    [midGap,gapSize,~,maxRejected] = findGaps_belowLightLine(TE.F,opt.MinGap);

    % --- report
    keep = keep + 1;
    results(keep).file    = fpath;                                 %#ok<AGROW>
    results(keep).datLoc  = thisDir;                               %#ok<AGROW>
    results(keep).a       = P.a;                                   %#ok<AGROW>
    results(keep).midGap  = midGap;                                %#ok<AGROW>
    results(keep).gapSize = gapSize;                               %#ok<AGROW>
    results(keep).midGapTHz  = midGap*1e-12;                       %#ok<AGROW>
    results(keep).gapSizeTHz = gapSize*1e-12;                      %#ok<AGROW>
    results(keep).gapRatio   = gapSize./midGap;                    %#ok<AGROW>
    results(keep).nBandsBelowLightLine = nnz(any(below,1));        %#ok<AGROW>
    if isfield(OB,'gapSize')
        results(keep).storedGapSize = OB.gapSize(:);               %#ok<AGROW>
    else
        results(keep).storedGapSize = [];                          %#ok<AGROW>
    end

    if opt.Verbose
        printOneResult(thisName,thisDir,P,OB,below,midGap,gapSize, ...
                       maxRejected,opt.MinGap,results(keep).storedGapSize);
    end

    if opt.Plot
        plotOneResult(thisName,P,OB,lightline,below,midGap,gapSize);
    end
end

if keep == 0
    fprintf('No usable optical band data in the files found.\n');
    return
end

%% --------------------------------------------------------------- summary
fprintf('%s\n',repmat('=',1,72));
fprintf('SUMMARY: %d file(s) with optical band data, MinGap = %.2f THz\n', ...
    keep,opt.MinGap*1e-12);
fprintf('%-44s %8s %9s %8s\n','design','gaps','best THz','ratio');
fprintf('%s\n',repmat('-',1,72));
for ii = 1:numel(results)
    [~,nm] = fileparts(results(ii).file);
    nm = strrep(nm,'_bds','');
    if numel(nm) > 43, nm = ['...',nm(end-39:end)]; end
    if isempty(results(ii).gapSize)
        fprintf('%-44s %8d %9s %8s\n',nm,0,'-','-');
    else
        [bestGap,bi] = max(results(ii).gapSize);
        fprintf('%-44s %8d %9.3f %8.4f\n',nm,numel(results(ii).gapSize), ...
            bestGap*1e-12,results(ii).gapRatio(bi));
    end
end
fprintf('%s\n',repmat('=',1,72));
end

%% ======================================================== local functions

function c = fullfileCell(d)
%FULLFILECELL dir() struct array -> cell array of full paths (column).
if isempty(d)
    c = {};
    return
end
c = arrayfun(@(x) fullfile(x.folder,x.name),d(:),'UniformOutput',false);
end

function [OB,P,why] = extractOpticalBand(S)
%EXTRACTOPTICALBAND Pull the optical band struct and its P out of a loaded file.
%
% Handles the shapes a _bds.mat can take: the optical wrapper written by
% solveOpticalBands (ds.opticalBand), a bare runner struct saved directly, and
% mechanical files (ds.full / ds.sym), which are rejected with a reason rather
% than half-analysed.
OB = []; P = []; why = '';

if ~isfield(S,'ds')
    why = 'no variable ds in file';
    return
end
ds = S.ds;

if isfield(ds,'opticalBand') && isstruct(ds.opticalBand)
    OB = ds.opticalBand;
elseif isfield(ds,'F') && isfield(ds,'k_norm')
    OB = ds;                     % a runner struct saved without the wrapper
elseif isfield(ds,'full') || isfield(ds,'sym')
    why = 'mechanical band data (ds.full/ds.sym), not optical';
    return
else
    why = 'ds has neither opticalBand nor F/k_norm';
    return
end

if ~isfield(OB,'F') || isempty(OB.F)
    OB = []; why = 'optical struct carries no band matrix F';
    return
end

% P travels with the results (runOpticalBand_*.m sets ds.P = P), which is what
% makes re-analysis possible without the original test script.
if isfield(OB,'P') && isstruct(OB.P)
    P = OB.P;
elseif isfield(ds,'P') && isstruct(ds.P)
    P = ds.P;
else
    OB = []; why = 'no P struct saved with the results, cannot rebuild the light line';
    return
end

% Fields opticalLightLine needs. Older files predate bandStructureDim; infer it
% from the k-path rather than guessing a default, since the two light-line
% forms differ.
if ~isfield(P,'bandStructureDim')
    if isfield(OB,'k_norm') && max(OB.k_norm(:)) > 1.5
        P.bandStructureDim = 3;
    else
        P.bandStructureDim = 1;
    end
end
if ~isfield(P,'bandStruct_2D')
    P.bandStruct_2D = isfield(OB,'k_norm') && max(OB.k_norm(:)) > 1.5;
end
end

function printOneResult(name,dirName,P,OB,below,midGap,gapSize,maxRejected, ...
                        minGap,storedGapSize)
%PRINTONERESULT Terminal report for a single saved band structure.
fprintf('%s\n',repmat('-',1,72));
fprintf('%s\n',name);
fprintf('  dir           : %s\n',dirName);

geomBits = {};
flds = {'a','w','r','th','r1','r2'};
for ii = 1:numel(flds)
    if isfield(P,flds{ii}) && isscalar(P.(flds{ii}))
        geomBits{end+1} = sprintf('%s=%.0fnm',flds{ii},P.(flds{ii})*1e9); %#ok<AGROW>
    end
end
if ~isempty(geomBits)
    fprintf('  geometry      : %s\n',strjoin(geomBits,'  '));
end

fprintf('  band matrix   : %d k-points x %d bands\n',size(OB.F,1),size(OB.F,2));
fprintf('  below light   : %d of %d bands have any point below the light line\n', ...
    nnz(any(below,1)),size(OB.F,2));

if isempty(gapSize)
    fprintf('  GAPS          : none >= %.2f THz',minGap*1e-12);
    if ~isempty(maxRejected)
        fprintf(' (widest region found was %.2f THz)',maxRejected*1e-12);
    end
    fprintf('\n');
else
    fprintf('  GAPS          : %d >= %.2f THz',numel(gapSize),minGap*1e-12);
    if ~isempty(maxRejected)
        fprintf(' (widest sub-threshold region %.2f THz)',maxRejected*1e-12);
    end
    fprintf('\n');
    fprintf('     %3s  %12s  %12s  %10s  %12s\n', ...
        '#','midgap THz','size THz','gap/mid','midgap nm');
    for k = 1:numel(gapSize)
        lam_nm = 299792458/midGap(k)*1e9;
        fprintf('     %3d  %12.4f  %12.4f  %10.4f  %12.1f\n', ...
            k,midGap(k)*1e-12,gapSize(k)*1e-12,gapSize(k)/midGap(k),lam_nm);
    end
end

% Recomputed vs stored: a mismatch means the file was written under different
% filter rules, which is worth knowing before comparing across old and new runs.
if ~isempty(storedGapSize)
    sameCount = numel(storedGapSize) == numel(gapSize);
    if sameCount && ~isempty(gapSize)
        sameCount = max(abs(sort(storedGapSize)-sort(gapSize))) < 1e6;  % 1 MHz
    end
    if ~sameCount
        fprintf(['  NOTE          : file stores %d gap(s), recomputed %d. The saved\n' ...
                 '                  values were scored under different filter rules;\n' ...
                 '                  the table above is today''s.\n'], ...
            numel(storedGapSize),numel(gapSize));
    end
end
end

function plotOneResult(name,P,OB,lightline,below,midGap,gapSize)
%PLOTONERESULT Band diagram with light line and gap shading.
[~,use2Dpath] = opticalLightLine(OB,P);

figure('Name',name,'NumberTitle','off'); hold on

% Gap shading first so the bands draw over it.
if use2Dpath, xspan = [0 3 3 0]; else, xspan = [0 1 1 0]; end
for k = 1:numel(gapSize)
    bgp = patch(xspan,(1e-12)*(midGap(k) + 0.5*[gapSize(k) gapSize(k) ...
        -gapSize(k) -gapSize(k)]),180/255*[1 1 1],'EdgeColor','none');
    alpha(bgp,0.5);
end

% Every solved mode in grey, then the ones under the light line - the modes the
% gap search actually scores - in black on top. The contrast is the point: it
% shows at a glance how much of the band structure the filter discarded.
kk = OB.k_norm(:);
Fabove = OB.F; Fabove(below)  = NaN;
Fbelow = OB.F; Fbelow(~below) = NaN;
% Markers only, no connecting line - see the note at the matching plot call in
% solveOpticalBands: a column is a frequency level, not a tracked mode, so a
% line between consecutive k-points would imply a continuity the data lacks.
plot(kk,Fabove*1e-12,'o','Color',[0.72 0.72 0.72],'MarkerSize',4);
plot(kk,Fbelow*1e-12,'ko','MarkerSize',5,'LineWidth',1.2);

% Light line on a dense grid, through the same function used for the filter.
if use2Dpath
    dense.k_norm = linspace(0,3,400)';
else
    dense.kx_norm = linspace(min(OB.kx_norm),max(OB.kx_norm),400)';
    dense.ky_norm = zeros(400,1);
    dense.k_norm  = dense.kx_norm;
end
plot(dense.k_norm,opticalLightLine(dense,P)*1e-12,'b-','LineWidth',1.2);

for k = 1:numel(midGap)
    plot(kk,midGap(k)*ones(size(kk))*1e-12,'.--r','LineWidth',0.5);
end

xlabel('k','FontSize',12);
ylabel('Frequency (THz)','FontSize',12);
amax = max(OB.F(:))*1e-12;
if ~isfinite(amax) || amax <= 0, amax = 1; end
if use2Dpath
    axis([0 3 0 amax]);
    set(gca,'XTick',[0; 1; 2; 3]);
    set(gca,'XTickLabel',{'\Gamma','X','M','\Gamma'},'fontsize',12);
else
    axis([0 1 0 amax]);
    set(gca,'XTick',[0; 1]);
    set(gca,'XTickLabel',{'\Gamma','X'},'fontsize',12);
end

ttl = sprintf('a = %.0f nm',P.a*1e9);
if isfield(P,'w'), ttl = [ttl,sprintf(', w = %.0f nm',P.w*1e9)]; end
if isfield(P,'r'), ttl = [ttl,sprintf(', r = %.0f nm',P.r*1e9)]; end
if ~isempty(gapSize)
    [bg,bi] = max(gapSize);
    ttl = [ttl,sprintf('  |  best gap %.2f THz (ratio %.3f)', ...
        bg*1e-12,bg/midGap(bi))];
else
    ttl = [ttl,'  |  no gap above threshold'];
end
title(ttl);
box on
hold off
end
