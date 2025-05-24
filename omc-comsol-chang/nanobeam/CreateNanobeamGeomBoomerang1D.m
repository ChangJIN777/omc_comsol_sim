% Function to create arrays of geometry parameters for nanobeam with holes
% Modified from createNanobeamCavity from FDTD scripts by Michael Burek
% Input:
% P: data structure consisting of following geometry parameters
%   a: lattice constant;  w: beam width; th: beam height;
%	hx: hole height; d: hole width; 
%   nholes: total no. of holes in 1/2 of beam;
%   ndef: # of holes in 1/2 the defect region;
%   maxdef: max defect ratio; 
%   oblong: hy is scaled by (1-maxdef)^(1+oblong), hx by (1-maxdef)^(1-oblong), 
%   e.g. oblong == 1 results in constant height but changing width in defect
%
% Output:
% P: data structure with updated geometry parameters
%    a_hole, hx_hole, hy_hole: list of lattice constants, hole heights
%    and hole widths
%    beamLenHalf, beamLen: lengths of half and full beams
%    geomHalf, geom: list of hole heights, hole widths, x- and y-
%    positions of holes for half and full beams
%
% Last updated by Cleaven Chia, 20170324

function [P] = CreateNanobeamGeomBoomerang1D(P)

% constants
a = P.a;        % lattice constant 
w = P.w;        % the width of the hole 
th = P.th;        % the height of the hole
r = P.r;        % the height of the hole 
wo = P.wo;        % the height of the hole in the lower portion
wi = P.wi;        % the height of the hole in the lower portion
ho = P.ho;        % the height of the hole in the lower portion
hi = P.hi;        % the height of the hole in the lower portion
d = P.d;        % the width of the hole in the lower portion
r1 = P.r1;      % the fillet radius of the edges of the hole 
r2 = P.r2;      % the fillet radius of the center of the hole 
nholes = P.nholes;
ndef = P.ndef;
maxdef = P.maxdef;
oblong = P.oblong;
holeatctr = P.holeatctr;
taperFunc = P.taperFunc;

%% cavity taper function
if strcmp(taperFunc,'quadratic')
    f = @(x) 1-x.^2;
elseif strcmp(taperFunc,'cubic')
    f = @(x) 1-3*x.^2+2*x.^3;
elseif strcmp(taperFunc,'linear')
    f = @(x) 1-x;
elseif strcmp(taperFunc,'none')
    f = @(x) 1;
else
    error('Invalid taper function specified')
end

%% end waveguide mirror taper function
if isfield(P,'wgmTaper') && isfield(P.wgmTaper,'func')
    if strcmp(P.wgmTaper.func,'quadratic')
        tf = @(x) 1-x.^2;
    elseif strcmp(P.wgmTaper.func,'cubic')
        tf = @(x) 1-3*x.^2+2*x.^3;
    elseif strcmp(P.wgmTaper.func,'linear')
        tf = @(x) 1-x;
    elseif strcmp(P.wgmTaper.func,'none')
        tf = @(x) 1;
    else
        error('Invalid waveguide mirror end taper function specified')
    end
end

% %% check if geometry is within fabrication tolerance
% if wid-hy < 150e-9
%     display(num2str(wid-hy));
%     error(['hy=',num2str(hy*1e9,'%.1f'),'nm, w=',num2str(wid*1e9,'%.1f'),...
%         'nm, w-hy=',num2str((wid-hy)*1e9,'%.1f'),'nm<200nm - gap too small for lithography']);
% % elseif ndef > 0 && (a*(1-maxdef)-hx*(1-maxdef)^(1-oblong) < 50e-9)
% elseif maxdef~= 1 && (a*(1-maxdef)-hx*(1-maxdef)^(1-oblong) < 50e-9)
%     error(['a_def=',num2str(a*(1-maxdef)*1e9,'%.1f'),'nm, hx_def=',num2str(hx*(1-maxdef)^(1-oblong)*1e9,'%.1f'),...
%         'nm, a_def-hx_def=',num2str((a*(1-maxdef)-hx*(1-maxdef)^(1-oblong))*1e9,'%.1f'),'nm<50nm - gap too small for lithography']);
% end
% 
% if (hx < 150e-9) || (hy < 150e-9)
%     error('Hole size too small for lithography')
% end


%% Definitions
%a_i = separation from ith hole to (i-1)th hole

% phononic mirror
if isfield(P,'phonMir') && P.nphonmir>0
    nPhM = P.nphonmir;
else
    nPhM = 0;
end

% end waveguide mirror taper
if isfield(P, 'wvgmir') && (P.wvgmir>0)
    nwgm = P.wvgmir;
else
    nwgm = 0;
end

% duplicate center hole
if isfield(P, 'ctrholeadd') && (P.ctrholeadd>0)
    nctr = P.ctrholeadd;
else
    nctr = 0;
end

% central hole - custom hole dimensions
if isfield(P, 'taperTo')
    if strcmp(P.taperTo,'custom')
        % we are only tapering h and d for the boomerang unit cells
        aT = P.a_ctr;
        hT = P.h_ctr;
        dT = P.d_ctr;
    else
        aT = 0;
        hT = 0;
        dT = 0;
    end
else
    aT = 0;
    hT = 0;
    dT = 0;
end

if holeatctr && isfield(P,'cavlen') && P.cavlen ~= 0
    cavlen = P.cavlen;
else
    cavlen = 0;
end

ntot = nholes + nwgm + nctr;
a_hole = zeros(1,ntot);
h_hole = zeros(1,ntot);
d_hole = zeros(1,ntot);

%% create defect and mirror cells in cavity
% max defect dimensions
% a_hole_maxdef = a*(1-maxdef);
% hx_hole_maxdef = hx*(1-maxdef)^(1-oblong);
% hy_hole_maxdef = hy*(1-maxdef)^(1+oblong);

for ki = 1:nholes+nctr
    % defect cells
    if ki <= ndef+nctr
        if ki <= nctr % duplicate central cavity cell
            k = 1; 
        else
            k = ki-nctr; % tapered defect cells
        end
        x = (k-1)/(ndef);
        h_hole(ki) = hT + (h-hT)*(1-maxdef*f(x))^(1-oblong);
        d_hole(ki) = dT + (d-dT)*(1-maxdef*f(x))^(1+oblong);
        
        % if hole at ctr, there is one fewer possible value for a compared
        % to h, d
        xa = (k-1-holeatctr)/(ndef-holeatctr);
        if xa < 0
            if holeatctr && ki > 1 
                a_hole(ki) = aT + (a-aT)*(1-maxdef);
            else
                % a_hole(ki) = NaN;
                a_hole(ki) = a;
            end
        else
            a_hole(ki) = aT + (a-aT)*(1-maxdef*f(xa)) + cavlen*(ki==2);
        end
    % mirror cells
    else
        h_hole(ki) = h;
        d_hole(ki) = d;
        a_hole(ki) = a;
    end
end

%% append phononic mirror, if defined
if nPhM>0 && isfield(P, 'phonMir')
    for k = 1:nPhM
        a_hole(nholes+nctr+k) = P.phonMir.a;
        h_hole(nholes+nctr+k) = P.phonMir.h;
        d_hole(nholes+nctr+k) = P.phonMir.d;
    end
end

%% append end waveguide mirror taper, if defined
if nwgm>0 && isfield(P, 'wgmTaper') %&& ~strcmp(P.wgmTaperFunc,'linearEndMaxDefect')
    % taper end - custom hole dimensions
    if isfield(P.wgmTaper, 'endtype')
        if strcmp(P.wgmTaper.endtype,'custom')
            aD = P.wgmTaper.a_end;
            hD = P.wgmTaper.h_end;
            dD = P.wgmTaper.d_end;
        elseif strcmp(P.wgmTaper.endtype,'maxdef')
            maxdef = P.wgmTaper.maxdef_end;
            oblong = P.wgmTaper.oblong_end;
            aD = a*(1-maxdef);
            hD = h*(1-maxdef)^(1-oblong);
            dD = d*(1-maxdef)^(1+oblong);
        else
            aD = 0;
            hD = 0;
            dD = 0;
        end
    else
        aD = 0;
        hD = 0;
        dD = 0;
    end
    
    nwvgmir = nwgm+1; %eliminate the first mirror hole and null hole at end of taper
    for k = 2:nwvgmir
        if strcmp(P.wgmTaper.endtype,'custom') || strcmp(P.wgmTaper.endtype,'maxdef')
            x = ((k-1)/(nwvgmir-1)); % not tapering to 0 hole size
        else
            x = ((k-1)/nwvgmir); % tapering to 0 hole size
        end
        
        a_hole(nholes+nctr+nPhM+k-1) = aD + (a-aD)*tf(x);
        h_hole(nholes+nctr+nPhM+k-1) = hD + (h-dD)*tf(x);
        d_hole(nholes+nctr+nPhM+k-1) = dD + (d-dD)*tf(x);
        
        % code for checking the fab compatibility of the structure
        % gapx = a_hole(nholes+nctr+nPhM+k-1) - 0.5*hx_hole(nholes+nctr+nPhM+k-1) ...
        %        - 0.5*hx_hole(nholes+nctr+nPhM+k-2);
        
        % % check if hole is within fabrication tolerance
        % if wid-d_hole(nholes+nctr+nPhM+k-1) < 100e-9
        %     display(num2str(wid-d_hole(nholes+nctr+nPhM+k-1))*1e9);
        %     error('Hole width in end taper too large relative to beam width');
        % elseif gapx < 50e-9
        %     display(gapx);
        %     error('Hole height in end taper too large relative to lattice constant');
        % end
        % 
        % if (hx_hole(nholes+nctr+nPhM+k-1) < 150e-9) ||...
        %    (d_hole(nholes+nctr+nPhM+k-1) < 150e-9)
        %     error('Hole size too small for lithograpd')
        % end
    end
    
end

%% assemble coordinates and diameters of holes
w_hole = a*ones(size(h_hole));

if holeatctr %&& (~isfield(P,'cavlen') || (isfield(P,'cavlen') && P.cavlen == 0))
    % for beam with hole at center of beam
    % half-beam
    xposHalf = [0,cumsum(a_hole(2:end))];   % excluding ctr hole
    yposHalf = zeros(size(a_hole));         % excluding ctr hole
    h_holeHalf = h_hole;
    d_holeHalf = d_hole;
    w_holeHalf = w_hole;
    a_holeHalf = a_hole;
    
    % full beam
    xposFull = [-1*fliplr(xposHalf),xposHalf(2:end)];
    yposFull = [yposHalf,yposHalf(2:end)];
    h_holeFull = [fliplr(h_hole) h_hole(2:end)];
    d_holeFull = [fliplr(d_hole) d_hole(2:end)];
    w_holeFull = [fliplr(w_hole) w_hole(2:end)];
    a_holeFull = [fliplr(a_hole) a_hole(2:end)];
    
% elseif holeatctr && (isfield(P,'cavlen') && P.cavlen > 0)
%     % for beam with hole at center of beam plus variable spacing to next
%     % hole
%     % half-beam
%     cavlen = P.cavlen;
%     xposHalf = [0,cumsum(a_hole(2:end))+cavlen];   % excluding ctr hole
%     yposHalf = zeros(size(a_hole));         % excluding ctr hole
%     hx_holeHalf = hx_hole;
%     hy_holeHalf = hy_hole;
%     
%     % full beam
%     xposFull = [-1*fliplr(xposHalf),xposHalf(2:end)];
%     yposFull = [yposHalf,yposHalf(2:end)];
%     hx_holeFull = [fliplr(hx_hole) hx_hole(2:end)];
%     hy_holeFull = [fliplr(hy_hole) hy_hole(2:end)];
    
elseif ~holeatctr && isfield(P,'cavlen') && P.cavlen > 0
    % for beam with cavity length in the middle of beam
    cavlen = P.cavlen;
    if (cavlen < 50e-9)
        display(num2str(cavlen));
        error('Cavity length too small to fabricate');
    end
    
    % half-beam
    xposHalf = cumsum([h_hole(1)/2+cavlen/2,a_hole(2:end)]);
    yposHalf = zeros(size(a_hole));
    h_holeHalf = h_hole;
    d_holeHalf = d_hole;
    w_holeHalf = w_hole;
    a_holeHalf = a_hole;
    
    % full beam
    xposFull = [-1*fliplr(xposHalf),xposHalf];
    yposFull = [yposHalf,yposHalf];
    h_holeFull = [fliplr(h_hole) h_hole];
    d_holeFull = [fliplr(d_hole) d_hole];
    w_holeFull = [fliplr(w_hole) w_hole];
    a_holeFull = [fliplr(a_hole) a_hole];
    
elseif ~holeatctr && isfield(P,'tapSup') && P.tapSup.on
    % for beam with tapered support in the middle of beam
    cavlen = 2*P.tapSup.tprlgth + 2*P.tapSup.strghtlgth + ...
             + P.tapSup.extrlgth1 + P.tapSup.extrlgth2 + P.tapSup.centLength;
%     display(num2str(cavlen));
    if (cavlen < 50e-9)
        display(num2str(cavlen));
        error('Cavity length too small to fabricate');
    end
    
    P.cavlen = cavlen;
    
    % half-beam
    xposHalf = cumsum([h_hole(1)/2+cavlen/2,a_hole(2:end)]);
    yposHalf = zeros(size(a_hole));
    h_holeHalf = h_hole;
    d_holeHalf = d_hole;
    w_holeHalf = w_hole;
    a_holeHalf = a_hole;
    
    % full beam
    xposFull = [-1*fliplr(xposHalf),xposHalf];
    yposFull = [yposHalf,yposHalf];
    h_holeFull = [fliplr(h_hole) h_hole];
    d_holeFull = [fliplr(d_hole) d_hole];
    w_holeFull = [fliplr(w_hole) w_hole];
    a_holeFull = [fliplr(a_hole) a_hole];
    
else
    % for beam with dielectric between two holes in the middle of beam
    
    % half-beam
    xposHalf = cumsum([a_hole(1)/2,a_hole(2:end)]); % cumulative sum 
    yposHalf = zeros(size(a_hole)); % placing the nano beam at the center
    h_holeHalf = h_hole;
    d_holeHalf = d_hole;
    w_holeHalf = w_hole;
    a_holeHalf = a_hole;
    
    % full beam
    xposFull = [-1*fliplr(xposHalf),xposHalf];
    yposFull = [yposHalf,yposHalf];
    h_holeFull = [fliplr(h_hole) h_hole];
    d_holeFull = [fliplr(d_hole) d_hole];
    w_holeFull = [fliplr(w_hole) w_hole];
    a_holeFull = [fliplr(a_hole) a_hole];

end

% final geometry parameters
P.a_hole = a_hole';
P.h_hole = h_hole';
P.d_hole = d_hole';
P.w_hole = a_hole';

P.beamLenHalf = xposHalf(nholes+nwgm)+a_hole(nholes)/2;
P.beamLen = xposFull(2*nwgm+nPhM+2*nholes-P.holeatctr)-xposFull(nPhM+1)+a_hole(nholes); % without phon mirror
P.beamLenPhMHalf = xposHalf(end)+a_hole(end)/2;
P.beamLenPhM = xposFull(end)-xposFull(1)+a_hole(end); % with phon mirror
P.geomHalf = [h_holeHalf' d_holeHalf' xposHalf' yposHalf' w_holeHalf' a_holeHalf'];
P.geom = [h_holeFull' d_holeFull' xposFull' yposFull' w_holeFull' a_holeFull'];
end