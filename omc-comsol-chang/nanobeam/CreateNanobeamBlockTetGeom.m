% Function to create arrays of geometry parameters for nanobeam with holes
% Modified from createNanobeamCavity from FDTD scripts by Michael Burek
% Input:
% P: data structure consisting of following geometry parameters
%   a: lattice constant;  w: beam width; th: beam height;
%	hx: hole height; hy: hole width; 
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

function [P] = CreateNanobeamBlockTetGeom(P)

% constants
a = P.a;
w = P.w;
hx = P.hx;
hy = P.hy;
nholes = P.nholes;
ndef = P.ndef;
% maxdef = P.maxdef;
% oblong = P.oblong;
holeatctr = P.holeatctr;    % here, holeatctr = tether in ctr, i.e. even total no. of holes
% taperFunc = P.taperFunc;

% % cavity taper function
% if strcmp(taperFunc,'quadratic')
%     f = @(x) 1-x.^2;
% elseif strcmp(taperFunc,'cubic')
%     f = @(x) 1-3*x.^2+2*x.^3;
% elseif strcmp(taperFunc,'linear')
%     f = @(x) 1-x;
% elseif strcmp(taperFunc,'none')
%     f = @(x) 1;
% else
%     error('Invalid taper function specified')
% end
% 
% % end waveguide mirror taper function
% if isfield(P,'wgmTaper') && isfield(P.wgmTaper,'func')
%     if strcmp(P.wgmTaper.func,'quadratic')
%         tf = @(x) 1-x.^2;
%     elseif strcmp(P.wgmTaper.func,'cubic')
%         tf = @(x) 1-3*x.^2+2*x.^3;
%     elseif strcmp(P.wgmTaper.func,'linear')
%         tf = @(x) 1-x;
%     elseif strcmp(P.wgmTaper.func,'none')
%         tf = @(x) 1;
%     else
%         error('Invalid waveguide mirror end taper function specified')
%     end
% end

% check if geometry is within fabrication tolerance
if hx >= a
    error('block dimensions larger than unit cell dimensions')
end
if hy >= w
    error('tether width larger than block width')
end

if isfield(P, 'taperTo')
    if strcmp(P.taperTo,'custom')
        aT = P.a_ctr;
        hxT = P.hx_ctr;
        hyT = P.hy_ctr;
        wT = P.w_ctr;
    else
        aT = 0;
        hxT = 0;
        wT = 0;
        hyT = 0;
    end
else
    aT = 0;
    hxT = 0;
    hyT = 0;
    wT = 0;
end

%% generate list of hx/hy/x/y
for ki = 1:nholes
    if ki <= ndef
        a_hole(ki) = aT;
        w_hole(ki) = wT;
        hx_hole(ki) = hxT;
        hy_hole(ki) = hyT;
    else
        a_hole(ki) = a;
        w_hole(ki) = w;
        hx_hole(ki) = hx;
        hy_hole(ki) = hy;
    end
end

% geometry for half beam
xposHalf = 0.5*([holeatctr*a_hole(1),cumsum(a_hole(1:end-1))+cumsum(a_hole(2:end))]);
yposHalf = 0*zeros(size(xposHalf));
hx_holeHalf = hx_hole;
hy_holeHalf = hy_hole;
w_holeHalf = w_hole;
a_holeHalf = a_hole;

% geometry for full beam
xposFull = [-1*fliplr(xposHalf),xposHalf(2-holeatctr:end)];
yposFull = [fliplr(yposHalf),xposHalf(2-holeatctr:end)];
hx_holeFull = [fliplr(hx_holeHalf),hx_holeHalf(2-holeatctr:end)];
hy_holeFull = [fliplr(hy_holeHalf),hy_holeHalf(2-holeatctr:end)];
w_holeFull = [fliplr(w_holeHalf),w_holeHalf(2-holeatctr:end)];
a_holeFull = [fliplr(a_holeHalf),a_holeHalf(2-holeatctr:end)];

%% final geometry parameters
P.a_hole = a_hole';
P.w_hole = w_hole';
P.hx_hole = hx_hole';
P.hy_hole = hy_hole';
P.beamLenHalf = xposHalf(nholes)+a_hole(nholes)/2;
P.beamLen = 2*xposHalf(nholes)+a_hole(nholes);
P.geomHalf = [hx_holeHalf' hy_holeHalf' xposHalf' yposHalf' w_holeHalf' a_holeHalf'];
P.geom = [hx_holeFull' hy_holeFull' xposFull' yposFull' w_holeFull' a_holeFull'];
end