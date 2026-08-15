function [lightline,use2Dpath] = opticalLightLine(OpticalBand,P)
%OPTICALLIGHTLINE Light line in Hz at every k-point of an optical band sweep.
%
% Returns the SAME curve solveOpticalBands draws on the band diagram, sampled
% at whatever k-points it is handed rather than on a fixed display grid, so the
% modes that pass a below-light-line test are exactly the modes a reader sees
% underneath the blue curves. Both solveOpticalBands and test_plotOpticalBands
% call this, which is the point: the filter and the figure must not drift apart.
%
% INPUTS
%   OpticalBand  ds struct from runOpticalBand_{1,2,3}D, or any struct with the
%                same fields. Needs k_norm on the 2D circuit; needs kx_norm and
%                ky_norm on the 1D path.
%   P            parameter struct. Needs a, bandStructureDim, and - for the 2D
%                circuit - bandStruct_2D.
%
% OUTPUTS
%   lightline    [nk x 1] light line in Hz, one entry per supplied k-point.
%   use2Dpath    true when the sweep is the 0->3 Gamma-X-M-Gamma circuit.
%
% On the 2D circuit the curve is assembled from three segments in the path
% parameter q = k_norm:
%
%   Gamma-X  (0<=q<1):  t = q/2,      f = t
%   X-M      (1<=q<2):  t = (q-1)/2,  f = sqrt(1/4 + (t/sqrt(3))^2)
%   M-Gamma  (2<=q<=3): t = (q-2)/2,  f = sqrt((1/2-t)^2 + (1/(2*sqrt(3))-t/sqrt(3))^2)
%
% all scaled by c/a. The three agree at the joins - 0.5 at X, 0.5774 at M, 0 at
% the closing Gamma - so the assembled curve is continuous.
%
% CAVEAT worth knowing before quoting a gap as physics: this is NOT the same as
% c*hypot(kx,ky)/(2a) built from COMSOL's own kx/ky parameter expressions. The
% two agree along Gamma-X and at X but diverge along X-M, differing by
% 2/sqrt(3) = 1.1547 at M, because they disagree about where the high-symmetry
% points sit. This function reproduces the DRAWN curve by design.
%
% Example:
%   S  = load('optical_boomerang_..._bds.mat');
%   ll = opticalLightLine(S.ds.opticalBand, S.ds.opticalBand.P);

narginchk(2,2);

c = 299792458;

use2Dpath = ismember(P.bandStructureDim,[2 3]) && ...
            isfield(P,'bandStruct_2D') && P.bandStruct_2D;

if use2Dpath
    q         = OpticalBand.k_norm(:);
    lightline = zeros(size(q));

    seg1 = q < 1;
    seg2 = q >= 1 & q < 2;
    seg3 = q >= 2;

    t1 = q(seg1)/2;
    t2 = (q(seg2)-1)/2;
    t3 = (q(seg3)-2)/2;

    lightline(seg1) = t1;
    lightline(seg2) = sqrt(1/4 + (t2./sqrt(3)).^2);
    lightline(seg3) = sqrt((1/2-t3).^2 + (1/(2*sqrt(3))-t3./sqrt(3)).^2);

    lightline = lightline*c/P.a;
else
    % 1D path: k_norm IS kx_norm and ky_norm is zero throughout, so the
    % magnitude form is exactly what the 1D plot branch draws.
    knorm     = hypot(OpticalBand.kx_norm(:), OpticalBand.ky_norm(:));
    lightline = c*knorm/2/P.a;   % /2 so knorm runs 0 to 0.5*(2*pi/a)
end
end
