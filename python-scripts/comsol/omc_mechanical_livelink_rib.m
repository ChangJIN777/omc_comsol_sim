function [bands, bandsAsym] = omc_mechanical_livelink_rib(a, w, th, t, s, r1, r2, nk, nbands)
% OMC_MECHANICAL_LIVELINK_RIB  Mechanical band structure of the 1D diamond rib
% OMC cell (spine + centred rib, the 'rib' cell of buildRibUnitCell_LN.m) via
% COMSOL LiveLink for MATLAB. Rib counterpart of omc_mechanical_livelink.m.
% Run with a COMSOL server + MATLAB connected:
%   comsol mphserver matlab            % (or "COMSOL with MATLAB" desktop)
% then in MATLAB (dimensions of test_rib_unitCell_explicitW.m):
%   [bands, bandsAsym] = omc_mechanical_livelink_rib(570e-9,1747e-9,283e-9, ...
%                            391e-9,254e-9,20e-9,20e-9,15,20);
%
% Inputs (SI metres), same names as the template parameters:
%   a   period along z (cell width)     t   rib width along z (t < a)
%   w   rib tip-to-tip width (y), = 2L  s   spine width (y), full (s < w)
%   th  slab thickness (x)              r1  root fillet, > 0 (t + 2*r1 < a)
%                                       r2  tip fillet, > 0  (2*r2 < t)
%       and r1 + r2 < (w - s)/2 (straight rib flank)
%   nk      number of k-points, default 15
%   nbands  number of bands kept per symmetry, default 20 (the template
%           solves 20 per mechanical study)
%
% Returns, along Gamma->X, [nk x nbands] eigenfrequencies (Hz):
%   bands      y-SYMMETRIC modes   (study std1 'mech sym')
%   bandsAsym  y-ANTISYMMETRIC modes (study std3 'mech antisym')
% Both symmetry studies of the template are solved at every k-point, as in
% omc_unitcell.mph. Each study returns the nbands eigenfrequencies nearest
% its 8 GHz shift.
%
% Assumes the template built by build_omc_unitcell_rib.m, with parameters
% a,w,th,t,s,r1,r2,kF and Eigenfrequency studies 'std1' and 'std3'.
%
% See also BUILD_OMC_UNITCELL_RIB, OMC_MECHANICAL_LIVELINK.

import com.comsol.model.*
import com.comsol.model.util.*

if nargin < 8, nk = 15; end
if nargin < 9, nbands = 20; end

checkRibGeometry(a, w, th, t, s, r1, r2);

model = mphopen('omc_unitcell_rib.mph');     % template in current dir

% Dataset of each study, looked up rather than assumed.
dsetSym  = datasetOfStudy(model, 'std1');
dsetAsym = datasetOfStudy(model, 'std3');

% set geometry parameters (SI meters)
model.param.set('a',  sprintf('%g[m]', a));
model.param.set('w',  sprintf('%g[m]', w));
model.param.set('th', sprintf('%g[m]', th));
model.param.set('t',  sprintf('%g[m]', t));
model.param.set('s',  sprintf('%g[m]', s));
model.param.set('r1', sprintf('%g[m]', r1));
model.param.set('r2', sprintf('%g[m]', r2));

svals = linspace(1e-3, 1.0, nk);         % Gamma -> X, s = k_z*a/pi in (0,1]
bands     = nan(nk, nbands);
bandsAsym = nan(nk, nbands);

for i = 1:nk
    kF = pi*svals(i)/a;                   % Floquet wavevector along z [1/m]
    model.param.set('kF', sprintf('%g[1/m]', kF));

    model.study('std1').run();            % y-symmetric eigenfrequency study
    f = mphglobal(model, 'freq', 'dataset', dsetSym);  % Hz
    f = sort(real(f));
    n = min(numel(f), nbands);
    bands(i, 1:n) = f(1:n).';

    model.study('std3').run();            % y-antisymmetric eigenfrequency study
    g = mphglobal(model, 'freq', 'dataset', dsetAsym); % Hz
    g = sort(real(g));
    m = min(numel(g), nbands);
    bandsAsym(i, 1:m) = g(1:m).';

    fprintf('k-point %2d/%2d  s=%.3f  sym f1=%.3f GHz  antisym f1=%.3f GHz\n', ...
        i, nk, svals(i), f(1)/1e9, g(1)/1e9);
end

% quick gap reports near 8 GHz
fprintf('\n[y-symmetric]');
gapReport(bands, 8e9);
fprintf('[y-antisymmetric]');
gapReport(bandsAsym, 8e9);
fprintf('[both symmetries]');
gapReportUnion({bands, bandsAsym}, 8e9);
end

function dset = datasetOfStudy(model, studyTag)
% Tag of the solution dataset that belongs to study studyTag.
tags = cell(model.result.dataset.tags);
for i = 1:numel(tags)
    try
        solTag = char(model.result.dataset(tags{i}).getString('solution'));
        if strcmp(char(model.sol(solTag).study), studyTag)
            dset = tags{i};
            return;
        end
    catch
        % not a solution dataset (cut plane, etc.)
    end
end
error('omc_mechanical_livelink_rib:dataset', ...
    'No dataset belongs to study %s in omc_unitcell_rib.mph.', studyTag);
end

function checkRibGeometry(a, w, th, t, s, r1, r2)
% Reject geometries the template cannot build, before any solve. The same
% rules as omc-comsol-chang/validateRibParams.m for d = 0, except that
% r1, r2 = 0 is rejected (COMSOL's Fillet errors on radius 0) and the
% straight-flank rule (a warning there) is enforced.
dims = [a w th t s r1 r2];
if ~isnumeric(dims) || numel(dims) ~= 7 || any(~isfinite(dims)) || any(dims <= 0)
    error('omc_mechanical_livelink_rib:badInput', ...
        'a, w, th, t, s, r1, r2 must be positive finite scalars in metres.');
end
if t >= a
    error('omc_mechanical_livelink_rib:ribTooWide', ...
        'Need t < a (got t=%g, a=%g); otherwise the ribs merge.', t, a);
end
if s >= w
    error('omc_mechanical_livelink_rib:spineTooWide', ...
        'Need s < w (got s=%g, w=%g); otherwise there is no rib.', s, w);
end
if 2*r2 >= t
    error('omc_mechanical_livelink_rib:tipFilletTooLarge', ...
        'Need 2*r2 < t (got r2=%g, t=%g).', r2, t);
end
if t + 2*r1 >= a
    error('omc_mechanical_livelink_rib:rootFilletTooLarge', ...
        'Need t + 2*r1 < a (got %g >= %g).', t + 2*r1, a);
end
if r1 + r2 >= (w - s)/2
    error('omc_mechanical_livelink_rib:filletsMeet', ...
        ['Need r1 + r2 < (w - s)/2 (got %g >= %g); otherwise the tip and ' ...
         'root fillets meet on the rib flank.'], r1 + r2, (w - s)/2);
end
end

function gapReport(bands, ftarget)
b = sort(bands, 2);
best = 0; lo = 0; hi = 0;
for j = 1:size(b,2)-1
    topLower = max(b(:,j));
    botUpper = min(b(:,j+1));
    if botUpper > topLower
        c = 0.5*(topLower+botUpper);
        ng = (botUpper-topLower)/c;
        if abs(c-ftarget)/ftarget < 0.5 && ng > best
            best = ng; lo = topLower; hi = botUpper;
        end
    end
end
fprintf('\nLargest gap near %.1f GHz: %.1f%%  (%.3f-%.3f GHz)\n', ...
        ftarget/1e9, 100*best, lo/1e9, hi/1e9);
end

function gapReportUnion(sets, ftarget)
% Complete gap over several symmetry sectors, by frequency interval: the
% widest interval with no eigenfrequency of ANY sector at ANY k-point.
% Band indices are not used, because each study returns the modes nearest
% its shift, so index j need not be the same band in two sectors. Only the
% frequency window every sector covers at every k is trusted: outside it a
% sector's modes were not computed.
fLo = -Inf;  fHi = Inf;  F = [];
for s = 1:numel(sets)
    B = sets{s};
    fLo = max(fLo, max(min(B, [], 2)));
    fHi = min(fHi, min(max(B, [], 2)));
    F = [F; B(isfinite(B))]; %#ok<AGROW>
end
F = unique(F(F >= fLo & F <= fHi));
best = 0; lo = 0; hi = 0;
for j = 1:numel(F)-1
    c = 0.5*(F(j) + F(j+1));
    ng = (F(j+1) - F(j))/c;
    if abs(c-ftarget)/ftarget < 0.5 && ng > best
        best = ng; lo = F(j); hi = F(j+1);
    end
end
fprintf('\nLargest gap near %.1f GHz: %.1f%%  (%.3f-%.3f GHz)  [trusted window %.3f-%.3f GHz]\n', ...
        ftarget/1e9, 100*best, lo/1e9, hi/1e9, fLo/1e9, fHi/1e9);
end
