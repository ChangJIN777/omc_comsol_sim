function [bands, bandsAsym] = omc_mechanical_livelink_rib(a, w, th, wi, wo, hi, ho, ai, nk, nbands)
% OMC_MECHANICAL_LIVELINK_RIB  Mechanical band structure of the 1D diamond rib
% OMC cell via COMSOL LiveLink for MATLAB. Rib counterpart of
% omc_mechanical_livelink.m. Run with a COMSOL server + MATLAB connected:
%   comsol mphserver matlab            % (or "COMSOL with MATLAB" desktop)
% then in MATLAB:
%   [bands, bandsAsym] = omc_mechanical_livelink_rib(1500e-9,3000e-9,160e-9, ...
%                            670e-9,1000e-9,955e-9,975e-9,150e-9,15,20);
%
% Inputs (SI metres), same names as the template parameters:
%   a   period along z                 wi  tooth length along z (wi < wo)
%   w   full beam width (y)            wo  void length along z  (wo < a)
%   th  beam thickness (x)             hi  tooth height along y (hi < ho)
%   ai  spine width (void starts       ho  void height along y
%       at y = ai/2)                       (ai/2 + ho < w/2)
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
% a,w,th,wi,wo,hi,ho,ai,kF and Eigenfrequency studies 'std1' and 'std3'.
%
% See also BUILD_OMC_UNITCELL_RIB, OMC_MECHANICAL_LIVELINK.

import com.comsol.model.*
import com.comsol.model.util.*

if nargin < 9,  nk = 15; end
if nargin < 10, nbands = 20; end

checkRibGeometry(a, w, th, wi, wo, hi, ho, ai);

model = mphopen('omc_unitcell_rib.mph');     % template in current dir

% Dataset of each study, looked up rather than assumed.
dsetSym  = datasetOfStudy(model, 'std1');
dsetAsym = datasetOfStudy(model, 'std3');

% set geometry parameters (SI meters)
model.param.set('a',  sprintf('%g[m]', a));
model.param.set('w',  sprintf('%g[m]', w));
model.param.set('th', sprintf('%g[m]', th));
model.param.set('wi', sprintf('%g[m]', wi));
model.param.set('wo', sprintf('%g[m]', wo));
model.param.set('hi', sprintf('%g[m]', hi));
model.param.set('ho', sprintf('%g[m]', ho));
model.param.set('ai', sprintf('%g[m]', ai));

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

function checkRibGeometry(a, w, th, wi, wo, hi, ho, ai)
% Reject geometries the template cannot build, before any solve.
dims = [a w th wi wo hi ho ai];
if ~isnumeric(dims) || numel(dims) ~= 8 || any(~isfinite(dims)) || any(dims <= 0)
    error('omc_mechanical_livelink_rib:badInput', ...
        'a, w, th, wi, wo, hi, ho, ai must be positive finite scalars in metres.');
end
if wi >= wo || hi >= ho
    error('omc_mechanical_livelink_rib:badVoid', ...
        ['Need wi < wo and hi < ho (got wi=%g, wo=%g, hi=%g, ho=%g); ' ...
         'otherwise the void polygon is degenerate or self-intersecting.'], ...
        wi, wo, hi, ho);
end
if wo >= a
    error('omc_mechanical_livelink_rib:voidTooLong', ...
        'Need wo < a (got wo=%g, a=%g); otherwise the void cuts the whole period.', ...
        wo, a);
end
if ai/2 + ho >= w/2
    error('omc_mechanical_livelink_rib:voidTooHigh', ...
        ['Need ai/2 + ho < w/2 (got %g >= %g); otherwise the void breaks ' ...
         'through the beam edge.'], ai/2 + ho, w/2);
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
