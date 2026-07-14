function bands = omc_mechanical_livelink(a, w, t, hx, hy, nk, nbands)
% OMC_MECHANICAL_LIVELINK  Mechanical band structure of the 1D diamond OMC cell
% via COMSOL LiveLink for MATLAB. Run with COMSOL 6.2 server + MATLAB connected:
%   comsol mphserver matlab            % (or "COMSOL with MATLAB" desktop)
% then in MATLAB:
%   bands = omc_mechanical_livelink(500e-9,700e-9,220e-9,150e-9,250e-9,15,12);
%
% Returns bands: [nk x nbands] eigenfrequencies (Hz) along Gamma->X.
% Assumes a parametric template built per comsol/README_template.md, with
% parameters a,w,t,hx,hy,kF and an Eigenfrequency study tag 'std1'.

import com.comsol.model.*
import com.comsol.model.util.*

if nargin < 6, nk = 15; end
if nargin < 7, nbands = 12; end

model = mphopen('omc_unitcell.mph');     % template in current dir

% set geometry parameters (SI meters)
model.param.set('a',  sprintf('%g[m]', a));
model.param.set('w',  sprintf('%g[m]', w));
model.param.set('t',  sprintf('%g[m]', t));
model.param.set('hx', sprintf('%g[m]', hx));
model.param.set('hy', sprintf('%g[m]', hy));

svals = linspace(1e-3, 1.0, nk);         % Gamma -> X, s = k_z*a/pi in (0,1]
bands = nan(nk, nbands);

for i = 1:nk
    kF = pi*svals(i)/a;                   % Floquet wavevector along z [1/m]
    model.param.set('kF', sprintf('%g[1/m]', kF));
    model.study('std1').run();            % eigenfrequency study
    f = mphglobal(model, 'freq', 'dataset', 'dset1');  % Hz
    f = sort(real(f));
    n = min(numel(f), nbands);
    bands(i, 1:n) = f(1:n).';
    fprintf('k-point %2d/%2d  s=%.3f  f1=%.3f GHz\n', i, nk, svals(i), f(1)/1e9);
end

% quick gap report near 8 GHz
gapReport(bands, 8e9);
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
