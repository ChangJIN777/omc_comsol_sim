function test_isFabricable()
%TEST_ISFABRICABLE  Cross-validate the pre-filter against the real geometry builder.
%
%   isFabricable exists to reject unfabricable taper designs without paying for
%   a solve. That is only sound if it agrees with CreateNanobeamGeom, which is
%   the code that actually throws. CreateNanobeamGeom is pure arithmetic on the
%   parameter arrays -- no COMSOL -- so the two can be compared directly.
%
%   The asymmetry matters:
%     FALSE POSITIVE  (we say fabricable, the builder throws) is UNSAFE -- it
%                     buys a failed multi-minute solve, exactly what the filter
%                     is supposed to prevent. Must be zero.
%     FALSE NEGATIVE  (we prune a design the builder accepts) is merely
%                     conservative, but it silently shrinks the search space,
%                     so it is asserted to be zero too.
%
%   Usage:  test_isFabricable

fprintf('=== isFabricable cross-validation ===\n\n');

P = struct('a',529e-9, 'w',750e-9, 'hx',196e-9, 'hy',578e-9, 'th',500e-9, ...
    'nholes',18, 'ndef',8, 'taperFunc','cubic', 'holeatctr',1, ...
    'useManualDefect',0, 'wvgmir',0, 'xsect','isoFit', ...
    'maxdef',0.22, 'oblong',1.96);

dAR0 = (P.hy/P.hx)*(1-0.22)^(2*1.96);

%% ---------- 1. the nominal design must be fabricable ----------
[ok0, m0, why0] = isFabricable(dAR0, 0.22, P);
fprintf('nominal design: ok=%d, margin=%.1f nm\n', ok0, m0*1e9);
assert(ok0, 'the nominal fabricated design was rejected: %s', why0{1});

%% ---------- 2. agreement over the current search box ----------
runScan('search box', P, ...
    linspace(dAR0*0.5, dAR0*1.5, 21), linspace(0.11, 0.33, 21));

%% ---------- 3. agreement over a WIDER region ----------
% The box scan above is expected to be entirely fabricable, which would make
% the comparison vacuous. This region is chosen to straddle the boundary so
% the filter has to actually discriminate.
fracRejected = runScan('wide region', P, ...
    linspace(dAR0*0.25, dAR0*2.5, 25), linspace(0.05, 0.60, 25));
assert(fracRejected > 0.02, ...
    ['the wide scan rejected almost nothing (%.1f%%), so it does not ' ...
     'actually test the discrimination'], 100*fracRejected);

%% ---------- 4. vectorization contract used by XConstraintFcn ----------
% bayesopt hands XConstraintFcn a table of many candidate rows and expects a
% logical column back.
dv = [0.9; 1.1; 1.4];
mv = [0.15; 0.22; 0.28];
okv = isFabricable(dv, mv, P);
assert(iscolumn(okv) && numel(okv) == 3 && islogical(okv), ...
    'vectorized call must return a logical column vector');
for k = 1:3
    assert(okv(k) == isFabricable(dv(k), mv(k), P), ...
        'vectorized result disagrees with the scalar call at row %d', k);
end
fprintf('vectorization: logical column, matches scalar calls\n');

%% ---------- 5. scalar/array broadcasting ----------
okb = isFabricable(dAR0, [0.15; 0.22; 0.30], P);
assert(numel(okb) == 3, 'scalar dAR should broadcast against an array maxdef');
okb2 = isFabricable([0.9; 1.2], 0.22, P);
assert(numel(okb2) == 2, 'scalar maxdef should broadcast against an array dAR');
fprintf('broadcasting  : scalar/array mixes expand correctly\n');

%% ---------- 6. non-physical parameters are rejected, not crashed ----------
bad = isFabricable([-1; 0.5; 0.5], [0.22; 0; 1], P);
assert(~any(bad), 'non-physical taper parameters must be rejected');
fprintf('robustness    : non-physical inputs rejected without error\n');

fprintf('\nALL ISFABRICABLE ASSERTIONS PASSED\n');
end

% =========================================================================
function fracRejected = runScan(label, P, dARs, mds)
%RUNSCAN  Compare isFabricable with CreateNanobeamGeom over a grid.

[DA, MD] = meshgrid(dARs, mds);
OK = isFabricable(DA, MD, P);

builderOK = false(size(DA));
for k = 1:numel(DA)
    Pk = P;
    Pk.maxdef = MD(k);
    Pk.oblong = log(DA(k)/(P.hy/P.hx)) / (2*log(1 - MD(k)));
    try
        CreateNanobeamGeom(Pk);
        builderOK(k) = true;
    catch
        builderOK(k) = false;
    end
end

falsePos = OK & ~builderOK;
falseNeg = ~OK & builderOK;
fracRejected = mean(~OK(:));

fprintf('%-12s %dx%d : filter %.1f%% ok, builder %.1f%% ok, agreement %.1f%%\n', ...
    label, numel(dARs), numel(mds), ...
    100*mean(OK(:)), 100*mean(builderOK(:)), 100*mean(OK(:) == builderOK(:)));

if any(falsePos(:))
    k = find(falsePos, 1);
    error('test_isFabricable:falsePositive', ...
        ['UNSAFE: isFabricable accepted a design CreateNanobeamGeom rejects ' ...
         '(%d of %d), e.g. dAR=%.4f maxdef=%.4f'], ...
        nnz(falsePos), numel(DA), DA(k), MD(k));
end
if any(falseNeg(:))
    k = find(falseNeg, 1);
    [~, ~, w] = isFabricable(DA(k), MD(k), P);
    error('test_isFabricable:falseNegative', ...
        ['isFabricable pruned a design CreateNanobeamGeom accepts (%d of %d), ' ...
         'e.g. dAR=%.4f maxdef=%.4f rejected by "%s"'], ...
        nnz(falseNeg), numel(DA), DA(k), MD(k), w{1});
end
end
