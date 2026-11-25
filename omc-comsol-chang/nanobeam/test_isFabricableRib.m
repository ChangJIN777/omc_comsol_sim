function test_isFabricableRib()
%TEST_ISFABRICABLERIB  Cross-validate the rib pre-filter against the real builder.
%
%   isFabricableRib rejects unbuildable / unfabricable rib tapers without a
%   solve. That is only sound if its BUILD_* rules agree with
%   createRibCavityGeom (which calls validateRibParams) -- the code that
%   actually throws. Both are pure MATLAB, so they are compared directly.
%
%     FALSE POSITIVE  (filter says valid, the builder throws) is UNSAFE: it
%                     buys a failed evaluation. Must be zero.
%     FALSE NEGATIVE  (filter prunes a design the builder accepts) silently
%                     shrinks the search space. Must be zero too.
%   The WARN_* rules are compared with the builder's warnings the same way.
%   The FAB_* floors are deliberately stricter than the builder, so they are
%   tested separately (section 5).
%
%   No COMSOL needed.  Usage:  test_isFabricableRib

% ribDefaultParams / validateRibParams / ribLength live in the root folder.
thisDir = fileparts(mfilename('fullpath'));
addpath(fileparts(thisDir), '-end');

% The builder warns a lot on the wide grids; lastwarn still records them.
warnState = warning('off', 'all');
restoreWarn = onCleanup(@() warning(warnState));

fprintf('=== isFabricableRib cross-validation ===\n\n');

% The design of test_nanobeamRibFEM.m.
P = struct('celltype','rib', 'xsect','rect', 'a',600e-9, 't',280e-9, ...
    's',100e-9, 'th',250e-9, 'r1',40e-9, 'r2',40e-9, 'w_mirror',1300e-9, ...
    'd',0, 'ncellHalf',14, 'ribatctr',1, 'ndef',5, 'taperFunc','quadratic', ...
    'taperConvention','centerspacing', 'ribTaperScheme','hole', ...
    'airrad',3.1e-6, 'airZOffset',0, 'mEndBC','lowreflecting', ...
    'maxdef',0.2, 'oblong',-2);
mirrorAR = P.w_mirror/P.t;
dAR0 = mirrorAR*(1 - 0.2)^(2*(-2));

%% ---------- 1. the nominal design, and the mapping ----------
[ok0, m0, why0, info0] = isFabricableRib(dAR0, 0.2, P);
fprintf('nominal design: ok=%d, margin=%.1f nm, oblong=%.4f, t_ctr=%.1f nm, w_ctr=%.1f nm\n', ...
    ok0, m0*1e9, info0.oblong, info0.t_ctr*1e9, info0.w_ctr*1e9);
assert(ok0, 'the nominal design was rejected: %s', why0{1});
assert(abs(info0.oblong - (-2)) < 1e-12, 'aspect ratio -> oblong mapping is off');
G = createRibCavityGeom(P);
assert(abs(info0.t_ctr - G.t_cell(1)) < 1e-15 && ...
       abs(info0.w_ctr - 2*G.L_cell(1)) < 1e-15 && ...
       abs(info0.a_ctr - G.a_cell(1)) < 1e-15, ...
    'centre-cell dimensions disagree with createRibCavityGeom');

%% ---------- 2. agreement over the optimizer's default search box ----------
runScan('search box', P, ...
    exp(linspace(log(mirrorAR/3), log(dAR0*1.5), 21)), linspace(0.05, 0.30, 21));

%% ---------- 3. agreement over a WIDER region ----------
% The box above is entirely valid, which would make the comparison vacuous.
% This region straddles the boundary on purpose.
fracRejected = runScan('wide region', P, ...
    exp(linspace(log(0.3), log(40), 25)), linspace(0.02, 0.70, 25));
assert(fracRejected > 0.02, ...
    ['the wide scan rejected almost nothing (%.1f%%), so it does not ' ...
     'actually test the discrimination'], 100*fracRejected);

%% ---------- 4. the other geometry conventions ----------
variants = { ...
    'cellwidth',      struct('taperConvention','cellwidth'); ...
    'gap at centre',  struct('ribatctr',0); ...
    'cubic',          struct('taperFunc','cubic'); ...
    'linear, ndef 8', struct('taperFunc','linear', 'ndef',8); ...
    'none, ndef 1',   struct('taperFunc','none', 'ndef',1); ...
    'ndef 0',         struct('ndef',0); ...
    'N 3 < ndef',     struct('ncellHalf',3)};
for v = 1:size(variants, 1)
    Pvar = P;
    fn = fieldnames(variants{v, 2});
    for j = 1:numel(fn)
        Pvar.(fn{j}) = variants{v, 2}.(fn{j});
    end
    runScan(variants{v, 1}, Pvar, ...
        exp(linspace(log(0.3), log(40), 15)), linspace(0.02, 0.70, 15));
end

%% ---------- 5. the fab floors ----------
[DA, MD] = meshgrid(exp(linspace(log(0.3), log(40), 25)), linspace(0.02, 0.70, 25));
Poff = P;
Poff.minRibFeature   = 0;
Poff.minRibGap       = 0;
Poff.minRibThickness = 0;
Poff.ribFabStrict    = 0;
[okOff, ~, ~, infoOff] = isFabricableRib(DA, MD, Poff);
assert(isequal(okOff, infoOff.valid), ...
    'with every WARN/FAB rule off, ok must equal the builder-validity mask');
[okOn, ~, ~, infoOn] = isFabricableRib(DA, MD, P);
assert(all(infoOn.valid(okOn)), 'a design passed the filter but is not buildable');
assert(nnz(okOn) <= nnz(okOff), 'the fab floors must never enlarge the feasible set');
% A 400 nm rib-gap floor: the narrowest gap of the nominal design is the
% mirror cell's a - t = 320 nm, so the shortfall is exactly 80 nm.
Pbig = P;
Pbig.minRibGap = 400e-9;
[okBig, mBig, whyBig] = isFabricableRib(dAR0, 0.2, Pbig);
assert(~okBig && strcmp(whyBig{1}, 'FAB_MIN_GAP'), ...
    'a 400 nm gap floor must reject the nominal design as FAB_MIN_GAP');
assert(abs(mBig - (P.a - P.t - 400e-9)) < 1e-15, ...
    'the fab margin must be the shortfall in metres');
fprintf('fab floors    : off == builder validity; on only shrinks; reason code + margin correct\n');

%% ---------- 6. vectorization contract used by XConstraintFcn ----------
% bayesopt hands XConstraintFcn a table of many candidate rows and expects a
% logical column back. The last row violates the root-fillet rule.
dv = [2; 8; 16; 1.2];
mv = [0.10; 0.20; 0.28; 0.34];
okv = isFabricableRib(dv, mv, P);
assert(iscolumn(okv) && numel(okv) == 4 && islogical(okv), ...
    'vectorized call must return a logical column vector');
assert(isequal(okv, [true; true; true; false]), ...
    'vectorized call returned the wrong feasibility pattern');
for k = 1:4
    assert(okv(k) == isFabricableRib(dv(k), mv(k), P), ...
        'vectorized result disagrees with the scalar call at row %d', k);
end
fprintf('vectorization : logical column, matches scalar calls\n');

%% ---------- 7. scalar/array broadcasting ----------
okb = isFabricableRib(dAR0, [0.10; 0.20; 0.30], P);
assert(numel(okb) == 3, 'scalar aspect ratio should broadcast against an array maxdef');
okb2 = isFabricableRib([4; 12], 0.2, P);
assert(numel(okb2) == 2, 'scalar maxdef should broadcast against an array aspect ratio');
fprintf('broadcasting  : scalar/array mixes expand correctly\n');

%% ---------- 8. non-physical parameters are rejected, not crashed ----------
[bad, mBad, whyBad] = isFabricableRib([-1; 5; 5; NaN], [0.2; 0; 1; 0.2], P);
assert(~any(bad) && all(isinf(mBad) & mBad < 0) && ...
       all(strcmp(whyBad, 'NONPHYSICAL_TAPER')), ...
    'non-physical taper parameters must be rejected with margin -Inf');
fprintf('robustness    : non-physical inputs rejected without error\n');

fprintf('\nALL ISFABRICABLERIB ASSERTIONS PASSED\n');
end

% =========================================================================
function fracRejected = runScan(label, P, dARs, mds)
%RUNSCAN  Compare isFabricableRib with createRibCavityGeom over a grid.

[DA, MD] = meshgrid(dARs, mds);
[~, ~, ~, info] = isFabricableRib(DA, MD, P);

mirrorAR  = P.w_mirror/P.t;
builderOK = false(size(DA));
builderWarn = false(size(DA));
for k = 1:numel(DA)
    Pk = P;
    Pk.maxdef = MD(k);
    Pk.oblong = log(DA(k)/mirrorAR) / (2*log(1 - MD(k)));
    lastwarn('');
    try
        createRibCavityGeom(Pk);
        builderOK(k) = true;
        [~, warnId] = lastwarn;
        builderWarn(k) = ~isempty(warnId);
    catch
        builderOK(k) = false;
    end
end

falsePos = info.valid & ~builderOK;
falseNeg = ~info.valid & builderOK;
warnDiff = builderOK & (info.warn ~= builderWarn);
fracRejected = mean(~info.valid(:));

fprintf('%-14s %dx%d : filter %5.1f%% valid, builder %5.1f%% valid, agreement %5.1f%%, warned %d\n', ...
    label, numel(dARs), numel(mds), 100*mean(info.valid(:)), ...
    100*mean(builderOK(:)), 100*mean(info.valid(:) == builderOK(:)), ...
    nnz(builderOK & builderWarn));

if any(falsePos(:))
    k = find(falsePos, 1);
    error('test_isFabricableRib:falsePositive', ...
        ['UNSAFE: isFabricableRib accepted a design createRibCavityGeom rejects ' ...
         '(%d of %d), e.g. aspect ratio=%.4f maxdef=%.4f'], ...
        nnz(falsePos), numel(DA), DA(k), MD(k));
end
if any(falseNeg(:))
    k = find(falseNeg, 1);
    error('test_isFabricableRib:falseNegative', ...
        ['isFabricableRib pruned a design createRibCavityGeom accepts (%d of %d), ' ...
         'e.g. aspect ratio=%.4f maxdef=%.4f'], ...
        nnz(falseNeg), numel(DA), DA(k), MD(k));
end
if any(warnDiff(:))
    k = find(warnDiff, 1);
    error('test_isFabricableRib:warnMismatch', ...
        ['WARN_* rules disagree with the builder''s warnings (%d of %d), ' ...
         'e.g. aspect ratio=%.4f maxdef=%.4f: filter %d, builder %d'], ...
        nnz(warnDiff), numel(DA), DA(k), MD(k), info.warn(k), builderWarn(k));
end
end
