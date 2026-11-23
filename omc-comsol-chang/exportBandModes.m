function modes = exportBandModes(model, P, ds)
%EXPORTBANDMODES  Sample the displacement of every band at chosen k-points.
%
%   modes = exportBandModes(model, P, ds)
%
%   Called by runBands when P.exportModes is set, right after the solve and
%   while the model is still alive. One mphinterp call per exported k-point.
%
%   P.exportModes   true, or a struct with any of
%       k        'X' (default) | 'GX' | 'all' | row indices into ds.F
%       n        point counts for ribModeGrid
%       divU     also sample the volumetric strain (default true)
%       dataset  parametric dataset tag (default 'pdset')
%
%   modes.U     3*nPts x nbands x nK complex single, rows [ux; uy; uz] in the
%               order of modes.grid; column j is band column j of ds.F at
%               that k (solver order, NOT sorted). NaN outside the solid.
%   modes.divU  nPts x nbands x nK ([] if not sampled)
%   modes.grid  ribModeGrid output
%   modes.kIdx, modes.k_norm, modes.F   the exported rows of ds.F
%   modes.a, modes.sector ([mbeveny mbevenz])
%   modes.freqCheck   max relative mismatch between ds.F and COMSOL's 'freq'
%               for the same (outersolnum, solnum); NaN if not evaluated.
%               Nonzero means the solnum order differs from ds.F's columns.
%
%   Uses only the dependent variables u, v, w (no physics-scoped names), so it
%   does not depend on the Solid Mechanics interface identifier.

opt = struct('k','X', 'n',struct(), 'divU',true, 'dataset','pdset');
if isstruct(P.exportModes)
    fn = fieldnames(P.exportModes);
    for k = 1:numel(fn), opt.(fn{k}) = P.exportModes.(fn{k}); end
end
[nk, nb] = size(ds.F);
if ischar(opt.k) || isstring(opt.k)
    switch lower(char(opt.k))
        case 'x',   kIdx = nk;
        case 'gx',  kIdx = [1 nk];
        case 'all', kIdx = 1:nk;
        otherwise
            error('exportBandModes:badK', ...
                'P.exportModes.k must be ''X'', ''GX'', ''all'' or row indices.');
    end
else
    kIdx = reshape(double(opt.k), 1, []);
end

G = ribModeGrid(P, opt.n);
nPts = size(G.xyz, 2);
exprs = {'u','v','w'};
if opt.divU, exprs{end+1} = 'd(u,X)+d(v,Y)+d(w,Z)'; end

nK = numel(kIdx);
U = complex(nan(3*nPts, nb, nK, 'single'));
divU = [];
if opt.divU, divU = complex(nan(nPts, nb, nK, 'single')); end
fChk = nan(1, nK);
for m = 1:nK
    ki = kIdx(m);
    out = cell(1, numel(exprs));
    [out{:}] = mphinterp(model, exprs, 'coord', G.xyz, 'edim', 3, ...
        'dataset', opt.dataset, 'outersolnum', ki, 'solnum', 1:nb);
    for c = 1:3
        U((c-1)*nPts + (1:nPts), :, m) = single(orient(out{c}, nb, nPts));
    end
    if opt.divU
        divU(:,:,m) = single(orient(out{4}, nb, nPts));
    end
    try
        fq = mphglobal(model, 'freq', 'dataset', opt.dataset, ...
            'outersolnum', ki, 'solnum', 1:nb);
        fChk(m) = max(abs(abs(fq(:)') - ds.F(ki,:))./max(ds.F(ki,:), 1));
    catch
    end
end
if any(fChk > 1e-6)
    warning('exportBandModes:orderMismatch', ...
        ['The sampled modes do not line up with the columns of ds.F (max ' ...
         'relative frequency mismatch %.2g).'], max(fChk));
end

modes = struct('U',U, 'divU',divU, 'grid',G, 'kIdx',kIdx, ...
    'k_norm',reshape(ds.k_norm(kIdx),1,[]), 'F',ds.F(kIdx,:), 'a',P.a, ...
    'sector',[P.mbeveny P.mbevenz], 'freqCheck',fChk);
end

% -------------------------------------------------------------------------

function A = orient(A, nb, nPts)
% mphinterp returns nSol x nPts; store as nPts x nSol. (.' : no conjugation)
if isequal(size(A), [nb nPts])
    A = A.';
elseif ~isequal(size(A), [nPts nb])
    error('exportBandModes:badSize', ...
        'mphinterp returned %d x %d, expected %d x %d.', size(A,1), size(A,2), nb, nPts);
end
end
