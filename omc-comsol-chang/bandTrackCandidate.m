function [cand, S] = bandTrackCandidate(sec)
%BANDTRACKCANDIDATE  Turn one solved symmetry sector into tracker input.
%
%   [cand, S] = bandTrackCandidate(sec)
%
%   sec    a sector struct from runBands / solveBands (e.g. ds.sym) solved
%          with P.exportModes, so that sec.modes exists and includes the X
%          point (the last row of sec.F), where bands are identified.
%
%   cand   modes at X, sorted by frequency: f (1 x nb), U (3*nPts x nb), w,
%          fp (features x nb), fpNames, px (x-mirror parity), order (sorted
%          column -> column of sec.F)
%   S      connectBandsAlongK output. Mode-shape connection is used when the
%          fields were exported at every k, frequency continuity otherwise.

if isempty(sec) || ~isstruct(sec) || ~isfield(sec,'modes') || isempty(sec.modes)
    error('bandTrackCandidate:noModes', ...
        ['The sector has no sampled modes. Solve with P.exportModes set; a ' ...
         'cached solveBands call returns an empty sector (load the _bds.mat).']);
end
md = sec.modes;
nk = size(sec.F, 1);
m = find(md.kIdx == nk, 1);
if isempty(m)
    error('bandTrackCandidate:anchorNotExported', ...
        'The X point (row %d of F) is not among the exported k-points.', nk);
end

[cand.f, ord] = sort(sec.F(nk,:));
cand.order = ord;
cand.U = double(md.U(:,ord,m));
cand.w = md.grid.w;
dv = [];
if ~isempty(md.divU), dv = md.divU(:,ord,m); end
FP = modeFingerprints(cand.U, md.grid, dv);
cand.fp = FP.vec;  cand.fpNames = FP.names;  cand.px = FP.px;

Uall = [];
if isequal(md.kIdx(:)', 1:nk), Uall = md.U; end
S = connectBandsAlongK(sec.F, Uall, md.grid.w, nk);
end
