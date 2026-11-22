function A = unitCellAdapter_hole()
%UNITCELLADAPTER_HOLE  'hole' celltype adapter for bayesoptUnitCell.
%
%   A = unitCellAdapter_hole();  cfg = A.defaultCfg();
%   ok = A.fabOKmetres(a, hx, hy, bw, th, cfg)    % metres in
%
%   Geometry (buildHoleUnitCell.m): beam of width beam_width (bw) and
%   thickness th, period a, one elliptic hole with DIAMETERS hx (along the
%   beam) and hy (across it). P defaults follow optimize_hole_unitCell.m:59-98.
%
%   See also BAYESOPTUNITCELL, BAYESOPT_HOLE, OPTIMIZE_HOLE_UNITCELL.
A.name     = 'hole';
A.varNames = {'a', 'hx', 'hy', 'bw', 'th'};
A.show = struct( ...
    'a',  struct('scale', 1e9, 'label', 'a_nm',         'unit', 'nm'), ...
    'hx', struct('scale', 1e9, 'label', 'hx_nm',        'unit', 'nm'), ...
    'hy', struct('scale', 1e9, 'label', 'hy_nm',        'unit', 'nm'), ...
    'bw', struct('scale', 1e9, 'label', 'beamWidth_nm', 'unit', 'nm'), ...
    'th', struct('scale', 1e9, 'label', 'th_nm',        'unit', 'nm'));
A.defaultCfg  = @holeDefaultCfg;
A.prepareCfg  = @holePrepareCfg;
A.derive      = @holeDerive;
A.fabOK       = @holeFabOKD;
A.fabOKmetres = @holeFabOK;
A.buildP      = @holeBuildP;
A.surrogate   = @holeSurrogate;
end

%% ========================================================================
function cfg = holeDefaultCfg()
cfg = bayesoptUnitCell('defaults');
nm  = 1e-9;

% Fixed values: optimize_hole_unitCell.m (a :65, hx0/hy0 :19-20, beam_width :68,
% th :71). TODO(confirm) the baseline: sweep_hole_unitCell.m uses a = 528.6,
% bw = 800, th = 500, hx0 = 195.6, hy0 = 651 nm (:24-25, :68, :71, :74) and
% test_hole_unitCell.m a = 528, hx = 197, hy = 578, bw = 750, th = 500 (:14-25).
% Bounds: NO existing script has an (hx, hy) box. These cover the
% fixed-area ranges of the two existing searches (148-437 nm in hx, 220-505 nm
% in hy) and the test point; holeFabOK prunes the rest. TODO(confirm).
cfg.vars.a  = struct('search', false, 'bounds', [450 700], 'fixed', 618e-9, 'step', nm);
cfg.vars.hx = struct('search', true,  'bounds', [100 500], 'fixed', 270e-9, 'step', nm);
cfg.vars.hy = struct('search', true,  'bounds', [150 600], 'fixed', 270e-9, 'step', nm);
cfg.vars.bw = struct('search', false, 'bounds', [550 850], 'fixed', 600e-9, 'step', nm);
cfg.vars.th = struct('search', false, 'bounds', [350 550], 'fixed', 380e-9, 'step', nm);

% Fixed ellipse area, as the existing optimizer and sweep: hx*hy = A0, so
% hy = A0/hx is DERIVED (not on the nm grid) and hy must not be searched.
% false -> hx and hy are both free (2-D). Free is the default (user decision
% 2026-11-22).
cfg.fixArea     = false;
cfg.A0          = [];         % [] -> hx.fixed*hy.fixed (optimize_hole_unitCell.m:27)
cfg.aspectRange = [];         % [rMin rMax] on hx/hy, or [] for none
                              % (optimize_hole_unitCell.m:28-29 uses [0.3 1.5])

% Process limits: defaults of nanobeam/isFabricable.m:50-52, same convention.
% NOTE nanobeam/CreateNanobeamGeom.m:142 uses 150 nm for the sidewall rule.
cfg.minHoleGap     = 50e-9;   % a - hx, solid bridge between adjacent holes
cfg.minSidewallGap = 148e-9;  % bw - hy, BOTH rails together
cfg.minFeature     = 50e-9;   % smallest hole diameter

cfg.xsect = 'isoFit';         % optimize_hole_unitCell.m:59
cfg.d_in  = 0;                % :69-70, radians
cfg.d_out = 0;
cfg.rxtal = 0;                % user decision 2026-11-22 ("0 deg for now");
                              % optimize_hole_unitCell.m:97 uses 45

% User decision 2026-11-22. (Existing scripts: 5 GHz optimize_hole_unitCell.m:32,
% 7 GHz sweep_hole_unitCell.m:40, 3 GHz test_hole_unitCell.m:57.)
cfg.targetFreq = 8e9;
cfg.sigma      = 0.5e9;       % detuning tolerance (user decision 2026-11-22)
cfg.maxFreq    = 15e9;        % TODO(confirm): no existing source
cfg.solveasym  = 0;           % :79

cfg.kpts     = 10;            % :75-76, :94, :98
cfg.nbands   = 5;
cfg.meshSize = 3;
cfg.max_dof  = 3e6;

cfg.runTag = 'bayesopt_hole_run1';
end

%% ========================================================================
function cfg = holePrepareCfg(cfg)
v = cfg.vars;
if cfg.fixArea
    if v.hy.search
        error('unitCellAdapter_hole:badFixArea', ...
            'cfg.fixArea derives hy = A0/hx, so set cfg.vars.hy.search = false.');
    end
    if isempty(cfg.A0)
        cfg.A0 = v.hx.fixed * v.hy.fixed;
    end
end
if ~isempty(cfg.aspectRange) && ...
        (numel(cfg.aspectRange) ~= 2 || cfg.aspectRange(1) >= cfg.aspectRange(2))
    error('unitCellAdapter_hole:badAspectRange', ...
        'cfg.aspectRange must be [rMin rMax] with rMin < rMax, or [].');
end
end

%% ========================================================================
function D = holeDerive(D, cfg)
if cfg.fixArea
    D.hy = cfg.A0 ./ D.hx;      % area invariant hx*hy = A0
end
end

%% ========================================================================
function ok = holeFabOKD(D, cfg)
ok = holeFabOK(D.a, D.hx, D.hy, D.bw, D.th, cfg);
end

%% ========================================================================
function [ok, margin] = holeFabOK(a, hx, hy, bw, th, cfg)
%HOLEFABOK  The ONE place that decides fabricability for the hole study.
%
%   Metres in, logical out, same size as the broadcast inputs. margin is the
%   worst signed slack in metres.
%     a - hx  >= minHoleGap        bridge between holes   (isFabricable.m:86)
%     bw - hy >= minSidewallGap    both rails together    (isFabricable.m:84)
%     hx, hy  >= minFeature                               (isFabricable.m:87-88)
%     isoFit: th > sqrt(3)*bw/6    the side facets of the pentagon keep a
%                                  positive height (buildHoleUnitCell.m:57)
%     other xsect: bw - 2*th*tan(d_out) > 0               (buildHoleUnitCell.m:59)
%     hy > th*tan(d_in)            hole taper factor > 0  (buildHoleUnitCell.m:37-39)
%   INFERRED from the builder and the nanobeam rules, not from a process spec.
tol  = 1e-15;
sEps = 1e-12;

probe  = a + hx + hy + bw + th;
margin = inf(size(probe));

margin = min(margin, (a - hx)  - cfg.minHoleGap     + tol);
margin = min(margin, (bw - hy) - cfg.minSidewallGap + tol);
margin = min(margin, hx - cfg.minFeature + tol);
margin = min(margin, hy - cfg.minFeature + tol);
if strcmp(cfg.xsect, 'isoFit')
    margin = min(margin, th - sqrt(3)*bw/6 - sEps);
else
    margin = min(margin, bw - 2*th*tan(cfg.d_out) - sEps);
end
margin = min(margin, hy - th*tan(cfg.d_in) - sEps);
if ~isempty(cfg.aspectRange)
    % rMin <= hx/hy <= rMax, written in metres so the margins stay comparable
    margin = min(margin, hx - cfg.aspectRange(1)*hy + tol);
    margin = min(margin, cfg.aspectRange(2)*hy - hx + tol);
end
if cfg.fixArea
    % The derived hy must stay inside its own box.
    bH = cfg.vars.hy.bounds * cfg.vars.hy.step;
    margin = min(margin, min(hy - bH(1), bH(2) - hy) + tol);
end

bad = ~isfinite(probe) | a <= 0 | hx <= 0 | hy <= 0 | bw <= 0 | th <= 0;
margin(bad) = -inf;
ok = margin >= 0;
end

%% ========================================================================
function P = holeBuildP(D, cfg)
%HOLEBUILDP  One design (SI scalars) -> the P struct solveBands expects.
%   Field for field optimize_hole_unitCell.m:59-98.
P.xsect    = cfg.xsect;
P.beamMat  = 'diamond';
P.celltype = 'hole';
P.unitcell = 'rectrangular';        % (sic) as every hole script
P.maxdef   = 0.15;                  % carried for parity; unused by the cell builder
P.oblong   = 0.7265;
P.a  = D.a;
P.hx = D.hx;
P.hy = D.hy;
P.beam_width = D.bw;
P.d_in  = cfg.d_in;
P.d_out = cfg.d_out;
P.th = D.th;
P.nperiod    = 1;
P.holeatedge = 0;
P.kpts   = cfg.kpts;
P.nbands = cfg.nbands;

% isoFit has no z mirror plane, and the builder always cuts the beam in half
% at y = 0 (buildHoleUnitCell.m:130-138). So: one plane, the y plane.
% solveBands sets mbeveny = +1 (solveBands.m:227); ds.sym is the y-even
% sector. With solveasym = 0, ds.full holds the same bands (:276).
P.TwoSymPlanes  = 0;
P.zSymCondition = 0;
P.mbeveny   = 0;
P.mbevenz   = 0;
P.solveasym = cfg.solveasym;

P.completeBandGaps = 1;
P.plotgeom    = 0;
P.savedat     = 1;
P.savebndplot = 0;
P.saveplots   = 0;
P.saveMPH     = 0;
P.bandStruct_2D    = 0;
P.bandStructureDim = 1;
P.airDiskH     = 1000e-9;           % read unconditionally (buildHoleUnitCell.m:16)
P.run_optical  = 0;
P.optical_freq = 200;
P.nbeam        = 2.386;
P.freq     = cfg.targetFreq;
P.meshSize = cfg.meshSize;
P.fixed_bc = 0;
P.anisoMat = 1;
P.rxtal    = cfg.rxtal;
P.max_dof  = cfg.max_dof;
end

%% ========================================================================
function [fCentre, gapFrac] = holeSurrogate(P, cfg)
%HOLESURROGATE  NOT PHYSICS. A Bragg-like 1/a frequency that drifts with the
%   hole size, and a gap bump with an interior optimum near hx/a = 0.45,
%   hy/bw = 0.70, minus a floor so most of the box is gapless. Centred on
%   cfg.targetFreq at the baseline so the dry run stays meaningful for any
%   target.
u = P.hx / P.a;
v = P.hy / P.beam_width;
fCentre = cfg.targetFreq*(618e-9/P.a)*(1 + 0.4*(u - 0.45))*(1 - 0.3*(v - 0.70));
gapFrac = 0.30*exp(-((u - 0.45)/0.15)^2)*exp(-((v - 0.70)/0.15)^2) - 0.03;
end
