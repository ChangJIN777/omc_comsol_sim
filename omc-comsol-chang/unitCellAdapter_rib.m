function A = unitCellAdapter_rib()
%UNITCELLADAPTER_RIB  'rib' celltype adapter for bayesoptUnitCell.
%
%   A = unitCellAdapter_rib();  cfg = A.defaultCfg();
%   ok = A.fabOKmetres(a, t, s, L, d, th, cfg)    % metres / radians in
%
%   Defaults, bounds and symmetry follow optimize_rib_unitCell.m; the geometry
%   rules mirror validateRibParams(P,'unitcell'). See rib_README.md, section 8.
%
%   See also BAYESOPTUNITCELL, BAYESOPT_RIB, OPTIMIZE_RIB_UNITCELL.
A.name     = 'rib';
A.varNames = {'a', 't', 's', 'L', 'd', 'th'};
A.show = struct( ...
    'a',  struct('scale', 1e9,    'label', 'a_nm',  'unit', 'nm'), ...
    't',  struct('scale', 1e9,    'label', 't_nm',  'unit', 'nm'), ...
    's',  struct('scale', 1e9,    'label', 's_nm',  'unit', 'nm'), ...
    'L',  struct('scale', 1e9,    'label', 'L_nm',  'unit', 'nm'), ...
    'd',  struct('scale', 180/pi, 'label', 'd_deg', 'unit', 'deg'), ...
    'th', struct('scale', 1e9,    'label', 'th_nm', 'unit', 'nm'));
A.defaultCfg  = @ribDefaultCfg;
A.prepareCfg  = @ribPrepareCfg;
A.derive      = @ribDerive;
A.fabOK       = @ribFabOKD;
A.fabOKmetres = @ribFabOK;
A.buildP      = @ribBuildP;
A.surrogate   = @ribSurrogate;
end

%% ========================================================================
function cfg = ribDefaultCfg()
cfg = bayesoptUnitCell('defaults');
nm  = 1e-9;

% Bounds are in units of step. Sources: optimize_rib_unitCell.m:47-60.
%   t >= 90 nm keeps 2*r2 < t with 40 nm fillets; t <= 200 nm keeps
%   r1 < (a - t)/2 at a = 300 nm; s >= 50 nm is the floor of
%   rib_optimize_sweep_diamond.m:84.
cfg.vars.a  = struct('search', true,  'bounds', [240 360], 'fixed', 300e-9, 'step', nm);
cfg.vars.t  = struct('search', true,  'bounds', [ 90 200], 'fixed', 100e-9, 'step', nm);
cfg.vars.s  = struct('search', true,  'bounds', [ 50 200], 'fixed', 100e-9, 'step', nm);
cfg.vars.L  = struct('search', true,  'bounds', [250 500], 'fixed', 370e-9, 'step', nm);
% d: grid of 1 degree, [0, 15] deg (optimize_rib_unitCell.m:59-60). Searching
% d, or fixing it > 0, switches the WHOLE study to r1 = r2 = 0 and the y-even
% sector (see ribPrepareCfg). Not searched by default (user decision,
% 2026-11-22): it is an etch outcome, like th.
cfg.vars.d  = struct('search', false, 'bounds', [0 15], 'fixed', 0, 'step', pi/180);
% th: 110 nm effective / 120 nm nominal (ribDefaultParams.m:77 and header).
% The bounds are NOT from an existing script. TODO(confirm) before searching th.
cfg.vars.th = struct('search', false, 'bounds', [100 130], 'fixed', 110e-9, 'step', nm);

% Hold the footprint A = a*s + t*(2L - s) at cfg.A0 by solving for L
% (optimize_rib_unitCell.m:62-67, :190, :350). L must then not be searched.
cfg.fixFootprint = false;
cfg.A0           = [];        % [] -> from the .fixed values above (94000 nm^2)

cfg.r1 = [];                  % [] -> 40 nm without a sidewall angle, else 0
cfg.r2 = [];                  %       (ribDefaultParams.m:100-106)
cfg.minFeature = 50e-9;       % min(a - t, s, th) >= 50 nm, rib_optimize_sweep_diamond.m:84
cfg.requireStraightFlank = true;   % r1 + r2 < L - s/2; validateRibParams only
                                   % WARNS (rib:filletsMeet, :109-113).
                                   % TODO(confirm) hard rule or not.

cfg.rxtal      = 0;           % beam along <100> (user decision 2026-11-22, "for now")
cfg.targetFreq = 11.2e9;      % optimize_rib_unitCell.m:72
cfg.sigma      = 0.5e9;       % detuning tolerance (user decision 2026-11-22)
cfg.maxFreq    = 25e9;        % TODO(confirm): no existing source
cfg.solveasym  = 0;           % optimize_rib_unitCell.m:148

cfg.kpts     = 10;            % optimize_rib_unitCell.m:135-136, :158, :160
cfg.nbands   = 10;
cfg.meshSize = 3;
cfg.max_dof  = 3e6;

cfg.runTag = 'bayesopt_rib_run1';
end

%% ========================================================================
function cfg = ribPrepareCfg(cfg)
%RIBPREPARECFG  Sidewall mode, fillet defaults, fixed-footprint consistency.
v = cfg.vars;
cfg.usesSidewall = v.d.search || v.d.fixed > 0;
if isempty(cfg.r1), cfg.r1 = 40e-9 * ~cfg.usesSidewall; end
if isempty(cfg.r2), cfg.r2 = 40e-9 * ~cfg.usesSidewall; end
if cfg.usesSidewall && (cfg.r1 > 0 || cfg.r2 > 0)
    error('unitCellAdapter_rib:filletWithSidewall', ...
        ['A sidewall angle is searched or set, but cfg.r1/cfg.r2 are nonzero. ', ...
         'Set them to 0 or [] (validateRibParams: rib:filletWithSidewall).']);
end
if v.d.fixed < 0 || v.d.fixed >= pi/2 || (v.d.search && v.d.bounds(2)*v.d.step >= pi/2)
    error('unitCellAdapter_rib:badSidewallAngle', ...
        'd is in RADIANS in [0, pi/2); cfg.vars.d.bounds are in units of cfg.vars.d.step.');
end
if cfg.fixFootprint
    if v.L.search
        error('unitCellAdapter_rib:badFixFootprint', ...
            'cfg.fixFootprint derives L, so set cfg.vars.L.search = false.');
    end
    if isempty(cfg.A0)
        cfg.A0 = v.a.fixed*v.s.fixed + v.t.fixed*(2*v.L.fixed - v.s.fixed);
    end
end
end

%% ========================================================================
function D = ribDerive(D, cfg)
if cfg.fixFootprint
    % A0 = a*s + t*(2L - s)  ->  L = ((A0 - a*s)/t + s)/2
    D.L = ((cfg.A0 - D.a.*D.s)./D.t + D.s)/2;
end
end

%% ========================================================================
function ok = ribFabOKD(D, cfg)
ok = ribFabOK(D.a, D.t, D.s, D.L, D.d, D.th, cfg);
end

%% ========================================================================
function [ok, margin] = ribFabOK(a, t, s, L, d, th, cfg)
%RIBFABOK  The ONE place that decides fabricability for the rib study.
%
%   Metres (d in radians) in, logical out, same size as the broadcast inputs,
%   so it serves the scan, the N-row XConstraintFcn table and the scalar check
%   in the objective. margin is the worst signed slack in metres.
%
%   Geometry rules = validateRibParams(P,'unitcell'):
%     t < a (:61)   L > s/2 (:68)   d in [0, pi/2) (:21)
%     d = 0: 2*r2 < t (:100), r1 < (a - t)/2 (:104), [r1 + r2 < L - s/2 (:109)]
%     d > 0: r1 = r2 = 0 (:77), t - 2*th*tan(d) > 0 (:83), s - 2*th*tan(d) > 0 (:87)
%   Process floors (cfg.minFeature): a - t, s, th (rib_optimize_sweep_diamond.m:84),
%   and t (inferred, not from a script).
tol  = 1e-15;       % absorbs representation error on >= rules
sEps = 1e-12;       % 1 pm: makes a rule STRICT, as validateRibParams is

probe  = a + t + s + L + d + th;        % errors clearly on incompatible sizes
margin = inf(size(probe));

margin = min(margin, (a - t) - cfg.minFeature + tol);       % gap between ribs
margin = min(margin, t  - cfg.minFeature + tol);
margin = min(margin, s  - cfg.minFeature + tol);
margin = min(margin, th - cfg.minFeature + tol);
margin = min(margin, (L - s/2) - sEps);                     % rib protrudes

isVert = (d == 0) & true(size(probe));

% d = 0: in-plane fillets
mFil = min(t - 2*cfg.r2, (a - t)/2 - cfg.r1) - sEps;
if cfg.requireStraightFlank
    mFil = min(mFil, (L - s/2) - (cfg.r1 + cfg.r2) - sEps);
end
mFil = mFil + zeros(size(probe));
mFil(~isVert) = inf;
margin = min(margin, mFil);

% d > 0: the walls lean in by th*tan(d) at the top face
inset = th .* tan(d);
mSide = min(t - 2*inset, s - 2*inset) - sEps + zeros(size(probe));
if cfg.r1 > 0 || cfg.r2 > 0
    mSide(:) = -inf;                    % fillets cannot be combined with d ~= 0
end
mSide(isVert) = inf;
margin = min(margin, mSide);

if cfg.fixFootprint
    % The derived L must stay inside its own box.
    bL = cfg.vars.L.bounds * cfg.vars.L.step;
    margin = min(margin, min(L - bL(1), bL(2) - L) + tol);
end

bad = ~isfinite(probe) | a <= 0 | t <= 0 | s <= 0 | L <= 0 | th <= 0 ...
    | d < 0 | d >= pi/2;
margin(bad) = -inf;
ok = margin >= 0;
end

%% ========================================================================
function P = ribBuildP(D, cfg)
%RIBBUILDP  One design (SI scalars) -> the P struct solveBands expects.
%   Solver/IO fields follow optimize_rib_unitCell.m:123-160.
%   The design is written BEFORE ribDefaultParams, so airradUC = 5*a follows
%   the current a (ribDefaultParams.m:97-99; the applyDesign rule at
%   optimize_rib_unitCell.m:353-355). L is explicit and P.w is never set, so
%   ribLength sees no L/w conflict (ribLength.m:20-36).
P = struct('a', D.a, 't', D.t, 's', D.s, 'L', D.L, 'd', D.d, 'th', D.th, ...
           'r1', cfg.r1, 'r2', cfg.r2, 'rxtal', cfg.rxtal);
P = ribDefaultParams(P);
P.celltype = 'rib';
P.anisoMat = 1;
P.kpts     = cfg.kpts;
P.nbands   = cfg.nbands;

% One sector into ds.sym (TwoSymPlanes = 0):
%   no sidewall: zSymCondition = 1 -> solveBands sets mbevenz = +1 (:224) and
%                keeps mbeveny = 1: y-even AND z-even, quarter model.
%   sidewall   : zSymCondition = 0 -> mbeveny = +1 (:227), mbevenz stays 0:
%                y-even, both z parities (validateRibParams.m:92 rejects a z
%                plane with d ~= 0). Applied to EVERY point of such a study,
%                d = 0 included, so the objective is one quantity throughout.
P.TwoSymPlanes  = 0;
P.zSymCondition = double(~cfg.usesSidewall);
P.mbeveny   = 1;
P.mbevenz   = 0;
P.solveasym = cfg.solveasym;

P.completeBandGaps = 1;
P.plotgeom    = 0;
P.savedat     = 1;
P.savebndplot = 0;
P.saveplots   = 0;
P.saveMPH     = 0;
P.bandStruct_2D = 0;
P.run_optical = 0;          % runBands refuses the rib air domain (runBands.m:138)
P.freq     = cfg.targetFreq;
P.meshSize = cfg.meshSize;
P.fixed_bc = 0;
P.max_dof  = cfg.max_dof;
end

%% ========================================================================
function [fCentre, gapFrac] = ribSurrogate(P, cfg)
%RIBSURROGATE  NOT PHYSICS. The bump of surrogateBands in
%   optimize_rib_unitCell.m:400-413 (optimum near t = 130, s = 80, L = 400 nm),
%   as a relative gap with a floor so part of the box is gapless, and a mild
%   penalty on the sidewall angle. Centred on cfg.targetFreq at the baseline
%   (a = 300, L = 370, t = 100 nm) so the dry run stays meaningful for any target.
L = ribLength(ribDefaultParams(P), P.a);
fCentre = cfg.targetFreq*(300e-9/P.a)*sqrt(370e-9/L)*(P.t/100e-9)^0.2;
gap = 2e9*exp(-((P.t - 130e-9)/60e-9)^2 - ((P.s - 80e-9)/50e-9)^2 ...
              - ((L - 400e-9)/120e-9)^2) * exp(-(P.d/(20*pi/180))^2);
gapFrac = gap/fCentre - 0.03;
end
