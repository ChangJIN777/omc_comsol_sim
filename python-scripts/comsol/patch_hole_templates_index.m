function report = patch_hole_templates_index(files, nDiamond)
% PATCH_HOLE_TEMPLATES_INDEX  Make the diamond refractive index of the hole
% templates a global parameter, in place, via COMSOL LiveLink for MATLAB.
%
%   patch_hole_templates_index()                 % both hole templates, n = 2.386
%   patch_hole_templates_index(files)            % char or cellstr of .mph files
%   patch_hole_templates_index(files, nDiamond)  % another index
%   report = patch_hole_templates_index(...)     % struct array, one per file
%   patch_hole_templates_index('-selftest')      % pure-MATLAB checks, no COMSOL
%
% omc_unitcell.mph (rectangular beam) and omc_unitcell_iso.mph (wedge) were
% made by hand in the GUI (README_template.md) and store the diamond relative
% permittivity as the number 5.76 (n = 2.40). There is no builder to rerun,
% so this function makes the two edits the rib builder
% (build_omc_unitcell_rib.m) already has, and nothing else:
%   1. global parameter   n_diamond = 2.386
%                         'diamond refractive index, telecom (configs/materials.yaml)'
%   2. diamond material   relpermittivity = n_diamond^2   (isotropic)
% The geometry, mesh, physics, studies, solver sequences, datasets and any
% stored solution are not touched and nothing is solved. A stored solution is
% NOT recomputed: it still belongs to the old index until the study is run.
%
% For each file:
%   * skipped, without saving, if it is already patched with this index
%     (idempotent); only the parameter is updated if it is patched with
%     another index;
%   * backed up first to <name>.bak.mph (an existing backup is kept and a
%     time-stamped one is written instead);
%   * the diamond material is found by its label ('Diamond...'), else as the
%     only material whose relative permittivity is not 1;
%   * after saving it is opened again and checked: parameter value, material
%     expression, and that the study / solution / dataset / mesh tags are the
%     ones it had before.
%
% Who needs the patched file: `comsol batch` (scripts/comsol_batch.py can only
% set global parameters), and anyone opening the template in the GUI or with
% mphopen. The Python/MPh drivers do not: they set the index at run time
% (src/comsol_client.py: set_diamond_index), and use the parameter when the
% template has it.
%
% Requires COMSOL with LiveLink for MATLAB connected (mphstart). Holds a
% COMSOL seat while it runs (a few seconds per file).
%
% See also BUILD_OMC_UNITCELL_RIB.

PARAM = 'n_diamond';
DESCR = 'diamond refractive index, telecom (configs/materials.yaml)';

if nargin >= 1 && (ischar(files) || isstring(files)) && strcmp(char(files), '-selftest')
    selfTest();
    report = [];
    return
end
here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(files)
    files = {fullfile(here, 'omc_unitcell.mph'), fullfile(here, 'omc_unitcell_iso.mph')};
end
if nargin < 2 || isempty(nDiamond)
    nDiamond = 2.386;
end
files = cellstr(files);
validateattributes(nDiamond, {'numeric'}, {'scalar', 'real', 'finite', '>=', 1}, ...
    mfilename, 'nDiamond');

import com.comsol.model.*
import com.comsol.model.util.*

report = struct('file', {}, 'action', {}, 'material', {}, 'before', {}, ...
                'backup', {}, 'verified', {});
for i = 1:numel(files)
    file = files{i};
    r = struct('file', file, 'action', 'missing', 'material', '', ...
               'before', '', 'backup', '', 'verified', false);
    if ~exist(file, 'file')
        warning('patch_hole_templates_index:missing', 'Not found, skipped: %s', file);
        report(end + 1) = r; %#ok<AGROW>
        continue
    end

    model = mphopen(file);
    tag = char(model.tag());
    try
        [comp, mat] = findDiamond(model);
        r.material = sprintf('%s/%s (%s)', comp, mat, ...
            char(model.component(comp).material(mat).label()));
        before = readState(model, comp, mat, PARAM);
        r.before = sprintf('%s = %s; relpermittivity = %s', PARAM, ...
            before.paramExpr, strjoin(before.eps, ','));
        r.action = planAction(before.hasParam, before.paramValue, before.eps, ...
                              nDiamond, PARAM);
        if ~strcmp(r.action, 'skip')
            r.backup = backupName(file);
            [ok, msg] = copyfile(file, r.backup);
            if ~ok
                error('patch_hole_templates_index:backup', ...
                    'Could not back up %s: %s', file, msg);
            end
            model.param.set(PARAM, num2str(nDiamond, 16), DESCR);
            if strcmp(r.action, 'patch')
                model.component(comp).material(mat).propertyGroup('def') ...
                    .set('relpermittivity', {[PARAM '^2']});
            end
            mphsave(model, file);
        end
    catch err
        ModelUtil.remove(tag);
        rethrow(err);
    end
    ModelUtil.remove(tag);

    % Verify on the saved file (also for a skipped one).
    model = mphopen(file);
    tag = char(model.tag());
    try
        after = readState(model, comp, mat, PARAM);
    catch err
        ModelUtil.remove(tag);
        rethrow(err);
    end
    ModelUtil.remove(tag);
    problems = verifyState(before, after, nDiamond, PARAM);
    if ~isempty(problems)
        error('patch_hole_templates_index:verify', ...
            '%s is NOT as expected after "%s": %s\nThe original is in %s.', ...
            file, r.action, strjoin(problems, '; '), r.backup);
    end
    r.verified = true;
    report(end + 1) = r; %#ok<AGROW>
    fprintf('%-8s %s\n         material %s\n         was: %s\n', ...
        upper(r.action), file, r.material, r.before);
    if ~isempty(r.backup)
        fprintf('         backup: %s\n', r.backup);
    end
end
end


% ---- COMSOL access -------------------------------------------------------

function [comp, mat] = findDiamond(model)
% Component and material tag of the diamond: by label, else the only
% material whose relative permittivity is not 1.
byLabel = {};
byEps = {};
comps = cell(model.component.tags());
for i = 1:numel(comps)
    mats = cell(model.component(comps{i}).material.tags());
    for j = 1:numel(mats)
        m = model.component(comps{i}).material(mats{j});
        if isDiamondLabel(char(m.label()))
            byLabel(end + 1, :) = {comps{i}, mats{j}}; %#ok<AGROW>
        end
        eps = readEps(m);
        if ~isempty(eps) && ~isUnity(eps)
            byEps(end + 1, :) = {comps{i}, mats{j}}; %#ok<AGROW>
        end
    end
end
if size(byLabel, 1) == 1
    hit = byLabel;
elseif isempty(byLabel) && size(byEps, 1) == 1
    hit = byEps;
else
    error('patch_hole_templates_index:material', ...
        ['Cannot identify the diamond material: %d labelled "Diamond...", ' ...
         '%d with relative permittivity other than 1.'], ...
        size(byLabel, 1), size(byEps, 1));
end
comp = hit{1};
mat = hit{2};
end


function eps = readEps(material)
% Relative permittivity expressions of a material as a cellstr ({} if the
% property is not set). The getter depends on how the value was entered.
eps = {};
pg = material.propertyGroup('def');
getters = {@() pg.getStringArray('relpermittivity'), ...
           @() pg.getStringMatrix('relpermittivity'), ...
           @() pg.getString('relpermittivity')};
for i = 1:numel(getters)
    try
        v = getters{i}();
        v = cellstr(string(v));
        v = v(:)';
        v = v(~cellfun(@isempty, v));
        if ~isempty(v)
            eps = v;
            return
        end
    catch
    end
end
end


function s = readState(model, comp, mat, param)
% Everything that is compared before / after the patch.
names = cell(model.param.varnames());
s.hasParam = any(strcmp(names, param));
s.paramExpr = '(none)';
s.paramValue = NaN;
if s.hasParam
    s.paramExpr = char(model.param.get(param));
    s.paramValue = model.param.evaluate(param);
end
s.eps = readEps(model.component(comp).material(mat));
s.params = sort(setdiff(names, {param}));
s.studies = cellstr(string(cell(model.study.tags())));
s.solutions = cellstr(string(cell(model.sol.tags())));
s.datasets = cellstr(string(cell(model.result.dataset.tags())));
s.meshes = cellstr(string(cell(model.component(comp).mesh.tags())));
geoms = cell(model.component(comp).geom.tags());
s.geomFeatures = {};
if ~isempty(geoms)
    s.geomFeatures = cellstr(string(cell( ...
        model.component(comp).geom(geoms{1}).feature.tags())));
end
end


% ---- pure logic (no COMSOL; covered by '-selftest') ------------------------

function tf = isDiamondLabel(label)
tf = strncmpi(strtrim(label), 'diamond', 7);
end


function tf = isUnity(eps)
% true for 1 or the 3x3 identity given as 9 expressions
v = str2double(eps);
if any(isnan(v))
    tf = false;                              % an expression, not a number
elseif isscalar(v)
    tf = v == 1;
elseif numel(v) == 9
    tf = isequal(v(:)', [1 0 0 0 1 0 0 0 1]);
else
    tf = all(v == 1);
end
end


function tf = usesParam(eps, param)
% Is the permittivity already written in terms of the parameter?
tf = ~isempty(eps) && any(~cellfun(@isempty, regexp(eps, ['\<' param '\>'], 'once')));
end


function action = planAction(hasParam, paramValue, eps, nDiamond, param)
% 'skip'    already patched with this index: nothing is written
% 'update'  already patched, other index: only the parameter changes
% 'patch'   parameter and material expression are written
if hasParam && usesParam(eps, param)
    if abs(paramValue - nDiamond) <= 1e-12 * max(1, abs(nDiamond))
        action = 'skip';
    else
        action = 'update';
    end
else
    action = 'patch';
end
end


function name = backupName(file, existsFcn, stamp)
% <name>.bak.mph next to the file; never overwrites an existing backup.
if nargin < 2
    existsFcn = @(f) exist(f, 'file') ~= 0;
end
if nargin < 3
    stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
end
[folder, base] = fileparts(file);
name = fullfile(folder, [base '.bak.mph']);
if existsFcn(name)
    name = fullfile(folder, [base '.bak_' stamp '.mph']);
end
end


function problems = verifyState(before, after, nDiamond, param)
% What is wrong with the saved file ({} = nothing).
problems = {};
if ~after.hasParam
    problems{end + 1} = ['parameter ' param ' is missing'];
elseif abs(after.paramValue - nDiamond) > 1e-12 * max(1, abs(nDiamond))
    problems{end + 1} = sprintf('%s = %.16g, expected %.16g', param, ...
        after.paramValue, nDiamond);
end
if ~usesParam(after.eps, param)
    problems{end + 1} = ['relpermittivity is "' strjoin(after.eps, ',') ...
        '", not ' param '^2'];
end
kept = {'params', 'studies', 'solutions', 'datasets', 'meshes', 'geomFeatures'};
for i = 1:numel(kept)
    if ~isequal(before.(kept{i}), after.(kept{i}))
        problems{end + 1} = [kept{i} ' changed']; %#ok<AGROW>
    end
end
end


function selfTest()
% Checks of the pure logic above; needs no COMSOL.
p = 'n_diamond';
assert(isDiamondLabel('Diamond (100)') && isDiamondLabel('  diamond'));
assert(~isDiamondLabel('Perfect vacuum') && ~isDiamondLabel('CVD diamond'));
assert(isUnity({'1'}) && isUnity({'1','0','0','0','1','0','0','0','1'}));
assert(~isUnity({'5.76'}) && ~isUnity({'5.76','0','0','0','5.76','0','0','0','5.76'}));
assert(~isUnity({'n_diamond^2'}));
assert(usesParam({'n_diamond^2'}, p) && usesParam({'0', 'n_diamond^2'}, p));
assert(~usesParam({'5.76'}, p) && ~usesParam({}, p) && ~usesParam({'n_diamond2^2'}, p));
% unpatched template
assert(strcmp(planAction(false, NaN, {'5.76'}, 2.386, p), 'patch'));
% parameter added by hand but the material still a number: still to patch
assert(strcmp(planAction(true, 2.386, {'5.76'}, 2.386, p), 'patch'));
% patched: idempotent; other index: parameter only
assert(strcmp(planAction(true, 2.386, {'n_diamond^2'}, 2.386, p), 'skip'));
assert(strcmp(planAction(true, 2.40, {'n_diamond^2'}, 2.386, p), 'update'));
% backup name, never overwriting
f = fullfile('c', 'omc_unitcell_iso.mph');
assert(strcmp(backupName(f, @(x) false, 'S'), fullfile('c', 'omc_unitcell_iso.bak.mph')));
assert(strcmp(backupName(f, @(x) true, 'S'), fullfile('c', 'omc_unitcell_iso.bak_S.mph')));
% verification
b = struct('hasParam', false, 'paramValue', NaN, 'eps', {{'5.76'}}, 'params', {{'a'}}, ...
    'studies', {{'std1'}}, 'solutions', {{'sol1'}}, 'datasets', {{'dset1'}}, ...
    'meshes', {{'mesh1'}}, 'geomFeatures', {{'blk1'}});
a = b;
a.hasParam = true;
a.paramValue = 2.386;
a.eps = {'n_diamond^2'};
assert(isempty(verifyState(b, a, 2.386, p)));
bad = a;
bad.eps = {'5.76'};
assert(isscalar(verifyState(b, bad, 2.386, p)));
bad = a;
bad.datasets = {};
bad.paramValue = 2.4;
assert(numel(verifyState(b, bad, 2.386, p)) == 2);
assert(numel(verifyState(b, b, 2.386, p)) == 2);
fprintf('patch_hole_templates_index: self-test passed\n');
end
