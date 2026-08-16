function out = calcFillingFactor(src)
%CALCFILLINGFACTOR Air/dielectric area ratio of a boomerang unit cell.
%
% Computes the in-plane filling factor from a saved *_bds.mat (or from a loaded
% ds/P struct), by reconstructing the same 2D cross-section that
% buildBoomerangUnitCell.m hands to COMSOL and measuring it with polyshape
% boolean operations. Nothing is estimated from a formula - the areas come from
% the actual clipped geometry.
%
%   fillingFactor = area(air hole) / area(dielectric)
%
% NOTE that this is the ratio you asked for, air over REMAINING SOLID, not the
% more commonly quoted air-over-cell. Both are returned; see airFraction. They
% relate as fillingFactor = airFraction/(1-airFraction).
%
% INPUTS
%   src   One of:
%           - path to a *_bds.mat file (char/string)
%           - cell array of such paths, giving a struct array back
%           - a ds struct as saved by solveOpticalBands (uses ds.P)
%           - a bare P struct with fields a, w, r
%
% OUTPUT (struct, or struct array when src is a cell array)
%   fillingFactor    area(air) / area(dielectric)          [dimensionless]
%   airFraction      area(air) / area(cell)                [dimensionless]
%   areaAir          air-hole area clipped to the cell     [m^2]
%   areaAirUnclipped air-hole area before the clip         [m^2]
%   areaDielectric   cell area minus the hole              [m^2]
%   areaCell         rhombic unit-cell area, a^2*sqrt(3)/2 [m^2]
%   armsOverhang     true when an arm reaches past the cell boundary. Under
%                    periodicity that arm merges with the neighbouring cell's,
%                    so the tiled structure is not the isolated boomerang the
%                    parameters describe - worth catching before a long solve.
%   a, w, armLength  geometry echoed back                  [m]
%   filletRadii      [r1 r2] as configured                 [m]
%   filletsIncluded  false - see LIMITATION below
%   source           the file path, or '' when given a struct
%
% GEOMETRY REPRODUCED (buildBoomerangUnitCell.m:45-64)
%   cell   rhombus [0 0; a/2 a*sqrt(3)/2; 3a/2 a*sqrt(3)/2; a 0], side a,
%          area a^2*sqrt(3)/2, centred on (3a/4, a*sqrt(3)/4).
%   hole   union of three w-by-r rectangles, each centred r/2 from the cell
%          centre and rotated 0/120/240 deg about its own centre, so the arms
%          point at 90/210/330 deg and each reaches r from the centre.
%   The builder's Compose formula is 'pol1-r1-r2-r3', i.e. the rectangles are
%   subtracted from the cell, which clips any arm that overhangs the boundary.
%   The intersect() here reproduces that clip, so an over-long arm is not
%   double-counted.
%
% LIMITATION - fillets are NOT included. The builder rounds the hole corners
% with radii P.r1 (inner) and P.r2 (outer) via addFillet; this function measures
% the un-filleted polygon, so areaAir is a slight OVER-estimate. Each rounded
% convex corner removes about (1-pi/4)*rf^2 ~ 0.215*rf^2. At the r1 = r2 = 10 nm
% used in test_Boomerang.m that is ~21 nm^2 per corner against a hole of order
% 1e5 nm^2 - a few parts in ten thousand, far below the mesh-convergence error
% on any band it would be compared against. It becomes worth revisiting if the
% fillet radii are ever raised into the tens of nm.
%
% Requires polyshape (base MATLAB, R2017b or newer). No toolbox needed.
%
% EXAMPLES
%   out = calcFillingFactor('./test/boomerang/08152026/optical_..._bds.mat');
%   fprintf('r = %.4f (air/dielectric), %.1f%% air by area\n', ...
%           out.fillingFactor, 100*out.airFraction);
%
%   % sweep a folder, alongside test_plotOpticalBands
%   f = dir('./test/boomerang/**/*_bds.mat');
%   R = calcFillingFactor(fullfile({f.folder},{f.name}));
%   scatter([R.fillingFactor],[R.airFraction]);

narginchk(1,1);

% Cell array -> recurse, so a whole sweep can be measured in one call.
if iscell(src)
    out = arrayfun(@(k) calcFillingFactor(src{k}),(1:numel(src)).');
    return
end

[P,srcName] = resolveParams(src);

% --- required geometry -------------------------------------------------------
need = {'a','w','r'};
missing = need(~isfield(P,need));
if ~isempty(missing)
    error('calcFillingFactor:missingFields', ...
        ['P is missing %s. This function reconstructs the boomerang cell and ' ...
         'needs the lattice constant a, arm width w and arm length r.'], ...
        strjoin(missing,', '));
end
if isfield(P,'celltype') && ~strcmpi(P.celltype,'boomerang')
    error('calcFillingFactor:wrongCelltype', ...
        ['celltype is ''%s'', but the geometry reconstructed here is the ' ...
         'boomerang cell from buildBoomerangUnitCell.m. Add a branch for ' ...
         '''%s'' before using this number.'],P.celltype,P.celltype);
end

a = P.a;
w = P.w;            % arm WIDTH
rArm = P.r;         % arm LENGTH - note P.r is not the ratio this function returns
if isfield(P,'r1'), rFil1 = P.r1; else, rFil1 = NaN; end
if isfield(P,'r2'), rFil2 = P.r2; else, rFil2 = NaN; end

% --- build the polygons ------------------------------------------------------
cellPgon = polyshape([0 a/2 a*(3/2) a], ...
                     [0 a*sqrt(3)/2 a*sqrt(3)/2 0]);

centre = [a*(3/4), a*sqrt(3)/4];        % cell centre, = builder's hole_pos
% Arm centres sit r/2 from the cell centre along 90/210/330 deg; the literals
% below are the builder's, kept in its form rather than re-derived so the two
% can be diffed by eye.
armCentres = [ centre(1),                      centre(2) + rArm/2 ;
               centre(1) - sqrt(3)*rArm/4,     centre(2) - rArm/4 ;
               centre(1) + sqrt(3)*rArm/4,     centre(2) - rArm/4 ];
armRotDeg  = [0; 120; 240];

holePgon = polyshape();                  % empty
for k = 1:3
    holePgon = union(holePgon, ...
        rectanglePgon(armCentres(k,:),w,rArm,armRotDeg(k)));
end

% Clip to the cell, reproducing the builder's 'pol1-r1-r2-r3' subtraction.
holeInCell = intersect(holePgon,cellPgon);

% --- areas -------------------------------------------------------------------
areaCell         = area(cellPgon);        % == a^2*sqrt(3)/2
areaAir          = area(holeInCell);
areaAirUnclipped = area(holePgon);        % before the clip to the cell
areaDielectric   = areaCell - areaAir;

% If the clip removed anything, at least one arm reaches past the cell boundary.
% That is not merely a bookkeeping detail: the cell is periodic, so an arm
% crossing the edge merges with its neighbour's arm in the tiled structure, and
% the thing being solved is no longer the isolated boomerang the parameters
% describe. Tolerance is relative to the cell so float noise on the boolean
% cannot trip it.
armsOverhang = (areaAirUnclipped - areaAir) > 1e-9*areaCell;

if areaDielectric <= 0
    error('calcFillingFactor:noDielectric', ...
        ['The hole fills the entire cell (air %.3g m^2 >= cell %.3g m^2), so ' ...
         'air/dielectric is undefined. Check a, w and r.'],areaAir,areaCell);
end

out = struct( ...
    'fillingFactor',  areaAir/areaDielectric, ...
    'airFraction',    areaAir/areaCell, ...
    'areaAir',        areaAir, ...
    'areaAirUnclipped',areaAirUnclipped, ...
    'areaDielectric', areaDielectric, ...
    'areaCell',       areaCell, ...
    'armsOverhang',   armsOverhang, ...
    'a',              a, ...
    'w',              w, ...
    'armLength',      rArm, ...
    'filletRadii',    [rFil1 rFil2], ...
    'filletsIncluded',false, ...
    'source',         srcName);
end

%% ======================================================== local functions

function pg = rectanglePgon(centreXY,width,height,rotDeg)
%RECTANGLEPGON COMSOL Rectangle with base='center' and rot, as a polyshape.
%
% COMSOL rotates about the base point, and base is 'center' for all three arms,
% so the rotation is about the rectangle's own centre. size is [width height]
% in the rectangle's local frame before rotation.
c = cosd(rotDeg);
s = sind(rotDeg);
local = [-width/2 -height/2 ;
          width/2 -height/2 ;
          width/2  height/2 ;
         -width/2  height/2 ];
rotated = local*[c s; -s c];             % rows are points, so post-multiply
pg = polyshape(rotated(:,1)+centreXY(1), rotated(:,2)+centreXY(2));
end

function [P,srcName] = resolveParams(src)
%RESOLVEPARAMS Accept a file path, a ds struct, or a bare P struct.
srcName = '';
if ischar(src) || isstring(src)
    srcName = char(src);
    if ~isfile(srcName)
        error('calcFillingFactor:noFile','Not a file: %s',srcName);
    end
    S = load(srcName);
    if isfield(S,'ds')
        P = extractP(S.ds,srcName);
    else
        error('calcFillingFactor:noDs', ...
            '%s contains no variable ds.',srcName);
    end
    return
end
if ~isstruct(src)
    error('calcFillingFactor:badInput', ...
        ['src must be a file path, a cell array of paths, or a struct; got ' ...
         '%s.'],class(src));
end
% A ds struct carries P; a bare P struct carries the geometry directly.
if isfield(src,'P') || isfield(src,'opticalBand')
    P = extractP(src,'');
else
    P = src;
end
end

function P = extractP(ds,srcName)
%EXTRACTP Dig P out of a ds struct saved by solveOpticalBands.
if isfield(ds,'P') && isstruct(ds.P)
    P = ds.P;                                   % wrapper level (current saves)
elseif isfield(ds,'opticalBand') && isstruct(ds.opticalBand) && ...
        isfield(ds.opticalBand,'P') && isstruct(ds.opticalBand.P)
    P = ds.opticalBand.P;                       % runner level (older saves)
else
    if isempty(srcName), srcName = 'the supplied struct'; end
    error('calcFillingFactor:noParams', ...
        ['%s carries no P struct, so the geometry cannot be reconstructed. ' ...
         'Files written before solveOpticalBands started saving ds.P are ' ...
         'affected; pass the P struct directly instead.'],srcName);
end
end
