function [h,ff] = plotBoomerangCell(src,varargin)
%PLOTBOOMERANGCELL Outline the dielectric and the air hole of a boomerang cell.
%
% Draws the in-plane cross-section buildBoomerangUnitCell.m builds: the rhombic
% unit cell, the three-armed air hole, and the geometric landmarks that decide
% whether the design is sane - cell centroid, hole centre, the inradius circle
% the arms must stay inside, and the circle the arm corners actually reach.
%
% The geometry comes from calcFillingFactor, so this plot and the filling-factor
% number can never disagree about what is being measured.
%
% INPUTS
%   src   Anything calcFillingFactor accepts: a path to a *_bds.mat, a ds
%         struct, or a bare P struct. Passing P plots a design that has not been
%         simulated yet, which is the point - look before you solve.
%
% NAME-VALUE OPTIONS
%   'Tile'   Draw the six neighbouring cells in outline to show how the pattern
%            repeats. Default false. Useful for seeing whether arms approach
%            each other across the cell boundary.
%   'Axes'   Target axes handle. Default [] = new figure.
%
% OUTPUTS
%   h   struct of graphics handles: cell, hole, centre, inradius, tipRadius.
%   ff  the calcFillingFactor result this plot was drawn from - areas, ratios,
%       armsOverhang, the polyshapes. Returned so a caller that needs both the
%       picture and the numbers pays for the geometry reconstruction ONCE:
%
%           [~,ff] = plotBoomerangCell(P,'Tile',true);
%           if ff.fillingFactor > 1.2, error('too much air'); end
%
%       test_Boomerang.m uses exactly that form for its pre-solve gate.
%
% WHERE THE HOLE SITS - the question this was written to answer:
% buildBoomerangUnitCell.m:67 sets hole_pos = [a*(1/2+1/4), a*sqrt(3)/4], and
% the cell polygon at line 48 has vertices (0,0), (a/2, a*sqrt(3)/2),
% (3a/2, a*sqrt(3)/2), (a, 0). The centroid of those four vertices is
% ((0+a/2+3a/2+a)/4, (0+a*sqrt(3)/2+a*sqrt(3)/2+0)/4) = (3a/4, a*sqrt(3)/4),
% which is hole_pos exactly. The hole is CENTRED in the cell - not at a corner,
% not offset. The plot marks both points on top of each other to make that
% visible rather than merely asserted.
%
% Note P.holeatedge ('1/0 for hole at edge/center of unit cell' in the test
% scripts) is NOT read by buildBoomerangUnitCell - the hole is centred whatever
% it is set to.
%
% Requires polyshape (base MATLAB, R2017b+). No toolbox needed.
%
% EXAMPLE
%   run('test_Boomerang.m');      % or just build P by hand
%   plotBoomerangCell(P,'Tile',true);

narginchk(1,5);

opt = struct('Tile',false,'Axes',[]);
for ii = 1:2:numel(varargin)
    name = validatestring(varargin{ii},fieldnames(opt),mfilename,'option');
    opt.(name) = varargin{ii+1};
end

ff = calcFillingFactor(src);

nm = 1e9;                      % everything is drawn in nm
a  = ff.a;
c  = ff.centre*nm;

if isempty(opt.Axes)
    figure('Name','Boomerang unit cell','NumberTitle','off','Color','w');
    ax = axes;
else
    ax = opt.Axes;
end
hold(ax,'on');

% --- neighbouring cells first, so they sit behind everything -----------------
if opt.Tile
    % Primitive vectors of the hexagonal lattice: a1 = (a,0) is the cell's own
    % edge, a2 = (a/2, a*sqrt(3)/2) the slanted one. Six neighbours share an
    % edge or a vertex with this cell.
    a1 = [a 0];
    a2 = [a/2 a*sqrt(3)/2];
    shifts = [1 0; -1 0; 0 1; 0 -1; 1 -1; -1 1];
    for k = 1:size(shifts,1)
        d = shifts(k,1)*a1 + shifts(k,2)*a2;      % metres, scaled to nm below
        plot(ax,scalePgon(translate(ff.cellPgon,d),nm), ...
            'FaceColor','none','EdgeColor',[0.80 0.80 0.80], ...
            'LineWidth',0.8,'HandleVisibility','off');
        % The neighbour holes are the informative part: they show how close the
        % arms come to merging across the boundary under periodicity.
        plot(ax,scalePgon(translate(ff.holePgon,d),nm), ...
            'FaceColor','none','EdgeColor',[0.93 0.70 0.65], ...
            'LineWidth',1.0,'HandleVisibility','off');
    end
end

% --- dielectric = cell minus hole -------------------------------------------
dielectric = subtract(ff.cellPgon,ff.holePgon);
h.dielectric = plot(ax,scalePgon(dielectric,nm), ...
    'FaceColor',[0.30 0.55 0.80],'FaceAlpha',0.55, ...
    'EdgeColor',[0.10 0.25 0.45],'LineWidth',1.5, ...
    'DisplayName','dielectric (diamond)');

% --- the air hole, outlined on top ------------------------------------------
h.hole = plot(ax,scalePgon(ff.holePgon,nm), ...
    'FaceColor',[1 1 1],'FaceAlpha',0.0, ...
    'EdgeColor',[0.85 0.20 0.10],'LineWidth',2, ...
    'DisplayName','air hole');

% --- the cell outline, so the clip is visible even where the hole touches ----
h.cell = plot(ax,scalePgon(ff.cellPgon,nm), ...
    'FaceColor','none','EdgeColor',[0.25 0.25 0.25], ...
    'LineWidth',1.2,'LineStyle','--','DisplayName','unit cell');

% --- landmarks ---------------------------------------------------------------
% Cell centroid and hole centre are the SAME point; drawn as two markers so a
% future geometry change that separates them shows up immediately instead of
% hiding behind a single dot.
h.centre = plot(ax,c(1),c(2),'k+','MarkerSize',14,'LineWidth',1.5, ...
    'DisplayName','cell centroid');
plot(ax,c(1),c(2),'ko','MarkerSize',9,'LineWidth',1.2, ...
    'DisplayName','hole centre (coincident)');

theta = linspace(0,2*pi,361);
inradius = a*sqrt(3)/4*nm;              % centre to nearest cell edge
tipR     = hypot(ff.w/2,ff.armLength)*nm;   % centre to farthest arm corner

h.inradius = plot(ax,c(1)+inradius*cos(theta),c(2)+inradius*sin(theta), ...
    ':','Color',[0.2 0.5 0.2],'LineWidth',1.3, ...
    'DisplayName',sprintf('cell inradius %.0f nm',inradius));
h.tipRadius = plot(ax,c(1)+tipR*cos(theta),c(2)+tipR*sin(theta), ...
    ':','Color',[0.85 0.45 0.0],'LineWidth',1.3, ...
    'DisplayName',sprintf('arm corner reach %.0f nm',tipR));

% --- annotation --------------------------------------------------------------
axis(ax,'equal');
grid(ax,'on');
box(ax,'on');
xlabel(ax,'x (nm)','FontSize',12);
ylabel(ax,'y (nm)','FontSize',12);
legend(ax,'Location','eastoutside');

margin = inradius - tipR;
title(ax,{ ...
    sprintf('a = %.0f nm, w = %.0f nm, r = %.0f nm', ...
        a*nm,ff.w*nm,ff.armLength*nm), ...
    sprintf(['hole centre = cell centroid = (%.1f, %.1f) nm  |  ' ...
             'fill %.4f  |  clearance %.1f nm'], ...
        c(1),c(2),ff.fillingFactor,margin)}, ...
    'FontSize',11);

hold(ax,'off');
end

%% ======================================================== local functions

function pgOut = scalePgon(pg,s)
%SCALEPGON polyshape with all vertices multiplied by s, for unit conversion.
% scale() exists in newer releases but takes a reference point and is not in
% every version this repo runs on; rebuilding from the vertex list is portable.
% Regions() is used so a multi-region shape survives intact.
rg = regions(pg);
if isempty(rg)
    pgOut = polyshape();
    return
end
pgOut = polyshape();
for k = 1:numel(rg)
    v = rg(k).Vertices;
    pgOut = union(pgOut,polyshape(v(:,1)*s,v(:,2)*s));
end
end
