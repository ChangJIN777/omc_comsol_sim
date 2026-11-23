function fig = plotRibBandTrack(res, ax)
%PLOTRIBBANDTRACK  Summary plot of a trackRibBands run.
%
%   fig = plotRibBandTrack(res)
%   fig = plotRibBandTrack(res, ax)      draw into existing axes
%
%   X-point frequency of the two tracked gap-edge bands against the design
%   index. Bridged and ambiguous designs carry an outline marker (square /
%   triangle) on both bands; lost and failed designs have no frequency and
%   are marked with a cross on the bottom axis and a dotted guide line.

T = res.table;
if nargin < 2 || isempty(ax)
    fig = figure('Color', 'w', 'Name', 'Rib band tracking');
    ax = axes(fig);
else
    fig = ancestor(ax, 'figure');
end
x = T.design;
st = string(T.status);
cLo = [42 120 214]/255;  cUp = [235 104 52]/255;  ink = [0.32 0.32 0.31];

hold(ax, 'on');
plot(ax, x, T.f_lo*1e-9, '-o', 'Color', cLo, 'MarkerFaceColor', cLo, ...
    'MarkerSize', 4, 'LineWidth', 1.5, 'DisplayName', 'band below the gap (at X)');
plot(ax, x, T.f_up*1e-9, '-o', 'Color', cUp, 'MarkerFaceColor', cUp, ...
    'MarkerSize', 4, 'LineWidth', 1.5, 'DisplayName', 'band above the gap (at X)');

marks = {'bridged',   's', 'bridged (midpoint solves)'
         'ambiguous', '^', 'ambiguous (best guess)'};
for m = 1:size(marks, 1)
    k = st == marks{m,1};
    if ~any(k), continue; end
    plot(ax, [x(k); x(k)], [T.f_lo(k); T.f_up(k)]*1e-9, marks{m,2}, ...
        'Color', ink, 'MarkerSize', 11, 'LineWidth', 1.2, 'LineStyle', 'none', ...
        'DisplayName', marks{m,3});
end

if numel(x) <= 25, set(ax, 'XTick', x); end
xlim(ax, [min(x) - 0.5, max(x) + 0.5]);
yl = ylim(ax);
ylim(ax, yl);
k = st == "lost" | st == "failed";
if any(k)
    xk = x(k);
    for m = 1:numel(xk)
        xline(ax, xk(m), ':', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    end
    plot(ax, xk, yl(1)*ones(size(xk)), 'x', 'Color', ink, 'MarkerSize', 10, ...
        'LineWidth', 1.5, 'LineStyle', 'none', 'Clipping', 'off', ...
        'DisplayName', 'lost / failed (no frequency)');
end
hold(ax, 'off');

grid(ax, 'on');
set(ax, 'GridAlpha', 0.12, 'Box', 'off', 'TickDir', 'out');
xlabel(ax, 'design index');
ylabel(ax, 'frequency at X (GHz)');
title(ax, sprintf('Tracked gap-edge bands: %d designs, %d bridge solves', ...
    height(T), res.nBridgeSolves));
legend(ax, 'Location', 'best', 'Box', 'off');
end
