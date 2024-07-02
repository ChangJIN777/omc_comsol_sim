function [model,P] = drawCross(model,P)
%DRAWCROSS Summary of this function goes here
%   Detailed explanation goes here
a = P.a;
hc = P.hc;
wc = P.wc;
th = P.th;
r1 = P.r1;
r2 = P.r2;
abssym = abs(P.mbevenz);

geom = model.component('comp1').geom.create('geom1', 3);
wp = geom.create('wp1', 'WorkPlane');
wp.set('unite', true);
if ~abssym
    wp.set('quickz',-th/2);
end

wp.geom.create('sq1', 'Square');
wp.geom.feature('sq1').set('base', 'center');
wp.geom.feature('sq1').set('size', a);
wp.geom.create('r1', 'Rectangle');
wp.geom.feature('r1').set('base', 'center');
wp.geom.feature('r1').set('size', [hc wc]);
wp.geom.create('r2', 'Rectangle');
wp.geom.feature('r2').set('base', 'center');
wp.geom.feature('r2').set('size', [wc hc]);
wp.geom.create('dif1', 'Difference');
wp.geom.feature('dif1').selection('input').set({'sq1'});
wp.geom.feature('dif1').selection('input2').set({'r1' 'r2'});
wp.geom.create('fil1', 'Fillet');
wp.geom.feature('fil1').set('radius', r1);
wp.geom.feature('fil1').selection('point').set('dif1(1)', [3 4 5 8 9 12 13 14]);
wp.geom.create('fil2', 'Fillet');
wp.geom.feature('fil2').set('radius', r2);
wp.geom.feature('fil2').selection('point').set('fil1(1)', [8 9 16 17]);

geom.create('ext1', 'Extrude');
geom.feature('ext1').selection('input').set({'wp1'});
if abssym
geom.feature('ext1').setIndex('distance', th/2, 0);
else
geom.feature('ext1').setIndex('distance', th, 0);
end
geom.run;
geom.run('fin');


%% Making selections

P.xEnd1 =  bndindex(geom, [-a/2 0 0], [1 0 0]);
P.xEnd2 = bndindex(geom, [ a/2 0 0], [1 0 0]);

P.yEnd1 =  bndindex(geom, [0 -a/2 0], [0 1 0]);
P.yEnd2 = bndindex(geom, [0 a/2 0], [0 1 0]);
% Note that this will return no indices if there is no boundary at z=0
P.zEnd = bndindex(geom, [0 0 0], [0 0 1]);

end

