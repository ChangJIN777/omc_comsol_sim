function outTag = addRibPrisms(geom, sfx, xc, xL, xR, t, L, s, th, d, halfY)
%ADDRIBPRISMS  One rib cell with sidewall angle d, as 3D tapered prisms.
%
%   outTag = addRibPrisms(geom, sfx, xc, xL, xR, t, L, s, th, d, halfY)
%
%   geom     3D geometry sequence (model.geom('...'))
%   sfx      tag suffix, unique per cell
%   xc       rib centre in x [m]; xL, xR spine segment ends in x [m]
%   t, L, s  rib width (x), rib half-length (axis to tip), spine width (y),
%            all measured at the BOTTOM face z = -th/2 [m]
%   th       slab thickness [m]; the slab spans z in [-th/2, th/2]
%   d        sidewall angle in radians (> 0). Every wall leans inward by
%            th*tan(d) between the bottom and the top face.
%   halfY    true -> build y >= 0 only (y = 0 symmetry plane)
%   outTag   tag of the united cell solid
%
%   This is the sidewall construction of the former buildRibUnitCell_LN
%   (trapezoid cross-sections intersected with the footprint), applied to the
%   centred-rib cell:
%     spine  = yz trapezoid  [s at the bottom -> s - 2*th*tan(d) at the top]
%              extruded along x over [xL, xR]                (LN: wp2/ext2)
%     rib    = xz trapezoid  [t -> t - 2*th*tan(d)] extruded along y   (LN: wp4/ext4)
%              INTERSECTED with
%              yz trapezoid  [tip at L -> L - th*tan(d)] extruded along x (LN: wp3/ext3)
%     cell   = Union(spine, rib), intbnd off
%   The union gives tapered spine walls between the ribs and tapered rib
%   flanks above the spine. That is the LN result, so the cross-section area
%   is a*s(u) + t(u)*(2L - s) at height u, with s(u) = s - 2u*tan(d) and
%   t(u) = t - 2u*tan(d).
%
%   No in-plane fillets: validateRibParams rejects r1/r2 > 0 with d > 0.
%   Work-plane conventions as in the LN builder: 'yz' local (y, z), extrudes
%   toward +x; 'xz' local (x, z), extrudes toward -y unless reversed.

k  = tan(d);
z0 = -th/2;
z1 =  th/2;
sT = s/2 - k*th;          % spine half-width at the top
tT = t/2 - k*th;          % rib half-width at the top
LT = L - k*th;            % rib tip at the top

%% spine: yz trapezoid extruded along +x from xL to xR
wpS = geom.feature.create(['wpS_',sfx], 'WorkPlane');
wpS.set('planetype', 'quick').set('quickplane', 'yz').set('quickx', xL);
polS = wpS.geom.feature.create('pol', 'Polygon');
polS.set('source', 'table');
if halfY
    polS.set('table', [0 z0; s/2 z0; sT z1; 0 z1]);
else
    polS.set('table', [-s/2 z0; s/2 z0; sT z1; -sT z1]);
end
extS = geom.feature.create(['extS_',sfx], 'Extrude');
extS.selection('input').set({['wpS_',sfx]});
extS.setIndex('distance', xR - xL, 0);

%% rib, x taper: xz trapezoid extruded along +y (reverse of the xz normal)
if halfY
    y0 = 0;   yLen = L;
else
    y0 = -L;  yLen = 2*L;
end
wpX = geom.feature.create(['wpX_',sfx], 'WorkPlane');
wpX.set('planetype', 'quick').set('quickplane', 'xz').set('quicky', y0);
polX = wpX.geom.feature.create('pol', 'Polygon');
polX.set('source', 'table');
polX.set('table', [xc-t/2 z0; xc+t/2 z0; xc+tT z1; xc-tT z1]);
extX = geom.feature.create(['extX_',sfx], 'Extrude');
extX.selection('input').set({['wpX_',sfx]});
extX.setIndex('distance', yLen, 0);
extX.set('reverse', true);

%% rib, tip taper: yz trapezoid extruded along +x across the rib
margin = 10e-9;
wpT = geom.feature.create(['wpT_',sfx], 'WorkPlane');
wpT.set('planetype', 'quick').set('quickplane', 'yz').set('quickx', xc - t/2 - margin);
polT = wpT.geom.feature.create('pol', 'Polygon');
polT.set('source', 'table');
if halfY
    polT.set('table', [0 z0; L z0; LT z1; 0 z1]);
else
    polT.set('table', [-L z0; L z0; LT z1; -LT z1]);
end
extT = geom.feature.create(['extT_',sfx], 'Extrude');
extT.selection('input').set({['wpT_',sfx]});
extT.setIndex('distance', t + 2*margin, 0);

ribInt = geom.feature.create(['ribInt_',sfx], 'Intersection');
ribInt.selection('input').set({['extX_',sfx] ['extT_',sfx]});

%% cell = spine U rib
outTag = ['ribCell_',sfx];
uni = geom.feature.create(outTag, 'Union');
uni.selection('input').set({['extS_',sfx] ['ribInt_',sfx]});
uni.set('intbnd', false);
end
