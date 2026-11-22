function outTag = addRibCell(wpGeom, sfx, xc, aL, aR, t, L, s, r1, r2, mirrorY)
%ADDRIBCELL  Draw one rib cell (vertical walls, d = 0) in a work-plane geometry.
%
%   outTag = addRibCell(wpGeom, sfx, xc, aL, aR, t, L, s, r1, r2, mirrorY)
%
%   wpGeom   work-plane geometry sequence (e.g. model.geom('g').feature('wp').geom)
%   sfx      tag suffix, unique per cell
%   xc       rib centre in x [m]
%   aL, aR   spine extent left/right of the rib centre [m]. aL = aR = a/2 for
%            the template's 'cellwidth' convention.
%   t, L, s  rib width, rib half-length (axis to tip), spine width [m]
%   r1, r2   concave root / convex tip fillet radii [m]. 0 skips the fillet
%            (COMSOL's Fillet errors on radius 0).
%   mirrorY  true -> also mirror to y < 0 and unite. false -> y >= 0 only
%            (a y-symmetry plane at y = 0).
%   outTag   tag of the object representing the finished cell
%
%   Reproduces the per-cell feature chain of comsol_templates/visible_rib_omc.m
%   (lines 26-64) with the same feature names, box extents and order:
%     r1_  rib [t x L], r3_  half-spine [(aL+aR) x s/2], both drawn at y0 = 4*L
%     union_                     Union, intbnd = false
%     fil2_ (+ fil2_<sfx>_box)   tip fillet r2. Box x in [-(aL+10nm), aR+10nm],
%                                y = y0 + L +-10 nm
%     fil1_ (+ fil1_<sfx>_box)   root fillet r1. Box x = +-(t/2+10nm),
%                                y = y0 + s/2 +-10 nm
%     mov_                       Move by (xc, -y0)
%     mir_, finalunion_          Mirror about y = 0 (keep) + Union, if mirrorY
%   NOTE the template's naming: 'fil2' is the TIP and 'fil1' the ROOT fillet.
%   In P, r2 = tip and r1 = root, matching buildCrossStrip (r1 inner, r2 tips).
%
%   WHY THE y0 = 4*L STAGING. A work-plane BoxSelection sees every object
%   already in the sketch. Drawing each cell far from the beam (and moving it
%   in only after its fillets are done) keeps the fillet boxes from catching
%   corners of cells placed earlier. validateRibParams warns when the
%   staging band would overlap the placed cells.

selDelta = 10e-9;
y0 = 4*L;

rib = wpGeom.feature.create(['r1_',sfx], 'Rectangle');
rib.set('pos', [-t/2 y0]).set('size', [t L]);
spine = wpGeom.feature.create(['r3_',sfx], 'Rectangle');
spine.set('pos', [-aL y0]).set('size', [aL+aR s/2]);

uni = wpGeom.feature.create(['union_',sfx], 'Union');
uni.set('intbnd', false);
uni.selection('input').set({['r1_',sfx] ['r3_',sfx]});
lastTag = ['union_',sfx];

if r2 > 0
    boxTip = wpGeom.feature.create(['fil2_',sfx,'_box'], 'BoxSelection');
    boxTip.set('xmin', -(aL+selDelta)).set('xmax', aR+selDelta);
    boxTip.set('ymin', y0+L-selDelta).set('ymax', y0+L+selDelta);
    boxTip.set('entitydim', 0).set('condition', 'inside');
    filTip = wpGeom.feature.create(['fil2_',sfx], 'Fillet');
    filTip.set('radius', r2);
    % The template export writes selection('pointinsketch'). 'point' is what
    % every other builder in this repo uses. TODO: confirm in mphlaunch that
    % both tip corners come out rounded.
    filTip.selection('point').named(['fil2_',sfx,'_box']);
    lastTag = ['fil2_',sfx];
end

if r1 > 0
    boxRoot = wpGeom.feature.create(['fil1_',sfx,'_box'], 'BoxSelection');
    boxRoot.set('xmin', -(t/2+selDelta)).set('xmax', t/2+selDelta);
    boxRoot.set('ymin', y0+s/2-selDelta).set('ymax', y0+s/2+selDelta);
    boxRoot.set('entitydim', 0).set('condition', 'inside');
    filRoot = wpGeom.feature.create(['fil1_',sfx], 'Fillet');
    filRoot.set('radius', r1);
    filRoot.selection('point').named(['fil1_',sfx,'_box']);
    lastTag = ['fil1_',sfx];
end

mov = wpGeom.feature.create(['mov_',sfx], 'Move');
mov.set('displx', xc).set('disply', -y0);
mov.selection('input').set({lastTag});
outTag = ['mov_',sfx];

if mirrorY
    mir = wpGeom.feature.create(['mir_',sfx], 'Mirror');
    mir.set('keep', true);
    mir.set('axis', [0 1]);      % 2D Mirror: axis = NORMAL of the line, so y -> -y
    mir.selection('input').set({outTag});
    fin = wpGeom.feature.create(['finalunion_',sfx], 'Union');
    fin.set('intbnd', false);
    fin.selection('input').set({['mir_',sfx] outTag});
    outTag = ['finalunion_',sfx];
end
end
