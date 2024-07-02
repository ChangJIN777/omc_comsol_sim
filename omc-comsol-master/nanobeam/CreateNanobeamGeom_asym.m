% function to create asymmetric nanobeam hole geometry from two parameter
% data structures PL, PR for left and right half of beam,
% and outputs the geometry in PR

function PR = CreateNanobeamGeom_asym(PR,PL)

% create geometries
PL = CreateNanobeamGeom(PL);
PR = CreateNanobeamGeom(PR);

% assemble final geometry parameters for asymmetric hole geometry
if PR.holeatctr % do not duplicate center hole from holesR
    hStartIdx = 2;
else
    hStartIdx = 1;
end

% total beam length
beamLen = PL.beamLenHalf + PR.beamLenHalf;

% assemble full geom array = [hx hy xpos ypos]
geomHalfL_xpos = -1*PL.geomHalf(:,3)+PL.beamLenHalf;
PL.geomHalf(:,3) = geomHalfL_xpos;
geomHalfR_xpos = PR.geomHalf(:,3)+PL.beamLenHalf;
PR.geomHalf(:,3) = geomHalfR_xpos;
geom = [flipud(PL.geomHalf(hStartIdx:end,:)); PR.geomHalf];

% save to param data struct
PR.beamLen = beamLen;
PR.geom = geom;
PR.beamLenHalfL = PL.beamLenHalf; % for translating holes in DrawNanobeam

end