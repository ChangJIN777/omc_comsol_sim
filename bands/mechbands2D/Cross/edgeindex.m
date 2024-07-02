% Function to return indices of edges collinear to specified edge.
% Modified from bndindex.
%
% Input arguments: g, p0, varargin
% g: geometry object in model (e.g. model.geom(<geomtag>) or model.geom(<geomtag>).obj(<objtag>) )
% p0: 1x3 array containing coordinates of point in line
% varargin: first element - 1x3 array specifying perpendicular vector to line
%           second element - double specifying maximum distance allowed from 
%               p0 for a point to be considered "in plane". This is to allow 
%               us to find indices of certain boundaries if they happen to be 
%               coplanar with boundaries that we do not want to select; 
%               OR 1x3 array specifying vector - normal is set to be cross 
%               product of varargin{1}-p0 and varargin{2}-p0
%
% Output: bnd
% bnd: boundaries adjacent to plane specified by p0 and varargin
%
% Cleaven Chia, 09/10/2016

function [edg] = edgeindex(g,p0,varargin)

n = varargin{1};    % perpendicular vector to line
dist = Inf;
tol = 1e-10;
if length(varargin) > 1 
    if length(varargin{2}) == 1
        dist = varargin{2}; % if distance is specified
    else
        n = cross(varargin{1}-p0,varargin{2}-p0);   % new normal
    end
end
if(dot(n,n)==0) 
	error('Normal not specified correctly.');
end

% get displacement vectors from p0
vc = g.getVertexCoord();   %returns 3 x NVertices array
vx = vc(1,:) - p0(1);
vy = vc(2,:) - p0(2);
vz = vc(3,:) - p0(3);

% want to find vertices in vc such that (vc-p0) dot n ~ 0
inplane = (abs(vx*n(1)+vy*n(2)+vz*n(3))<max(abs(vx*n(1)+vy*n(2)+vz*n(3)))*tol) ...
          &(sqrt(vx.^2+vy.^2+vz.^2) < dist);
        %returns 1 x NVertices array of 0s and 1s
v_inds = 1:g.getNVertices();
v_inplane_inds = v_inds(inplane);

%get intersection of faces shared by vertices in plane
edg_adjto_vtx = g.getAdj(0,1);  % returns (NVertices+1) x 1 cell; 1st elem is 
                                % null; (n+1)th elem contains (NBndAdj x 1)
                                % array with indices of faces adj to nth
                                % vertex
vtx_adjto_edg = g.getAdj(1,0);  % returns (NBoundaries+1) x 1 cell; 1st elem is 
                                % null; (n+1)th elem contains (NVtxAdj x 1)
                                % array with indices of vertices adj to nth
                                % boundary

% for vertices in plane, get possible boundaries adjacent to these vertices
bnds_to_search = [];
for vi = v_inplane_inds
    bnds_to_search = union(bnds_to_search, edg_adjto_vtx{vi+1});
end

% for each possible boundary, get vertices adjacent to boundary
edg = [];
for bi = transpose(bnds_to_search)
    vtx_bi = vtx_adjto_edg{bi+1};   %gets vertices adj to bnd
    vtx_common = intersect(v_inplane_inds, vtx_bi); %gets vertices common to those in plane
    % if boundary is in plane, all vertices of vtx_bi should be in
    % vtx_common
    % if boundary is not in plane, then vtx_common will be shorter than
    % vtx_bi
    if isequal(sort(vtx_common), vtx_bi)
        edg(end+1) = bi;
    end
end