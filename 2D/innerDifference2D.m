function D = innerDifference2D(vertices, bounds, triangles)
%innerDifference2D computes the inner difference matrix
%   Computes the inner difference matrix with 0 Dirichlet boundary
%   conditions on the input mesh.
%
% Parameters
%   vertices: The N-by-2 array of vertices in the mesh
%   bounds: The N-by-2 array of boundary edges in the mesh
%       where vertices(bounds(j,1),:) and vertices(bounds(j,2),:)
%       are the coordinates of the endpoints of the jth
%       boundary edge
%   triangles: The N-by-3 array of triangles in the mesh where
%       vertices(trinalges(j,1),:), vertices(trinalges(j,2),:),
%       and vertices(trinalges(j,3),:) are the coordinates of
%       the jth triangle.
%
% Outputs:
%   D: The inner difference matrix.

% Setup
% Keep the x and y derivative matrices separate for ease of indexing
numTriangles = size(triangles, 1);
numVertices = size(vertices, 1);
Dx = zeros(numTriangles, numVertices);
Dy = zeros(numTriangles, numVertices);

% Fill the matrix
for index = 1:numTriangles
    verts = vertices(triangles(index,:),:);
    denom = (verts(2,1)-verts(1,1))*(verts(3,2)-verts(1,2))...
           -(verts(2,2)-verts(1,2))*(verts(3,1)-verts(1,1));
    Dx(index,triangles(index,1)) = (verts(3,2)-verts(2,2))/denom;
    Dx(index,triangles(index,2)) = (verts(1,2)-verts(3,2))/denom;
    Dx(index,triangles(index,3)) = (verts(2,2)-verts(1,2))/denom;
    Dy(index,triangles(index,1)) = (verts(2,1)-verts(3,1))/denom;
    Dy(index,triangles(index,2)) = (verts(3,1)-verts(1,1))/denom;
    Dy(index,triangles(index,3)) = (verts(1,1)-verts(2,1))/denom;
end

% Delete the columns for boundary vertices
boundVertices = reshape(bounds,[],1);
Dx(:,boundVertices) = [];
Dy(:,boundVertices) = [];

% Make the matrix sparse for efficiency
D = sparse([Dx; Dy]);

end

