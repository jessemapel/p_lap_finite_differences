function D = outerDifference2D(vertices, bounds, triangles)
%outerDifference2D computes the outer difference matrix
%   Computes the outer difference matrix with 0 Dirichlet boundary
%   conditions on the input set of points.
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
%   D: The outer difference matrix

% Setup
% Keep the x and y derivative matrices separate for ease of indexing
numTriangles = size(triangles, 1);
numVertices = size(vertices, 1);
Dx = zeros(numVertices, numTriangles);
Dy = zeros(numVertices, numTriangles);

% Fill the matrix
% Fit a quadratic over the midpoints of the adjacent triangles
% then compute the derivative that way.
% If there are not enough adjacent triangles, then fit a plane.
for index = 1:numVertices
    if ismember(index,bounds)
        continue
    end
    [adjacentTriangles,~] = find(triangles==index);
    numAdjacent = size(adjacentTriangles,1);
    vertex = vertices(index,:);
    degree = 2;
    if numAdjacent < 6
        degree = 1;
    end
    
    observations = zeros(numAdjacent,degree*3);
    weightVec = zeros(1,numAdjacent);
    for triIndex = 1:numAdjacent
        triangle = triangles(adjacentTriangles(triIndex),:);
        verts = vertices(triangle,:);
        midPoint = sum(verts)./3;
        observations(triIndex,1) = 1;
        observations(triIndex,[2, 3]) = midPoint;
        if (degree == 2)
            observations(triIndex,[4, 5]) = midPoint.^2;
            observations(triIndex,6) = midPoint(1)*midPoint(2);
        end
        
        adjacentVerts = find(triangle~=index);
        legOne = verts(adjacentVerts(1),:)-vertex;
        legTwo = verts(adjacentVerts(2),:)-vertex;
        weightVec(triIndex) = atan2(norm(legOne(2)*legTwo(1)- legOne(1)*legTwo(2)),dot(legOne,legTwo));
    end
    normalMatrix = observations'*diag(weightVec)*observations;
    % If the normal matrix is close to uninvertible, fit a plane
    if and(rcond(normalMatrix) < 1e-7, degree > 1)
        degree = 1;
        observations(:,4:end) = [];
        normalMatrix = observations'*diag(weightVec)*observations;
    end
    solution = normalMatrix\(observations'*diag(weightVec));
    Dx(index,adjacentTriangles) = solution(2,:);
    Dy(index,adjacentTriangles) = solution(3,:);
    if (degree == 2)
        Dx(index,adjacentTriangles) = Dx(index,adjacentTriangles)...
                                    + 2*vertex(1)*solution(4,:)...
                                    + vertex(2)*solution(6,:);
        Dy(index,adjacentTriangles) = Dy(index,adjacentTriangles)...
                                    + 2*vertex(2)*solution(5,:)...
                                    + vertex(1)*solution(6,:);
    end
end

% Delete the rows for boundary vertices
boundVertices = reshape(bounds,[],1);
Dx(boundVertices,:) = [];
Dy(boundVertices,:) = [];

% Make the matrix sparse for efficiency
D = sparse([Dx Dy]);

end

