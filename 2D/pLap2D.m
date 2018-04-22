function [ lambda, u, iterations ] = pLap2D( p, maxIterations, ...
                                             vertices, bounds, triangles, ...
                                             threshold, startNumber )
%pLap2D Compute an eigen-value of the 2D p-Laplacian
%   Compute an eigen-value and eigenfunction of the 2D p-Laplacian
%   using iterative Newton's method with finite diferences. First an
%   eigen-function and eigen-value of the Laplacian are computed, then
%   that is used as a starting point for Newton's method. Iteration stops
%   when either the size of the update falls below the threshold, or the
%   maximum number of iterations is reached.
%
% Parameters
%   p: The p value used in the p-Laplacian
%   maxIterations: The maximum number of iterations of Newton's method
%   vertices: The points to approximate the eigenfunction at
%   bounds: The boundary edges as indices into vertices
%   triangles: The mesh triangles as indices into vertices
%   threshold: The convergence threshold
%   startNumber: The number of the eigen-function, eigen-value pair to use
%
% Outputs
%   lambda: The eigen-value
%   u: The eigen-function approximated at points
%   iterations: The number of iterations taken

% Create difference matrices
D_in = innerDifference2D(vertices, bounds, triangles);
D_out = outerDifference2D(vertices, bounds, triangles);

% Compute the vertex weights for numeric integration
numPoints = size(vertices,1);
legOne = vertices(triangles(:,2),:)-vertices(triangles(:,1),:);
legTwo = vertices(triangles(:,3),:)-vertices(triangles(:,1),:);
areas = (legOne(:,1).*legTwo(:,2)-legOne(:,2).*legTwo(:,1))/2;
vertWeights = zeros(numPoints,1);
for index=1:numPoints
    [adjacent,~] = find(triangles==index);
    vertWeights(index) = sum(areas(adjacent));
end
vertWeights(reshape(bounds,[],1)) = [];

% Create the object function and jacobian functions
obj_func = @(v)objectFunction(v,p,D_in,D_out,vertWeights);
jac_func = @(v)jacobianFunction(v,p,D_in,D_out,vertWeights);
converge_func = @(v) norm(v(end)) < threshold;

% Use Laplacian eigenvalue and function as initial guess
[V, D] = eigs(D_out * D_in, startNumber, 'sm');
scale = ( 3*p / sum(vertWeights.*(abs(V(:, startNumber)).^ p)) )^(1/p);
u_0 = [scale * V(:, startNumber);
       D(startNumber, startNumber)];

% Do Newton's Method
[ u_final, iterations ] = iterativeNewton(u_0, maxIterations, obj_func, jac_func, converge_func);

u = u_final(1:end-1);
lambda = u_final(end);

end

