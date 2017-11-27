function D = innerDifference1D(points)
%innerDifference1D computes the inner difference matrix
%   Computes the inner difference matrix with 0 Dirichlet boundary
%   conditions on the input set of points.
%
% Parameters
%   points: The set of points that the function is approximated at
%
% Outputs:
%   D: The inner difference matrix

% Setup the matrix
numPoints = size(points, 2);
D = zeros(numPoints - 1, numPoints - 2);

% Boundary conditions
D(1, 1) = 1 / (points(2) - points(1));
D(end, end) = -1 / (points(end) - points(end - 1));

% Fill the matrix
for index = 2:numPoints-2
    dist = points(index + 1) - points(index);
    D(index, index - 1) = -1 / dist;
    D(index, index) = 1 / dist;
end

% Make the matrix sparse for efficiency
D = sparse(D);

end

