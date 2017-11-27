function D = outerDifference1D(points)
%outerDifference1D computes the outer difference matrix
%   Computes the outer difference matrix with 0 Dirichlet boundary
%   conditions on the input set of points.
%
% Parameters
%   points: The set of points that the function is approximated at
%
% Outputs:
%   D: The outer difference matrix

% Setup the matrix
numPoints = size(points, 2);
D = zeros(numPoints - 2, numPoints - 1);

% Fill the matrix
for index = 1:numPoints-2
    dist = points(index + 1) - points(index);
    D(index, index) = -1 / dist;
    D(index, index + 1) = 1 / dist;
end

% Make the matrix sparse for efficiency
D = sparse(D);

end

