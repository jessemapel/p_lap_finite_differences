function [ lambda, u, iterations ] = pLap1D( p, maxIterations, points, ...
                                             threshold, startNumber )
%pLap1D Compute an eigen-value of the 1D p-Laplacian
%   Compute an eigen-value and eigenfunction of the 1D p-Laplacian
%   using iterative Newton's method with finite diferences. First an
%   eigen-function and eigen-value of the Laplacian are computed, then
%   that is used as a starting point for Newton's method. Iteration stops
%   when either the size of the update falls below the threshold, or the
%   maximum number of iterations is reached.
%
% Parameters
%   p: The p value used in the p-Laplacian
%   maxIterations: The maximum number of iterations of Newton's method
%   points: The points to approximate the eigenfunction at
%   threshold: The convergence threshold
%   startNumber: The number of the eigen-function, eigen-value pair to use
%
% Outputs
%   lambda: The eigen-value
%   u: The eigen-function approximated at points
%   iterations: The number of iterations taken

% Create difference matrices
D_in = innerDifference1D(points);
D_out = outerDifference1D(points);


end

function [G] = objectFunc(u, lambda, p, D_in, D_out, points)
%objectFunc Computes the object function

% compute several things for later
numPoints = length(points);
u_abs = abs(u);

% compute the p-Laplacian
u_in = D_in * u;
u_lap = D_out * (abs(u_in).^(p-2) .* u_in);

% compute the norm
dist = (points(2:end) - points(1:end-1));
trap = spdiags([dist dist], [-1, 0], numPoints - 1, numPoints - 2);
u_norm = sum(trap * (u_abs .^ p)) / (p * 2);

% put together the object function
G = [u_lap + lambda * u_abs.^(p-2) .* u;
     u_norm];

end

function [J] = jacobianFunc(u, lambda, p, D_in, D_out, points)
%jacobianFunc Compute the Jacobian of the object function


end

