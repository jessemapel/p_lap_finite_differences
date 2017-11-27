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


end

