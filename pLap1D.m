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

% Create the object function and jacobian functions
obj_func = @(v)objectFunction(v,p,D_in,D_out,points);
jac_func = @(v)jacobianFunction(v,p,D_in,D_out,points);
converge_func = @(v) norm(v(1:end-1)) < threshold;

% Use Laplacian eigenvalue and function as initial guess
[V, D] = eigs(D_out * D_in, startNumber, 'sm');
u_0 = [V(:, startNumber);
       D(startNumber, startNumber)];

% Do Newton's Method
[ u_final, iterations ] = iterativeNewton(u_0, maxIterations, obj_func, jac_func, converge_func);

u = u_final(1:end-1)
lambda = u_final(end)

u_final(end)
plot(points, [0; u_final(1:end-1); 0])

end

