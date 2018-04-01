function [ u, iterations ] = iterativeNewton( u_0, maxIterations, ...
                                              objectFunc, jacobianFunc, ...
                                              convergenceTest)
%iterativeNewton Uses Newton's Method to find a zero of a function.
%   Uses iterative Newton's Method to find a zero of an object function.
%
% Parameters
%   u_0: The initial starting vector
%   maxIterations: The maximum number of iterations to take
%   objectFunc: A function that takes a vector and outputs the value of
%               the object function at that vector
%   jacobianFunc: A function that takes a vector and outputs the Jacobian
%                 matrix of the object function at that vector
%   convergenceTest: A function that takes the update step and determines
%                    if iteration has converged
%
% Output
%   u: The final vector that is a zero of the object function
%   iterations: The number of iterations taken to converge

% Loop setup
u = u_0;
iterations = 0;
go = true;

% Iteration loop
while go
    % Evaluate functions
    G = objectFunc(u);
    J = jacobianFunc(u);
    
    % Take step
    chi = J\G;
    u = u - chi;
    
    % Check for convergence
    if convergenceTest(chi)
        go = false;
    end
    
    % Check for max iterations
    iterations = iterations + 1;
    if iterations >= maxIterations
        go = false;
    end
end

end

