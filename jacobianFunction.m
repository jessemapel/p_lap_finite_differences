function J = jacobianFunction(u_lambda, p, D_in, D_out, points)
%jacobianFunction Computes the jacobian of the object function
u = u_lambda(1:end-1);
lambda = u_lambda(end);

numPoints = length(points);

u_in = D_in * u;
u_scaled = spdiags(abs(u_in).^(p-2), 0, numPoints-1, numPoints-1);
D_lap = D_out * u_scaled * D_in;
D_obj = D_lap + lambda * spdiags(abs(u).^(p-2), 0, numPoints-2, numPoints-2);

dists = transpose(points(2:end) - points(1:end-1));
trap_dist = 0.5 * transpose(dists(2:end) + dists(1:end-1));
D_con = transpose(abs(u).^(p-2) .* u) ./ trap_dist;
J = [D_obj abs(u).^(p-2).*u;
     D_con 0];

end

