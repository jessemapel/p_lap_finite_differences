function J = jacobianFunction(u, lambda, p, D_in, D_out, points)
%jacobianFunction Computes the jacobian of the object function

numPoints = length(points);

u_in = D_in * u;
u_scaled = spdiag(abs(u_in).^(p-2));
D_lap = D_out * u_scaled;
D_obj = D_lap + lambda * spdiag(abs(u).^(p-2));

dists = zeros(1, numPoints);
for index = 2:numPoints-2
    dists(index) = points(index + 1) - points(index);
dists(1) = dist(2);
dists(end) = dists(end-1);
D_con = transpose(abs(u).^(p-2) .* u) ./ dists;

J = [D_obj abs(u).^(p-2).*u;
     D_con 0];

end

