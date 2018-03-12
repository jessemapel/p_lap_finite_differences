function G = objectFunction(u_lambda, p, D_in, D_out, points)
%objectFunction Computes the object function
u = u_lambda(1:end-1);
lambda = u_lambda(end);

numPoints = length(points);
u_abs = abs(u);

% compute the p-Laplacian
u_in = D_in * u;
u_lap = D_out * (abs(u_in).^(p-2) .* u_in);

% compute the norm
dist = transpose(points(2:end) - points(1:end-1));
trap = spdiags([dist dist], [-1, 0], numPoints - 1, numPoints - 2);
u_norm = sum(trap * (u_abs .^ p)) / (p * 2);

% put together the object function
G = [u_lap + lambda * u_abs.^(p-2) .* u;
     u_norm];
end

