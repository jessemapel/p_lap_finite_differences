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
u_norm = trapz(points(2:end-1), (u_abs .^ p)) / p;

% put together the object function
G = [u_lap + lambda * u_abs.^(p-2) .* u;
     u_norm - 1];
end

