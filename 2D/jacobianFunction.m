function J = jacobianFunction(u_lambda, p, D_in, D_out, vertWeights)
%jacobianFunction Computes the jacobian of the object function
u = u_lambda(1:end-1);
lambda = u_lambda(end);

numPoints = length(u);
numTri = size(D_in,1)/2;

u_in = D_in * u;
u_norm = vecnorm([u_in(1:numTri) u_in(numTri+1:end)],2,2);
main_diag = [u_norm;u_norm].^(p-2)+(p-2)*[u_norm;u_norm].^(p-4).*u_in;
off_diag = (p-2)*u_norm.^(p-4).*u_in(1:numTri).*u_in(numTri+1:end);
D_norm = spdiags([off_diag off_diag],[-numTri numTri],2*numTri,2*numTri)...
       + spdiags(main_diag,0,2*numTri,2*numTri);
D_lap = D_out * D_norm * D_in;
D_obj = D_lap + lambda * spdiags(abs(u).^(p-2), 0, numPoints, numPoints);

D_con = transpose(abs(u).^(p-2) .* u) .* vertWeights / 3;
J = [D_obj abs(u).^(p-2).*u;
     D_con 0];

end

