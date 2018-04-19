function G = objectFunction( u_lambda, p, D_in, D_out,...
                             vertices, bounds, triangles)
%objectFunction Computes the object function
u = u_lambda(1:end-1);
lambda = u_lambda(end);

numIntPoints = length(u);
numPoints = size(vertices,1);
u_abs = abs(u);

% compute the p-Laplacian
u_in = D_in * u;
u_scale = transpose(vecnorm([u_in(1:numIntPoints)';u_in(numIntPoints+1:end)']));
u_lap = D_out * ([u_scale;u_scale].^(p-2) .* u_in);

% compute the norm
legOne = vertices(triangles(:,2),:)-vertices(triangles(:,1),:);
legTwo = vertices(triangles(:,3),:)-vertices(triangles(:,1),:);
areas = (legOne(:,1).*legTwo(:,2)-legOne(:,2).*legTwo(:,1))./2;
vertWeights = zeros(numPoints,1);
for index=1:numPoints
    [adjacent,~] = find(triangles==index);
    vertWeights(index) = sum(areas(adjacent));
end
vertWeights(reshape(bounds,[],1)) = [];
u_norm = sum(vertWeights.*(u_abs.^ p))/(3*p);

% put together the object function
G = [u_lap + lambda * u_abs.^(p-2) .* u;
     u_norm - 1];
end

