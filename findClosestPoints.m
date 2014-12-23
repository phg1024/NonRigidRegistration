% Find out point to vertex correspondence given a set of points and a mesh
% based on point to surface distance
function correspondence = findClosestPoints(points, S)

nverts = size(S.vertices, 1);
nfaces = size(S.faces, 1);
npoints = size(points, 1);

% compute the normals of each face and each vertex
faceNormal = zeros(nfaces, 3);
vertNormal = zeros(nverts, 3);
for i=1:nfaces
    fvec = S.faces(i,:);
    i1 = fvec(1); i2 = fvec(2); i3 = fvec(3);
    v1 = S.vertices(i1, :);
    v2 = S.vertices(i2, :);
    v3 = S.vertices(i3, :);
    fn = cross(v2-v1, v3-v1);
    w1 = acos(dot(v2-v1, v3-v1)/(norm(v2-v1)*norm(v3-v1)));
    w2 = acos(dot(v3-v2, v1-v2)/(norm(v3-v2)*norm(v1-v2)));
    w3 = acos(dot(v1-v3, v2-v3)/(norm(v1-v3)*norm(v2-v3)));
    faceNormal(i,:) = fn / norm(fn);
    vertNormal(i1, :) = vertNormal(i1, :) + fn * w1;
    vertNormal(i2, :) = vertNormal(i2, :) + fn * w2;
    vertNormal(i3, :) = vertNormal(i3, :) + fn * w3;
end

for i=1:nverts
    vertNormal(i,:) = vertNormal(i,:) / norm(vertNormal(i,:));
end

% find correspondence
dist = zeros(nverts, npoints);
for i=1:npoints
    dist(:, i) = sum((S.vertices - repmat(points(i,:), nverts, 1)) .^2, 2);
end
dist = abs(dist);

% for each column, find the minimum row
[~, correspondence] = min(dist);

end
