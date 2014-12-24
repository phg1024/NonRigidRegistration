% Find out point to vertex correspondence given a set of points and a mesh
% based on point to surface distance
function [correspondence, weights] = findClosestPoints(points, S)

nverts = size(S.vertices, 1);
nfaces = size(S.faces, 1);
npoints = size(points, 1);

% compute the normals of each face and each vertex
vertNormal = zeros(nverts, 3);
for i=1:nfaces
    fvec = S.faces(i,:);
    i1 = fvec(1); i2 = fvec(2); i3 = fvec(3);
    v1 = S.vertices(i1, :);
    v2 = S.vertices(i2, :);
    v3 = S.vertices(i3, :);
    fn = cross(v2-v1, v3-v1);
    fa = 0.5 * norm(fn);
    v2mv1 = normalize(v2-v1); v3mv1 = normalize(v3-v1);   
    w1 = acos(clamp(dot(v2mv1, v3mv1), -1, 1));
    v3mv2 = normalize(v3-v2); v1mv2 = normalize(v1-v2);
    w2 = acos(clamp(dot(v3mv2, v1mv2), -1, 1));
    v1mv3 = normalize(v1-v3); v2mv3 = normalize(v2-v3);
    w3 = acos(clamp(dot(v1mv3, v2mv3), -1, 1));
    
    vertNormal(i1, :) = vertNormal(i1, :) + fn * fa;
    vertNormal(i2, :) = vertNormal(i2, :) + fn * fa;
    vertNormal(i3, :) = vertNormal(i3, :) + fn * fa;
end

for i=1:nverts
    vertNormal(i,:) = normalize(vertNormal(i,:));
end

% find correspondence
dist = zeros(nverts, npoints);
for i=1:npoints
    dist(:, i) = sum((S.vertices - repmat(points(i,:), nverts, 1)) .^ 2, 2);
end
dist = abs(dist);

% for each column, find the minimum row
[~, correspondence] = min(dist);

% find the weights for each point
weights = abs(sum((vertNormal(correspondence, :) .* normalizeMatrix((S.vertices(correspondence, :) - points))), 2));

end

function v = normalize(v)
nv = norm(v);
if nv > 1e-12
    v = v / nv;
end
end

function m = normalizeMatrix(m)
nrows = size(m, 1);
for i=1:nrows
    m(i,:) = normalize(m(i,:));
end
end

function v = clamp(v, lower, upper)
v = min(max(v, lower), upper);
end
