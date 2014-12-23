function points = samplePointsFromMesh(mesh, landmarks, npoints)

% filter out the front part
valid_indices = find(mesh.vertices(:,3)>-0.25);

% randomly choose some points
pindices = valid_indices(randperm(length(valid_indices), npoints));
pindices = setdiff(pindices, landmarks);

points = mesh.vertices(pindices, :);

end