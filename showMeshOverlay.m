function showMeshOverlay(mesh, ref, name)
aligned = 1;

if ~aligned
    center1 = mean(mesh.vertices);
    nverts = size(mesh.vertices, 1);
    FV.vertices = mesh.vertices - repmat(center1, nverts, 1);
else
    FV.vertices = mesh.vertices;
end
FV.faces = mesh.faces;
patch(FV, 'facecolor', [0.5 0.5 0.95], 'edgecolor', 'none', 'vertexnormalsmode', 'auto', 'facealpha', 0.5);

if ~aligned
    center2 = mean(ref.vertices);
    nverts = size(ref.vertices, 1);
    FV.vertices = ref.vertices - repmat(center2, nverts, 1);
else
    FV.vertices = ref.vertices;
end
FV.faces = ref.faces;
patch(FV, 'facecolor', [0.95 0.5 0.5], 'edgecolor', 'none', 'vertexnormalsmode', 'auto', 'facealpha', 0.5);
hold on;
camlight;
lighting gouraud;
material dull;
axis vis3d;
axis equal;
title(name);
end