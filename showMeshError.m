function showMeshError(mesh, ref, name, range, showColorBar)
FV.vertices = mesh.vertices;
FV.faces = mesh.faces;
nverts = size(FV.vertices, 1);
colors = zeros(nverts, 3);
dists = zeros(nverts, 1);

aligned = 1;
if aligned
    center1 = zeros(1, 3);
    center2 = zeros(1, 3);
else
    center1 = mean(mesh.vertices);
    center2 = mean(ref.vertices);
end

for i=1:nverts
    dists(i) = norm((mesh.vertices(i,:)-center1) - (ref.vertices(i,:)-center2));
end
if nargin < 4
    maxDist = max(dists); minDist = 0;
else
    maxDist = range(1); minDist = range(2);
end
distDiff = maxDist - minDist;
for i=1:nverts
    val = min((dists(i)-minDist) / distDiff, 1.0);
    colors(i, :) = interpolate(val);
end

patch(FV, 'FaceVertexCData', colors, 'facecolor', 'interp', 'edgecolor', 'none', ...
    'vertexnormalsmode', 'auto'); 
camlight;
lighting gouraud;
material dull;
axis vis3d;
axis equal;
if nargin < 5
    showColorBar = true;
end
if showColorBar
    colorbar;
    colormap jet(256);
    caxis([minDist maxDist]);
end
title(name);
end

function c = interpolate(val)
nrows = 256;
ridx = max(ceil(val*nrows), 1);
cmap = jet(nrows);
c = cmap(ridx,:);
end