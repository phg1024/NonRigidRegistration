function showMeshWithPointCloud(mesh, points, name, saveit)
if nargin<4
    saveit = false;
end

FV.vertices = mesh.vertices;
FV.faces = mesh.faces;
patch(FV,'facecolor',[0.5 0.5 0.5], 'edgecolor', 'none', 'vertexnormalsmode', 'auto'); camlight;
lighting gouraud;
material dull;
axis vis3d;
axis equal;
title(name);

hold on;
plot3(points(:,1), points(:,2), points(:,3), 'x');

if saveit
    set(gcf,'units','normalized','outerposition',[0 0 1 1])    
    savefig([name, '.fig']);
    print('-dpng','-r300', name);    
end
end