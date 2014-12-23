function showMesh(mesh, name, saveit)
if nargin<3
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

if saveit
    set(gcf,'units','normalized','outerposition',[0 0 1 1])    
    savefig([name, '.fig']);
    print('-dpng','-r300', name);    
end
end