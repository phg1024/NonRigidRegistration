%% driver routine for 3D Laplacian deformation
close all; clear all;

% read target mesh and template mesh
template_mesh_path = 'C:\Users\Peihong\Desktop\Data\FaceWarehouse_Data_0\Tester_1\Blendshape\shape_0.obj';


% experiment with a target mesh by sampling point over the mesh
target_mesh_path = 'C:\Users\Peihong\Desktop\Data\FaceWarehouse_Data_0\Tester_1\TrainingPose\pose_1.obj';


S = triangulateMesh(loadMesh(template_mesh_path));
T = triangulateMesh(loadMesh(target_mesh_path));