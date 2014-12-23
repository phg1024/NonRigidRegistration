%% driver routine for 3D Laplacian deformation
close all; clear all;

% read target mesh and template mesh
template_mesh_path = 'C:\Users\Peihong\Desktop\Data\FaceWarehouse_Data_0\Tester_1\Blendshape\shape_0.obj';


% experiment with a target mesh by sampling point over the mesh
target_mesh_path = 'C:\Users\Peihong\Desktop\Data\FaceWarehouse_Data_0\Tester_1\TrainingPose\pose_1.obj';


S = triangulateMesh(loadMesh(template_mesh_path));
T = triangulateMesh(loadMesh(target_mesh_path));

landmarks = load('landmarks.mat');
landmarks = landmarks.landmarks + 1;

%figure;showMeshWithLandmarks(S, landmarks, 'source');
%figure;showMeshWithLandmarks(T, landmarks, 'target');

% sample points from the target mesh as point cloud
npoints = 1024;
point_cloud = samplePointsFromMesh(T, landmarks, npoints);

figure;showMeshWithPointCloud(T, point_cloud, 'Point Cloud');

lm_points = T.vertices(landmarks,:);

tic;
Td = laplacianDeformation(S, T, landmarks, lm_points, point_cloud);
toc;

figure;showMeshWithLandmarks(S, landmarks, 'source');
figure;showMeshWithLandmarks(T, landmarks, 'target');
figure;showMeshWithLandmarks(Td, landmarks, 'deformed');
figure;showMeshOverlay(Td, T, 'overlay');
figure;showMeshError(Td, T, 'error');
