%% driver routine for 3D Laplacian deformation
close all; clear all;

% read target mesh and template mesh
template_mesh_path = 'C:\Users\Peihong\Desktop\Data\FaceWarehouse_Data_0\Tester_106\Blendshape\shape_0.obj';


% experiment with a target mesh by sampling point over the mesh
target_mesh_path = 'C:\Users\Peihong\Desktop\Data\FaceWarehouse_Data_0\Tester_106\TrainingPose\';

S = triangulateMesh(loadMesh(template_mesh_path));

landmarks = load('landmarks.mat');
landmarks = landmarks.landmarks + 1;

figure;showMeshWithLandmarks(S, landmarks, 'source');

for i=1:18
tgt_mesh_file = [target_mesh_path, 'pose_', num2str(i), '.obj']
T = triangulateMesh(loadMesh(tgt_mesh_file));

%figure;showMeshWithLandmarks(S, landmarks, 'source');
%figure;showMeshWithLandmarks(T, landmarks, 'target');

% sample points from the target mesh as point cloud
npoints = 3000;
point_cloud = samplePointsFromMesh(T, landmarks, npoints);

figure;showMeshWithPointCloud(T, point_cloud, 'Point Cloud');

lm_points = T.vertices(landmarks,:);

tic;
Td = laplacianDeformation(S, landmarks, lm_points, point_cloud);
toc;

figure;showMeshWithLandmarks(T, landmarks, ['target ', num2str(i)]); print('-dpng','-r300',['target ', num2str(i)]);
figure;showMeshWithLandmarks(Td, landmarks, ['deformed ', num2str(i)]); print('-dpng','-r300',['deformed ', num2str(i)]);
figure;showMeshOverlay(Td, T, ['overlay ', num2str(i)]); print('-dpng','-r300',['overlay ', num2str(i)]);
figure;showMeshError(Td, T, ['error ', num2str(i)]); print('-dpng','-r300',['error ', num2str(i)]);

end
