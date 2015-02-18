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

start_mesh = 1;
nmeshes = 1;
meshes_per_batch = 1;

finished = false;
batchidx = 1;
while ~finished

first_mesh = (batchidx-1)*meshes_per_batch+1;
last_mesh = min(batchidx*meshes_per_batch, nmeshes);

if first_mesh > nmeshes
    break;
end

parfor i=first_mesh:last_mesh
tgt_mesh_file{i} = [target_mesh_path, 'pose_', num2str(i), '.obj']
T{i} = triangulateMesh(loadMesh(tgt_mesh_file{i}));

%figure;showMeshWithLandmarks(S, landmarks, 'source');
%figure;showMeshWithLandmarks(T, landmarks, 'target');

% sample points from the target mesh as point cloud
npoints = 16384;
%point_cloud = samplePointsFromMesh(T, landmarks, npoints);
point_cloud{i} = samplePointsFromMesh2(T{i}, npoints, 1e-2);

figure;showMeshWithPointCloud(T{i}, point_cloud{i}, 'Point Cloud');

lm_points{i} = T{i}.vertices(landmarks,:);

tic;
Td{i} = laplacianDeformation(S, landmarks, lm_points{i}, point_cloud{i});
toc;
end

for i=first_mesh:last_mesh
figure;showMeshWithLandmarks(T{i}, landmarks, ['target ', num2str(i)]); print('-dpng','-r300',['target ', num2str(i)]);
figure;showMeshWithLandmarks(Td{i}, landmarks, ['deformed ', num2str(i)]); print('-dpng','-r300',['deformed ', num2str(i)]);
figure;showMeshOverlay(Td{i}, T{i}, ['overlay ', num2str(i)]); print('-dpng','-r300',['overlay ', num2str(i)]);
figure;showMeshError(Td{i}, T{i}, ['error ', num2str(i)]); print('-dpng','-r300',['error ', num2str(i)]);
end

batchidx = batchidx + 1;
end
