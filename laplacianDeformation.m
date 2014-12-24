function Td = laplacianDeformation(S, landmarks, lm_points, point_cloud)

Td = S;

[nverts, ~] = size(S.vertices);
[nfaces, ~] = size(S.faces);

% find the neighbor information for every vertex on the source mesh
N = cell(nverts, 1);
for i=1:nfaces
    f = S.faces(i,:);
    
    N{f(1)} = union(N{f(1)}, [f(2) f(3)]);
    N{f(2)} = union(N{f(2)}, [f(1) f(3)]);
    N{f(3)} = union(N{f(3)}, [f(1) f(2)]);
end

v = S.vertices;

% compute delta_i
delta = zeros(nverts, 3);
for i=1:nverts
    Ni = N{i};
    nsum = 0;
    for j=1:length(Ni)
        nsum = nsum + v(Ni(j), :);
    end
    delta(i,:) = v(i,:) - nsum / length(Ni);
end

% compute the A matrices
A = cell(nverts, 1);
for i=1:nverts
    Vi = [v(i,:)', skewMatrix(v), eye(3)];
    Ai= zeros(3*length(Ni), 7);
    Ni = N{i};
    for j=1:length(Ni)
        jstart = (j-1)*3+1; jend=j*3;
        vij = v(Ni(j),:)';
        Ai(jstart:jend, :) = [vij, skewMatrix(vij), eye(3)];
    end
    A{i} = [Vi;Ai];
end

w_icp_step = size(lm_points, 1) / size(point_cloud, 1); E_total = realmax;
%w_dist_step = 100.0 * size(lm_points, 1) / size(Td.vertices, 1);
w_icp = 0.0;
w_data = 1.0;
%w_dist = 10.0*w_dist_step;
w_dist = 1.0;
iters = 0; eps1 = 1e-5; eps2 = 1e-8;
itmax = 10;
while iters < itmax    
    iters = iters + 1;
    fprintf('iteration %d\n', iters);
    
    % find icp correspondence
    [icp_correspondence, icp_weights] = findClosestPoints(point_cloud, Td);
    
%     figure; plot(icp_weights);
%     
%     figure; hold on;
%     showMeshWithPointCloud(Td, point_cloud, 'Point Cloud');
%     plot3(Td.vertices(icp_correspondence,1), Td.vertices(icp_correspondence,2), Td.vertices(icp_correspondence,3), 'ro');
%     hold off;
    
    % compute error
    % data term
    E_data = mean(sum((Td.vertices(landmarks,:) - lm_points).^2, 2));
    % icp term
    E_icp = mean(sum((Td.vertices(icp_correspondence,:) - point_cloud).^2, 2));
    % total error
    E0 = E_total;
    E_total = E_data + E_icp;
    E(iters) = E_total;
    E_diff = abs(E_total - E0);
    fprintf('E = %.6f\n', E_total);
    if E_total < eps1 || E_diff < eps2
        fprintf('total error smaller than threshold. stop.\n');
        break;
    end
    
    % icp term
    disp('assembling icp term');
    nicp = size(point_cloud, 1); 
    M_icp_i = zeros(1, nicp*3); M_icp_j = zeros(1, nicp*3); M_icp_v = zeros(1, nicp*3);
    b_icp = reshape((point_cloud .* repmat(icp_weights, 1, 3))', 3*nicp, 1);
    idx = 1;
    for i=1:nicp
        dstart = (icp_correspondence(i)-1)*3;
        wi = icp_weights(i);
        M_icp_i(idx) = idx; M_icp_j(idx) = dstart+1; M_icp_v(idx) = wi; idx = idx + 1;
        M_icp_i(idx) = idx; M_icp_j(idx) = dstart+2; M_icp_v(idx) = wi; idx = idx + 1;
        M_icp_i(idx) = idx; M_icp_j(idx) = dstart+3; M_icp_v(idx) = wi; idx = idx + 1;
    end
    M_icp = sparse(M_icp_i, M_icp_j, M_icp_v, 3*nicp, nverts*3);
    
    % data term
    disp('assembling data term');
    ndata = size(lm_points, 1);
    M_data_i = zeros(1, ndata*3); M_data_j = zeros(1, ndata*3); M_data_v = zeros(1, ndata*3);
    b_data = reshape(lm_points', 3*ndata, 1);
    idx = 1;
    for i=1:ndata
        dstart = (landmarks(i)-1)*3;
        M_data_i(idx) = idx; M_data_j(idx) = dstart+1; M_data_v(idx) = 1; idx = idx + 1;
        M_data_i(idx) = idx; M_data_j(idx) = dstart+2; M_data_v(idx) = 1; idx = idx + 1;
        M_data_i(idx) = idx; M_data_j(idx) = dstart+3; M_data_v(idx) = 1; idx = idx + 1;
    end
    M_data = sparse(M_data_i, M_data_j, M_data_v, 3*ndata, nverts*3);
    
    disp('assembling distortion term');
    M_dist_i = zeros(1, nverts*27); M_dist_j = zeros(1, nverts*27); M_dist_v = zeros(1, nverts*27);
    b_dist = zeros(nverts*3, 1);
    idx = 1;
    for i=1:nverts
        Ti = makeDmatrix(delta(i,:)) * ((A{i}'*A{i}) \ A{i}');
        Ni = N{i};
        
        istart = (i-1)*3;
        % vi
        % laplacian part
        M_dist_i(idx) = istart+1; M_dist_j(idx) = istart+1; M_dist_v(idx) = 1; idx = idx + 1;
        M_dist_i(idx) = istart+2; M_dist_j(idx) = istart+2; M_dist_v(idx) = 1; idx = idx + 1;
        M_dist_i(idx) = istart+3; M_dist_j(idx) = istart+3; M_dist_v(idx) = 1; idx = idx + 1;
        
        % deformation part
        M_dist_i(idx) = istart+1; M_dist_j(idx) = istart+1; M_dist_v(idx) = -Ti(1, 1); idx = idx + 1;
        M_dist_i(idx) = istart+1; M_dist_j(idx) = istart+2; M_dist_v(idx) = -Ti(1, 2); idx = idx + 1;
        M_dist_i(idx) = istart+1; M_dist_j(idx) = istart+3; M_dist_v(idx) = -Ti(1, 3); idx = idx + 1;
        
        M_dist_i(idx) = istart+2; M_dist_j(idx) = istart+1; M_dist_v(idx) = -Ti(2, 1); idx = idx + 1;
        M_dist_i(idx) = istart+2; M_dist_j(idx) = istart+2; M_dist_v(idx) = -Ti(2, 2); idx = idx + 1;
        M_dist_i(idx) = istart+2; M_dist_j(idx) = istart+3; M_dist_v(idx) = -Ti(2, 3); idx = idx + 1;
        
        M_dist_i(idx) = istart+3; M_dist_j(idx) = istart+1; M_dist_v(idx) = -Ti(3, 1); idx = idx + 1;
        M_dist_i(idx) = istart+3; M_dist_j(idx) = istart+2; M_dist_v(idx) = -Ti(3, 2); idx = idx + 1;
        M_dist_i(idx) = istart+3; M_dist_j(idx) = istart+3; M_dist_v(idx) = -Ti(3, 3); idx = idx + 1;
        
        % N(vi)
        w = -1.0 / length(Ni);
        for j=1:length(Ni)
            jstart = (Ni(j)-1)*3;
            
            % laplacian part
            M_dist_i(idx) = istart+1; M_dist_j(idx) = jstart+1; M_dist_v(idx) = w; idx = idx + 1;
            M_dist_i(idx) = istart+2; M_dist_j(idx) = jstart+2; M_dist_v(idx) = w; idx = idx + 1;
            M_dist_i(idx) = istart+3; M_dist_j(idx) = jstart+3; M_dist_v(idx) = w; idx = idx + 1;
            
            % deformation part
            M_dist_i(idx) = istart+1; M_dist_j(idx) = jstart+1; M_dist_v(idx) = -Ti(1, j*3+1); idx = idx + 1;
            M_dist_i(idx) = istart+1; M_dist_j(idx) = jstart+2; M_dist_v(idx) = -Ti(1, j*3+2); idx = idx + 1;
            M_dist_i(idx) = istart+1; M_dist_j(idx) = jstart+3; M_dist_v(idx) = -Ti(1, j*3+3); idx = idx + 1;
            
            M_dist_i(idx) = istart+2; M_dist_j(idx) = jstart+1; M_dist_v(idx) = -Ti(2, j*3+1); idx = idx + 1;
            M_dist_i(idx) = istart+2; M_dist_j(idx) = jstart+2; M_dist_v(idx) = -Ti(2, j*3+2); idx = idx + 1;
            M_dist_i(idx) = istart+2; M_dist_j(idx) = jstart+3; M_dist_v(idx) = -Ti(2, j*3+3); idx = idx + 1;
            
            M_dist_i(idx) = istart+3; M_dist_j(idx) = jstart+1; M_dist_v(idx) = -Ti(3, j*3+1); idx = idx + 1;
            M_dist_i(idx) = istart+3; M_dist_j(idx) = jstart+2; M_dist_v(idx) = -Ti(3, j*3+2); idx = idx + 1;
            M_dist_i(idx) = istart+3; M_dist_j(idx) = jstart+3; M_dist_v(idx) = -Ti(3, j*3+3); idx = idx + 1;
        end
    end
    M_dist = sparse(M_dist_i, M_dist_j, M_dist_v, nverts*3, nverts*3);
    
    disp('sovling least square problem');
    M = [M_data * w_data; M_icp * w_icp; M_dist * w_dist];
    b = [b_data * w_data; b_icp * w_icp; b_dist * w_dist];
    
    v_new = (M'*M) \ (M'*b);
    v_new = reshape(v_new, 3, nverts)';
    disp('done');
    
    Td.vertices = v_new;
    
    % increase w_icp
    w_icp = iters * w_icp_step;
    % decrease w_dist
    %w_dist = w_dist - w_dist_step;
end

figure;plot(log10(E), '-x');title('Error');xlabel('iteration');ylabel('log(error)');

end

function D = makeDmatrix(d)
d = reshape(d, 3, 1);
D = [d, skewMatrix(d), zeros(3)];
end

function S = skewMatrix(v)
S = [0 v(3) -v(2); -v(3) 0 v(1); v(2) -v(1) 0];
end