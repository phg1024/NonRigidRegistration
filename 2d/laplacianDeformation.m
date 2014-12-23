function pp = laplacianDeformation(p, q, pivots, mode)

% minimize E_fitting + E_distortion

npts = size(p, 1);

% compute laplacians
delta = zeros(npts, 2);
delta(1,:) = p(1,:) - p(2,:);
for i=2:npts-1
    delta(i,:) = p(i,:) - 0.5*(p(i-1,:) + p(i+1,:));
end
delta(npts,:) = p(npts,:) - p(npts-1,:);

C = cell(npts, 1);
C{1} = [makeCMatrix(p(1,:)); makeCMatrix(p(2,:))];
for i=2:npts-1
    C{i} = [makeCMatrix(p(i,:)); makeCMatrix(p(i-1,:)); makeCMatrix(p(i+1,:))];
end
C{npts} = [makeCMatrix(p(npts,:)); makeCMatrix(p(npts-1,:))];

D = cell(npts, 1);
for i=1:npts
    dx = delta(i,1); dy = delta(i,2);
    D{i} = [dx dy 0 0; dy -dx 0 0];
end

% use only a fraction of the data
if nargin < 3
    k = 9;
    indices = [1, setdiff(randperm(npts, k), [1, npts]), npts];
else
    indices = pivots;
end
k = length(indices);

% fittint term
Ai = zeros(1, k*2); Aj = zeros(1, k*2); Av = zeros(1, k*2);
idx = 1;
for i=1:k
    Ai(idx) = idx; Aj(idx) = (indices(i)-1)*2+1; Av(idx) = 1; idx = idx + 1;
    Ai(idx) = idx; Aj(idx) = indices(i)*2; Av(idx) = 1; idx = idx + 1;
end

% distortion term
Bi = zeros(1, 3*npts-2); Bj = zeros(1, 3*npts-2); Bv = zeros(1, 3*npts-2);
idx = 1;
for i=1:npts
    if i==1
        Bi(idx) = (i-1)*2+1; Bj(idx) = 1; Bv(idx) = 1; idx = idx + 1;
        Bi(idx) = (i-1)*2+1; Bj(idx) = 3; Bv(idx) = -1; idx = idx + 1;
        
        Bi(idx) = i*2; Bj(idx) = 2; Bv(idx) = 1; idx = idx + 1;
        Bi(idx) = i*2; Bj(idx) = 4; Bv(idx) = -1; idx = idx + 1;
    elseif i==npts
        Bi(idx) = (i-1)*2+1; Bj(idx) = 2*npts-3; Bv(idx) = -1; idx = idx + 1;
        Bi(idx) = (i-1)*2+1; Bj(idx) = 2*npts-1; Bv(idx) = 1; idx = idx + 1;
        
        Bi(idx) = i*2; Bj(idx) = 2*npts-2; Bv(idx) = -1; idx = idx + 1;
        Bi(idx) = i*2; Bj(idx) = 2*npts; Bv(idx) = 1; idx = idx + 1;
    else
        Bi(idx) = (i-1)*2+1; Bj(idx) = (i-2)*2+1; Bv(idx) = -0.5; idx = idx + 1;
        Bi(idx) = (i-1)*2+1; Bj(idx) = (i-1)*2+1; Bv(idx) = 1; idx = idx + 1;
        Bi(idx) = (i-1)*2+1; Bj(idx) = i*2+1; Bv(idx) = -0.5; idx = idx + 1;
        
        Bi(idx) = i*2; Bj(idx) = (i-1)*2; Bv(idx) = -0.5; idx = idx + 1;
        Bi(idx) = i*2; Bj(idx) = i*2; Bv(idx) = 1; idx = idx + 1;
        Bi(idx) = i*2; Bj(idx) = (i+1)*2; Bv(idx) = -0.5; idx = idx + 1;
    end
end

A = sparse(Ai, Aj, Av, k*2, npts*2);
B = sparse(Bi, Bj, Bv, npts*2, npts*2);

if nargin < 4
    mode = 1;
end

if mode == 1
    for i=1:npts
        T = D{i} * ((C{i}'*C{i})\C{i}');
        if i==1
            B(1:2,1:4) = B(1:2,1:4) - T;
        elseif i==npts
            B(end-1:end, end-1:end) = B(end-1:end, end-1:end) - T(:,1:2);
            B(end-1:end, end-3:end-2) = B(end-1:end, end-3:end-2) - T(:,3:4);
        else
            B((i-1)*2+1:i*2,(i-1)*2+1:i*2) = B((i-1)*2+1:i*2,(i-1)*2+1:i*2) - T(:,1:2);
            B((i-1)*2+1:i*2,(i-2)*2+1:(i-1)*2) = B((i-1)*2+1:i*2,(i-2)*2+1:(i-1)*2) - T(:,3:4);
            B((i-1)*2+1:i*2,i*2+1:(i+1)*2) = B((i-1)*2+1:i*2,i*2+1:(i+1)*2) - T(:,5:6);
        end
    end
end

M = [A; B];
if mode == 1
    b = [reshape(q(indices,:)',k*2,1); zeros(npts*2, 1)];
else
    b = [reshape(q(indices,:)', k*2,1); reshape(delta', 2*npts, 1)];
end

pp = (M'*M)\(M'*b);
pp = reshape(pp, 2, npts)';
end

function C = makeCMatrix(p)
x = p(1); y = p(2);
C = [x y 1 0; y -x 0 1];
end