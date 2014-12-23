% generate two sets of points
p = [ [0:0.1:10]', sin([0:0.1:10]*3)' ];

%q = [ [0:0.1:10]', 2.0 * exp(-(p(:,1)-4.0).^2/8.0) + p(:,2)];
q = [ [0:0.1:10]', p(:,1) + p(:,2)];

pivots = [1, setdiff(randperm(size(p,1), 3), [1, size(p,1)]), size(p,1)];
pp = laplacianDeformation(p, q, pivots,2);

figure; hold on;
plot(p(:,1), p(:,2), '-');
plot(q(:,1), q(:,2), '-');
plot(pp(:,1), pp(:,2), '-.');
plot(q(pivots,1), q(pivots,2), 'bs');

