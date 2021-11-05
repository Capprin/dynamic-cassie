%% optimize_so3 test on a reaction wheel satellite
% create samples in ea. shape dimension
num_samples=31; bound=pi/4;
dim = linspace(-bound,bound,num_samples);
samples = {dim,dim};
% create local connection
z = zeros(num_samples);
o = ones(num_samples);
% the following are valid positive configurations
A_orig = {z, z; -o, z; z, -o};
%A_orig = {-o z; z z; z -o};
%A_orig = {-o z; z -o; z z};
% optimize
[grid, X, Y, Z, A_opt] = optimize_so3(samples, A_orig);

%% compare lengths
X_orig = [A_orig{1,1}(:)';A_orig{1,2}(:)'];
Y_orig = [A_orig{2,1}(:)';A_orig{2,2}(:)'];
Z_orig = [A_orig{3,1}(:)';A_orig{3,2}(:)'];
n_orig = norm([mean(vecnorm(X_orig)) mean(vecnorm(Y_orig)) mean(vecnorm(Z_orig))]);
disp(['Original norm-average metric: ', num2str(n_orig)]);

X_opt = [A_opt{1,1}(:)';A_opt{1,2}(:)'];
Y_opt = [A_opt{2,1}(:)';A_opt{2,2}(:)'];
Z_opt = [A_opt{3,1}(:)';A_opt{3,2}(:)'];
n_opt = norm([mean(vecnorm(X_opt)) mean(vecnorm(Y_opt)) mean(vecnorm(Z_opt))]);
disp(['Optimal norm-average metric: ', num2str(n_opt)]);

%% plot outputs
figure(1);clf;

subplot(2,3,1);
quiver(grid{1},grid{2},A_orig{1,1},A_orig{1,2}, 'Color', [0 0 0]);
axis('square');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('A_x, original coordinates');
subplot(2,3,2);
quiver(grid{1},grid{2},A_orig{2,1},A_orig{2,2}, 'Color', [0 0 0]);
axis('square');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('A_y, original coordinates');
subplot(2,3,3);
quiver(grid{1},grid{2},A_orig{3,1},A_orig{3,2}, 'Color', [0 0 0]);
axis('square');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('A_z, original coordinates');

subplot(2,3,4);
quiver(grid{1},grid{2},A_opt{1,1},A_opt{1,2}, 'Color', [0 0 0]);
axis('square');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('A_x, optimal coordinates');
subplot(2,3,5);
quiver(grid{1},grid{2},A_opt{2,1},A_opt{2,2}, 'Color', [0 0 0]);
axis('square');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('A_y, optimal coordinates');
subplot(2,3,6);
quiver(grid{1},grid{2},A_opt{3,1},A_opt{3,2}, 'Color', [0 0 0]);
axis('square');
xlabel('\alpha_1'); ylabel('\alpha_2');
title('A_z, optimal coordinates');
