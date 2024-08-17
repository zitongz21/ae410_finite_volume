clear all; clc;
ylabel_str = {'$\rho$', '$u$', '$p$', '$e$'};
% W_L = [1.0; 0.0; 100.0;];   W_R = [1; 0.0; 0.1;];

L1 = 19;   L4 = 1;   R = 287;    T = 293;
gamma = 1.4;    Nx = 2000;
T_t = 0.10;     Nt = 1000;
L_state = [1e7, 293, 0]; R_state = [1e4, 293, 0];
BC = {'wall', 'wall'; 0, 0;}; method = 'theBrit';


W_L = [1e7/R/T; 0; 1e7;];	W_R =[1e4/R/T; 0; 1e4;];
X = linspace(-1,19, Nx+1);  T = linspace(0, T_t, Nt+1)';
D = zeros(Nt+1, Nx+1); U = zeros(Nt+1, Nx+1); P = zeros(Nt+1, Nx+1);
D(1, X< 0) = W_L(1); U(1, X< 0) = W_L(2); P(1, X< 0) = W_L(3);
D(1, X>=0) = W_R(1); U(1, X>=0) = W_R(2); P(1, X>=0) = W_R(3);

figure('Color', [1,1,1]);
for i = 1:Nt
    [W, ~, ~] = ...
        RP_Euler_Solver_ZZ(W_L, W_R, gamma, X, T(i+1));
    for j = 1:3; subplot(3,1,j); plot(W(j,:)); end; drawnow
    D(i+1,:) = W(1,:); U(i+1,:) = W(2,:); P(i+1,:) = W(3,:);
end

X_mat = ones(Nt+1,1) * X; T_mat = T * ones(1, Nx+1);
fig = figure('Color', [1,1,1]);
subplot(3,1,1); contourf( X_mat, T_mat, log(D)); colormap jet;
subplot(3,1,2); contourf( X_mat, T_mat, U); colormap jet;
subplot(3,1,3); contourf( X_mat, T_mat, log(P)); colormap jet;