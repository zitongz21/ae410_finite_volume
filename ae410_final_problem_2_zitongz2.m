function ae410_final_problem_2_zitongz2()
disp('Calculating Problem 2 (a)...'); problem_2a();
disp('Calculating Problem 2 (b)...'); problem_2b();

function problem_2a()
% W :: \rho, u, p   L_state :: p, T, u
W_L = [1.0; 0.0; 1.0;];   W_R = [0.125; 0.0; 0.1;];
L1 = 1; L4 = 1; R = 287; gamma = 1.4; T_t = 0.25;
output_opt.test_csv = 'test_1.csv';
output_opt.figure_name = 'Toro Test 1: ROE''s Method';
output_opt.RP_flag = 1; output_opt.plot_flag = 1;
Nt = [100,500]; Nx = [100,500];
toro_test_zitongz2(W_L, W_R, T_t, L1, L4, R, gamma, Nt, Nx, 'theBrit', ...
    output_opt);

% calculate norm and plot rate of convergence
N = [500, 1000, 2500, 5000, 10000];
output_opt.plot_flag = 0;
[L2, Ln] = toro_test_zitongz2(W_L, W_R, T_t, L1, L4, R, gamma, N, N, ...
    'theBrit', output_opt); % calculate norm

figure('Color', [1,1,1]); % plot rate of convergence
for p = 1:3
    ax(p) = subplot(1,3,p);
    hold(ax(p), 'on'); grid(ax(p), 'on'); box(ax(p), 'on');
    set(ax(p), 'FontSize', 16, 'FontName', 'Times New Roman', ...
        'xscale', 'log', 'yscale', 'log');
    plot(ax(p), 1./N, L2(p,:), 'LineWidth', 1.2);
    plot(ax(p), 1./N, Ln(p,:), 'LineWidth', 1.2);
    plot(ax(p), 1./N, 1./N, 'LineWidth', 1.6, 'Color', [0,0,0]);
    legend(ax(p), {'$L-2$ norm', '$L-\infty$ norm', ...
        '$1^{st}$-order accurate'}, 'FontSize',14, 'interpreter','latex');
end


function problem_2b()
% Problem 2 (b)
W_L = [1.0;-2.0; 0.4]; W_R = [1.0; 2.0; 0.4;]; 
L1 = 1; L4 = 1; R = 287; gamma = 1.4; T_t = 0.15; Nt = 10000; Nx = 500;
output_opt.test_csv = 'test_2.csv';
output_opt.figure_name = 'Toro Test 2: ROE''s Method';
output_opt.RP_flag = 1; output_opt.plot_flag = 1;
% This too shall fail
toro_test_zitongz2(W_L, W_R, T_t, L1, L4, R, gamma, Nt, Nx, 'theBrit', ...
    output_opt); 

% Take it easy, mate.

output_opt.figure_name = 'Toro Test 2 Modified: ROE''s Method';
W_L = [1.0;-0.5; 0.4]; W_R = [1.0; 0.5; 0.4;]; T_t = 0.15;
Nt = [100,500]; Nx = [100,500];	output_opt.test_csv = '';
toro_test_zitongz2(W_L, W_R, T_t, L1, L4, R, gamma, Nt, Nx, 'theBrit', ...
    output_opt);

% calculate norm and plot rate of convergence
N = [500, 1000, 2500, 5000, 10000];
output_opt.plot_flag = 0;
[L2, Ln] = toro_test_zitongz2(W_L, W_R, T_t, L1, L4, R, gamma, N, N, ...
    'theBrit', output_opt); % calculate norm

figure('Color', [1,1,1]); % plot rate of convergence
for p = 1:3
    ax(p) = subplot(1,3,p);
    hold(ax(p), 'on'); grid(ax(p), 'on'); box(ax(p), 'on');
    set(ax(p), 'FontSize', 16, 'FontName', 'Times New Roman', ...
        'xscale', 'log', 'yscale', 'log');
    plot(ax(p), 1./N, L2(p,:), 'LineWidth', 1.2);
    plot(ax(p), 1./N, Ln(p,:), 'LineWidth', 1.2);
    plot(ax(p), 1./N, 1./N, 'LineWidth', 1.6, 'Color', [0,0,0]);
    legend(ax(p), {'$L-2$ norm', '$L-\infty$ norm', ...
        '$1^{st}$-order accurate'}, 'FontSize',14, 'interpreter','latex');
end