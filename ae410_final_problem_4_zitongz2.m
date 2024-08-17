function ae410_final_problem_4_zitongz2()
disp('Problem 4: Godunov''s Method');
disp('Calculating Toro Test 1...');         tic; %toro_test_1(); toc;
disp('Calculating Toro Test 2...');         tic; %toro_test_2(); toc;
disp('Calculating Shock Tube Problem...');	tic;
L1 = 19;   L4 = 1;   R = 287;    gamma = 1.4;    Nx = 1000;
T_t = 0.10;          Nt = 10000;
L_state = [1e7, 293, 0];    R_state = [1e4, 293, 0];
BC = {'wall', 'wall'; 0, 0;}; method = 'theRussian';
[rho, u, p, tmp] = shock_tube(L1, L4, R, gamma, Nx, Nt, T_t, ...
    L_state, R_state, BC, method);      toc;

disp('Plotting solution contour...');
X = linspace(0, L1+L4, Nx+1); X_c = 0.5 * ( X(1:end-1) + X(2:end) );
T = linspace(0, T_t, Nt+1)';
X_mat = ones(Nt+1,1) * X_c; T_mat = T * ones(1, Nx);
Z = {p, rho, tmp, u};

figure('Color', [1,1,1], 'Position', [1 41 1920 963], ...
    'Name', 'Solution Contour (Godunov''s Method)');
x_ind = X_c <= 20; t_ind = 1:10:Nt+1;
edge_color = {'none', 'none', 'none', 'none'};
c_scale = {'log', 'log', 'log', ''};
cmap = {'purple_seq', 'green_seq', 'black_body', 'cool_warm'};
title_cell = {'Pressure ($Pa$)', 'Density ($kg/m^3$)', ...
    'Temperature ($K$)', 'Velocity ($m/s$)'};
c_lim = {[], [], [], [-900, 900]};
for j = 1:4
    ax_a(j) = subplot(2,2,j); hold(ax_a(j), 'on'); box(ax_a(j), 'on');
    x_t_contour_zitongz2(ax_a(j), X_mat, T_mat, Z{j}, x_ind, t_ind, ...
            c_lim{j}, title_cell{j}, edge_color{j}, c_scale{j}, cmap{j});
    scatter(ax_a(j), 19.99, 0.0174, 's', 'filled', ...
        'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', [0,0,0]);
    scatter(ax_a(j), 18.51, 0.0209, 'o', 'filled', ...
        'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', [0,0,0]);
    scatter(ax_a(j), 0.005, 0.00299,'d', 'filled', ...
        'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', [0,0,0]);
end

[rho_roe, u_roe, p_roe, ~] = shock_tube(L1, L4, R, gamma, ...
    Nx, Nt, T_t, L_state, R_state, BC, 'theBrit');
T_roe = linspace(0, T_t, Nt+1)';
X_roe = linspace(0, L1+L4, Nx+1);
X_c_roe = 0.5 * ( X_roe(1:end-1)+X_roe(2:end) );

fig(1) = figure('Color', [1,1,1], 'Position', [1 41 1920 960], ...
    'Name', 'ROE - Godunov, t = 0.01:0.05');
fig(2) = figure('Color', [1,1,1], 'Position', [1 41 1920 960], ...
    'Name', 'ROE - Godunov, t = 0.05:0.10');
T_plot = 0.01:0.01:0.10; 
var_str = {'\rho', 'u', 'p'};

for it = 1:5;
    t_id{1} = 1000*it + 1;   t_id{2} = 1000*(it+5) + 1;
    t_str{1} = ['(t=', num2str(T_plot(it)),      's)'];
    t_str{2} = ['(t=', num2str(T_plot(it)+0.05), 's)'];
    for v = 1:3
        for f = 1:2
            ax = subplot(3,5,it+(v-1)*5, 'Parent', fig(f));
            hold(ax, 'on'); grid(ax, 'on');
            set(ax, 'FontName', 'Times New Roman', 'FontSize', 12);
            xlabel(ax, '$x (m)$', 'interpreter', 'latex');
            title(ax, ['$' var_str{v} t_str{f} '$'], 'FontSize', 14, ...
                'interpreter', 'latex');
            switch v
                case 1; y_roe = rho_roe(t_id{f},:); y = rho(t_id{f},:);
                case 2; y_roe = u_roe(t_id{f},  :); y = u(t_id{f},  :);
                case 3; y_roe = p_roe(t_id{f},  :); y = p(t_id{f},  :);
            end
            plot(ax, X_c_roe, y_roe, 'LineWidth', 1.5);
            plot(ax, X_c, y, 'LineWidth', 1.5);
            legend(ax, {'R', 'G'}, 'FontSize', 12, ...
                'Location', 'best');
        end
    end
end


function toro_test_1()
W_L = [1.0; 0.0; 1.0;];   W_R = [0.125; 0.0; 0.1;];
L1 = 1; L4 = 1; R = 287; gamma = 1.4; T_t = 0.25;
output_opt.test_csv = 'test_1.csv';
output_opt.figure_name = 'Toro Test 1: Godunov''s Method';
output_opt.RP_flag = 1; output_opt.plot_flag = 1;
Nt = [100,500]; Nx = [100,500];
toro_test_zitongz2(W_L, W_R, T_t, L1, L4, R, gamma, Nt, Nx, ...
    'theRussian', output_opt);

function toro_test_2()
% Take it easy, mate.
W_L = [1.0;-0.5; 0.4]; W_R = [1.0; 0.5; 0.4;]; T_t = 0.15;
L1 = 1; L4 = 1; R = 287; gamma = 1.4; T_t = 0.25;
Nt = [100,500]; Nx = [100,500];	output_opt.test_csv = '';
output_opt.RP_flag = 1; output_opt.plot_flag = 1;
output_opt.figure_name = 'Toro Test 2 Modified: Godunov''s Method';
toro_test_zitongz2(W_L, W_R, T_t, L1, L4, R, gamma, Nt, Nx, ...
    'theRussian', output_opt);

function [rho, u, p, tmp] = shock_tube(L1, L4, R, gamma, Nx, Nt, T_t, ...
    L_state, R_state, BC, method)
    [q1, q2, q3] = shcok_tube_1d_solver(...
        L1, L4, R, gamma, Nx, T_t, Nt, L_state, R_state, method, BC);
    rho = q1; u = q2./q1;
    p = ( q3 - 0.5 * q2.^2 ./ q1) * (gamma-1);
    clear q1; clear q2; clear q3;
    tmp = p./rho/R;