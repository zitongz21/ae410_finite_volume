function ae410_final_problem_3_zitongz2()
L1 = 19;   L4 = 1;   R = 287;    gamma = 1.4;    Nx = 2000;
T_t = 0.10;          Nt = 20000;
L_state = [1e7, 293, 0]; R_state = [1e4, 293, 0];
BC = {'wall', 'wall'; 0, 0;}; method = 'theBrit';

disp('Solving shock tube problem with ROE''s method...');
[rho, u, p, tmp] = shock_tube(L1, L4, R, gamma, Nx, Nt, T_t, ...
    L_state, R_state, BC, method);
X = linspace(0, L1+L4, Nx+1); X_c = 0.5 * ( X(1:end-1) + X(2:end) );
T = linspace(0, T_t, Nt+1)';
X_mat = ones(Nt+1,1) * X_c; T_mat = T * ones(1, Nx);
Z = {p, rho, tmp, u};

%% part (a)
disp('Problem 3 (a) in progress. Plotting solution contour...');
figure('Color', [1,1,1], 'Position', [1 41 1920 963], ...
    'Name', 'Solution Contour');
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
end

figure('Color', [1,1,1], 'Position', [5 532 1911 446], ...
    'Name', 'Solution Contour Zoomed in');
% expansion facn region, linear pressure, density and temperature;
x_ind = X_c <= 1; t_ind = 1:2:Nt/5+1;
c_lim = {[], [], [], [0, 1]};
cmap = {'purple_seq', 'green_seq', 'black_body', 'warm'};
for j = 1:4
    ax_b(j) = subplot(1,4,j); hold(ax_b(j), 'on'); box(ax_b(j), 'on');
    x_t_contour_zitongz2(ax_b(j), X_mat, T_mat, Z{j}, x_ind, t_ind, ...
            c_lim{j}, title_cell{j}, edge_color{j}, '', cmap{j});
    scatter(ax_b(j), 0.005, 0.00299, 's', 'filled', ...
        'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', [0,0,0]);
end
disp(['Problem 3 (a) Complete. ' ...
    'Press any key to clear all plots and continue.']); pause();

%% part (b)
disp('Problem 3 (b) in progress. Calculating primary shock speed...');
X_shock = ones(Nt+1,1);
for i = 2:Nt+1
    [~, i_shock] = min(diff(u(i,:)));
    X_shock(i) = X(i_shock+1);
end
figure('Color', [1,1,1], 'Position', [245 395 1600 600], ...
    'Name', 'Identifying and Locating the Primary Shock');
for i = 1:2; 
    ax(i) = subplot(1,2,i); hold(ax(i), 'on'); grid(ax(i), 'on');
    set(ax(i), 'FontSize', 14, 'FontName', 'Times New Roman');
end
plot(ax(1), X_c, u(1001,:), 'LineWidth', 1.8); 
plot(ax(1), X(2:end-1), diff(u(1001,:)), 'LineWidth', 1.5); 
title(ax(1), ['$t = ', num2str(T_t/Nt*1000), 's$'], ...
    'FontSize', 18, 'interpreter', 'latex');
legend(ax(1), {'$u$', '$\Delta u$'}, ...
    'FontSize', 18, 'interpreter', 'latex');
xlabel(ax(1), '$X (m)$', 'interpreter', 'latex'); 
scatter(ax(2), X_shock, T, 16);
title(ax(2), 'Shock Location', 'FontSize', 16);
xlabel(ax(2), '$X (m)$', 'interpreter', 'latex'); 
ylabel(ax(2), '$t (s)$', 'interpreter', 'latex');
t_gating = T >= 5e-5 & T<= 0.016;
s_shock = polyfit(T(t_gating), X_shock(t_gating), 1); s_shock = s_shock(1);
x_fit = 1:20; t_fit = (x_fit-1)./s_shock;
plot(ax(2), x_fit, t_fit, 'LineWidth', 2.4);
text(9, mean(t_fit), ['$ t = (x-1)/', num2str(s_shock) ...
    '$'], 'FontSize', 16,  'interpreter', 'latex', 'Parent', ax(2), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
legend(ax(2), {'Shock Location', 'Linear Fit of Primary shock'}, ...
    'FontSize', 14);
xlabel(ax(2), '$X (m)$', 'interpreter', 'latex'); 
ylabel(ax(2), '$t (s)$', 'interpreter', 'latex');

clear('rho', 'u', 'p', 'tmp');

fun = @(M_s) shock_tube_equation(M_s, gamma, R, L_state(2), R_state(2), ...
    L_state(1), R_state(1));
M_s = fzero(fun, 1);    c4 = sqrt(gamma*R*L_state(2));
T_t = 0.02; 
N_X = [250, 500, 1000, 2000, 4000, 6000]; shock_speed = zeros(1,6);
for n = 1:6
    Nt = 2*N_X(n);
    [rho, u, ~, ~] = shock_tube(L1, L4, R, gamma, N_X(n), Nt, T_t, ...
        L_state, R_state, BC, method);
    X = linspace(0, L1+L4, N_X(n)+1);   X_shock = ones(Nt+1,1);
    for i = 2:Nt+1; 
        [~, i_shock] = min(diff(u(i,:))); X_shock(i) = X(i_shock+1);
    end
    T = linspace(0, T_t, Nt+1)';      t_gating = T >= 5e-5 & T<= 0.016;  
    ss = polyfit(T(t_gating), X_shock(t_gating), 1);
    shock_speed(n) = ss(1); clear('u');
end
ss_extrap = interp1((L1+L4)./N_X, shock_speed./c4, 0,'pchip');

disp([shock_speed./c4,ss_extrap]);
figure('Color', [1,1,1], 'Position', [100 260 720 600], ...
    'Name', 'Convergence of Primary Shock Mach Number');
ax = axes('FontSize', 14, 'FontName', 'Times New Roman');
grid(ax, 'on'); hold(ax, 'on');
plot(ax, [(L1+L4)./N_X,0], M_s * ones(1,7), 'LineWidth', 1.2);
line_c = [222,45,38]/255;
plot((L1+L4)./N_X, shock_speed./c4, '-o', 'LineWidth',1.6, 'Color',line_c);
plot([(L1+L4)./N_X(end), 0], [shock_speed(end)./c4, ss_extrap], ...
    'LineWidth', 1.6, 'LineStyle', ':', 'Color', line_c);
title(ax, 'Convergence of Primary Shock Mach Number', 'FontSize', 16);
ylim(ax, [3.14, 3.24]);
xlabel(ax, '$\Delta x\ (m)$', 'FontSize', 16, 'interpreter', 'latex');
ylabel(ax, '$M_s$', 'FontSize', 16, 'interpreter', 'latex');
legend(ax ,{'Theoretical','Solution','Extrapolated Solution'}, ...
    'FontSize',14, 'location', 'northwest');
disp(['Problem 3 (b) Complete. ' ...
    'Press any key to clear all plots and continue.']); pause();

%% part (c)
disp('Problem 3 (c) in progress. Calculating test time...');
rho_model = rho(:,X==L4+L1/2); diff_rho = diff(rho_model);
T_c = 0.5*(T(1:end-1) + T(2:end));
diff_gating = diff_rho > 0.0001;
t_shock = max(T_c(diff_gating & T_c < 0.01));
t_disc  = min(T_c(diff_gating & T_c > 0.01));
figure('Color', [1,1,1], 'Name', 'Test Time');
for i = 1:2; 
    ax(i) = subplot(1,2,i); grid(ax(i), 'on'); hold(ax(i), 'on');
    set(ax(i), 'FontSize', 14, 'FontName', 'Times New Roman');
    xlabel(ax(i), '$t(s)$', 'FontSize', 16, 'interpreter', 'latex');
end
test_gating_1 = T>=t_shock & T<=t_disc;
test_gating_2 = T>=t_shock & T<=t_disc;
legend_str = {'$t\in[0,0.02]$', 'Test Time'};
plot(ax(1), T, rho_model, 'LineWidth', 1.2);
plot(ax(1), T(test_gating_1), rho_model(test_gating_1), 'LineWidth', 2.4);
title(ax(1), 'Density at Test Model''s Location', 'FontSize', 16);
ylabel(ax(1), '$\rho_m (kg/m^3)$', 'FontSize', 16, 'interpreter', 'latex');
legend(ax(1), legend_str, 'FontSize', 14, 'interpreter', 'latex', ...
    'location', 'northwest');
plot(ax(2), T_c, diff_rho, 'LineWidth', 1.2);
plot(ax(2), T_c(test_gating_2), diff_rho(test_gating_2), 'LineWidth', 2.4);
title(ax(2), 'Density Change per Time Step at Test Model''s Location', ...
    'FontSize', 16);
ylabel(ax(2), '$\Delta\rho_m (kg/m^3)$', 'FontSize', 16, ...
    'interpreter', 'latex');
legend(ax(2), legend_str, 'FontSize', 14, 'interpreter', 'latex', ...
    'location', 'northwest');
disp({  't_shock',      num2str(t_shock, '%.6f'); ...
        't_contact',    num2str(t_disc, '%.6f'); ...
        'Test time',    num2str(t_disc - t_shock, '%.6f');});
disp(['Problem 3 (b) Complete. ' ...
    'Press any key to clear all plots and continue.']); pause();

%% shock tube
function [rho, u, p, tmp] = shock_tube(L1, L4, R, gamma, Nx, Nt, T_t, ...
    L_state, R_state, BC, method)
    [q1, q2, q3] = shcok_tube_1d_solver(...
        L1, L4, R, gamma, Nx, T_t, Nt, L_state, R_state, method, BC);
    rho = q1; u = q2./q1;
    p = ( q3 - 0.5 * q2.^2 ./ q1) * (gamma-1);
    clear q1; clear q2; clear q3;
    tmp = p./rho/R;

%% Primary shock Mach number
function epsilon = shock_tube_equation(M_s, gamma, R, T4, T1, p4, p1)
    c1 = sqrt(gamma*R*T1);  c4 = sqrt(gamma*R*T4);
    g1 = gamma - 1;     g2 = gamma + 1;
    epsilon = ( 2 * gamma * M_s.^2 - g1 ) / g2 * ...
        ( 1 - g1/g2*c1/c4*(M_s - 1./M_s)).^(-2*gamma/g1) - p4/p1;