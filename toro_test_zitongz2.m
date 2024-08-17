function [L2, Ln] = toro_test_zitongz2(W_L, W_R, T_t, L1, L4, R, gamma, ...
    Nt, Nx, method, output_opt)
test_csv = output_opt.test_csv;
RP_flag = output_opt.RP_flag;
plot_flag = output_opt.plot_flag;

% W :: \rho, u, p   L_state :: p, T, u
L_state = [W_L(3), W_L(3)/R/W_L(1), W_L(2)];
R_state = [W_R(3), W_R(3)/R/W_R(1), W_R(2)];
BC = {'laissezfaire', 'laissezfaire'; 0, 0;};

c_fv = [222,45,38; 254,178,76; 49,163,84; 197,27,138;]/255; % plot color
n_fv = min([numel(Nx), numel(Nt)]); % number of finite element solutions
L2 = zeros(3, n_fv);        Ln = zeros(3, n_fv);            % norms
nx_str = 'N_x = ';          nt_str = 'N_t = ';
for i = 1:n_fv
    x{i} = linspace(-1, 1, Nx(i)+1);
    x_c{i} = 0.5 * ( x{i}(1:end-1) + x{i}(2:end) );
    W_FV{i} = zeros(3, Nx(i));
    [~, ~, ~, W_FV{i}(1,:), W_FV{i}(2,:), W_FV{i}(3,:)] = ...
        shcok_tube_1d_solver( L1, L4, R, gamma, Nx(i), T_t, Nt(i), ...
        L_state, R_state, method, BC);
    legend_str{i} = ['$FV_' num2str(i) '$'];
    nx_str = [nx_str, num2str(Nx(i)), ','];
    nt_str = [nt_str, num2str(Nt(i)), ','];
    if i>4; c_fv(i,:) = rand(1,3); end;
    % norm calculations
    W_rp{i} = RP_Euler_Solver_ZZ(W_L, W_R, gamma, x_c{i}, T_t);
    L2(:,i) = mean(   (W_rp{i}(1:3,:) - W_FV{i}).^2, 2 ) .^0.5;
    Ln(:,i) = max( abs(W_rp{i}(1:3,:) - W_FV{i}),[], 2 );
end
if n_fv == 1; legend_str{1} = 'FV'; end
% load posted solution
if ~isempty(test_csv); 
    test_BL = xlsread(test_csv);	legend_str{end+1} = '$ BL_1 $';
end

% exact RP solver
if RP_flag == 1
    x_rp = linspace(-1,1,101);      legend_str{end+1} = '$ BL_2 $';
    W = RP_Euler_Solver_ZZ(W_L, W_R, gamma, x_rp, T_t);
end

% plotting solution
if plot_flag == 1;
figure('Color', [1,1,1], 'Name', output_opt.figure_name);
title_str = {'Density', 'Velocity', 'Pressure'};
y_lab = {'$\rho$', '$u$', '$p$'};
for a = 1:3;
    ax(a) = subplot(3,1,a); 
    hold(ax(a), 'on'); grid(ax(a), 'on'); box(ax(a), 'on');
    set(ax(a), 'FontSize', 14, 'FontName', 'Times New Roman');
    for i = 1:n_fv
        plot(ax(a), x_c{i}, W_FV{i}(a,:), ...
            'LineWidth', 2.4-0.2*n_fv, 'Color', c_fv(i,:));
    end
    if ~isempty(test_csv)
        plot(test_BL(:,1), test_BL(:,a+1), 'LineWidth', 1.8, ...
            'Color', [0,0,0]/255);
    end
    if plot_flag == 1;
        plot(x_rp, W(a,:), 'LineWidth', 1.2, 'Color', [43,140,190]/255);
    end
    title(ax(a), [title_str{a} ' $(t=' num2str(T_t) ...
        ',\ ', nx_str(1:end-1) ',\ ' nt_str(1:end-1) ')$'], ...
        'FontSize',16, 'interpreter', 'latex');
    xlabel(ax(a), '$ X\ (\mathrm{m}) $', ...
        'FontSize', 14, 'interpreter', 'latex');
    ylabel(ax(a), y_lab{a}, 'FontSize', 14, 'interpreter', 'latex');
    legend(ax(a), legend_str, 'FontSize', 14, 'interpreter', 'latex');
end
end