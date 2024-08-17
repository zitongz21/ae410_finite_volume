function [q1, q2, q3, rho_t, u_t, p_t] = shcok_tube_1d_solver(...
    L1, L4, R, gamma, Nx, T_t, Nt, L_state, R_state, method, BC)
%%% -- input -- %%%
% L1      :: driven section length
% L4      :: driver section length
% R       :: gas constant
% gamma   :: specific heat ratio
% Nx      :: number of cells
% T_t     :: final physical time
% Nt      :: number of time-steps
% L_state :: driver section parameters, {p_L, T_L, u_L};
% R_state :: driven section parameters, {p_R, T_R, u_R};
% method  :: method of calculating numerical fluxes
%       Godunov's method - 'godunov' or 'theRussian'
%       Roe's method     - 'roe' or 'theBrit'
% BC      :: boundary condition type and wall velocity
%%% -- output -- %%%
% q1      :: matrix of conservd variable, density
% q2      :: matrix of conservd variable, momentum
% q3      :: matrix of conservd variable, energy
% rho_t   :: density at final time
% u_t     :: velocity at final time
% p_t     :: pressure at final time

g1 = gamma - 1;
p0_L = L_state(1);                      p0_R = R_state(1);
u0_L = L_state(3);                      u0_R = R_state(3);
rho0_L = p0_L/( R*L_state(2) );         rho0_R = p0_R/( R*R_state(2) );
rhoE0_L = p0_L/g1 + 0.5*rho0_L*u0_L^2;  
rhoE0_R = p0_R/g1 + 0.5*rho0_R*u0_R^2;

dt = T_t / Nt;
X = linspace(0, L1 + L4, Nx + 1);   dx = (L1+L4)/Nx;
X_c = 0.5 * ( X(1:end-1) + X(2:end) );

q1 = zeros(Nt+1, Nx);	q1(1,:) = ones(1, Nx) * rho0_R;
q2 = zeros(Nt+1, Nx);	q2(1,:) = ones(1, Nx) * rho0_R * u0_R;
q3 = zeros(Nt+1, Nx);	q3(1,:) = ones(1, Nx) * rhoE0_R;
l_gating = X_c < L4;
q1(1,l_gating) = rho0_L;
q2(1,l_gating) = rho0_L * u0_L;
q3(1,l_gating) = rhoE0_L;

% W_L_rp = [rho0_L; u0_L; p0_L]; W_R_rp = [rho0_R; u0_R; p0_R];
% X_RP = linspace(-1,1,21); O_RP = zeros(1,21);
% monitor_flag = 0;

epsilon_roe = 1e-2; dt_o_dx = dt / dx;
for it = 1:Nt
    if strcmpi(method, 'roe') || strcmpi(method, 'theBrit')
        F = roe_flux([q1(it,:);q2(it,:);q3(it,:);], ...
            gamma, epsilon_roe, BC);
    elseif strcmpi(method, 'godunov') || strcmpi(method, 'theRussian')
        if mod(it, Nt/100) == 0; disp([num2str(it*100/Nt), '%']); end
        F = godunov_flux([q1(it,:);q2(it,:);q3(it,:);], dx, dt, gamma, BC);
    end
    
    cfl_num = max(q2(it,:)./q1(it,:)) * dt_o_dx;
    if cfl_num > 0.5
        disp(cfl_num);
    end
    
    q1(it+1,:) = q1(it,:) - dt_o_dx * diff(F(1,:));
    q2(it+1,:) = q2(it,:) - dt_o_dx * diff(F(2,:));
    q3(it+1,:) = q3(it,:) - dt_o_dx * diff(F(3,:));
end

[rho_t, u_t, p_t] = q_to_w([q1(end,:); q2(end,:); q3(end,:)], g1);

function F = roe_flux(q, gamma, epsilon_roe, BC)
g1 = gamma - 1;
[rho, u, p] = q_to_w(q, g1);
rho_ghost = [rho(1), rho, rho(end)];    delta_rho = diff(rho_ghost);
p_ghost =   [  p(1), p,     p(end)];    delta_p   = diff(p_ghost);
if strcmpi(BC{1,1}, 'reflective') || strcmpi(BC{1,1}, 'wall')
    u_ghost =   [ -u(1) + 2*BC{2,1}, u];
else
    u_ghost =   [  u(1),             u];
end
if strcmpi(BC{1,2}, 'reflective') || strcmpi(BC{1,2}, 'wall')
    u_ghost =   [ u_ghost, -u(end) + 2*BC{2,2} ];
else
    u_ghost =   [ u_ghost,  u(end)             ];
end
delta_u   = diff(u_ghost);

[~, q2, q3] = w_to_q([rho_ghost; u_ghost; p_ghost], g1);
rho_roe = ( rho_ghost(1:end-1) .* rho_ghost(2:end) ).^0.5;
roe_denom = rho_ghost(1:end-1).^0.5 + rho_ghost(2:end).^0.5;
u_roe = (   rho_ghost(1:end-1).^0.5 .* u_ghost(1:end-1) + ...
            rho_ghost(2:end)  .^0.5 .* u_ghost(2:end)   ) ./ roe_denom;

H_ghost = ( q3 + p_ghost )./rho_ghost;
H_roe = (   rho_ghost(1:end-1).^0.5 .* H_ghost(1:end-1) + ...
            rho_ghost(2:end)  .^0.5 .* H_ghost(2:end)   ) ./ roe_denom;
c_roe = sqrt( g1 * ( H_roe - 0.5 * u_roe.^2 ) );
if ~isreal(c_roe)
    disp('Imag. Speed of sound.');
end
I = ones(1,size(rho,2)+1);
lam{1} = u_roe - c_roe;     v{1} = [I; lam{1}; H_roe - u_roe.*c_roe];
lam{2} = u_roe;             v{2} = [I; lam{2};   0.5 * u_roe.*u_roe];
lam{3} = u_roe + c_roe;     v{3} = [I; lam{3}; H_roe + u_roe.*c_roe];
alfa{2} = delta_rho - delta_p./c_roe.^2;

for i = 1:2:3 % entropy fix
    lam_fix = 0.5 * (lam{i}.^2/epsilon_roe + epsilon_roe);
    fix_gating = abs(lam{i}) < epsilon_roe;
    lam{i}(fix_gating) = lam_fix(fix_gating);
    alfa{i} = (delta_p + (i-2)*c_roe.*rho_roe.*delta_u)./(2*c_roe.^2);
end

F(1,:) = 0.5 * ( q2(1,1:end-1) + q2(1,2:end) );
F(2,:) = 0.5 * ( ...
    rho_ghost(1,1:end-1).*u_ghost(1,1:end-1).^2 + p_ghost(1,1:end-1) + ...
    rho_ghost(1,2:end)  .*u_ghost(1,2:end).^2   + p_ghost(1,2:end) );
F(3,:) = 0.5 * ( ...
    ( q3(1,1:end-1) + p_ghost(1,1:end-1) ) .* u_ghost(1,1:end-1) + ...
    ( q3(1,2:end)   + p_ghost(1,2:end)   ) .* u_ghost(1,2:end)  );
for i = 1:3
    for p = 1:3
        F(i,:) = F(i,:) - 0.5 * abs(lam{p}) .* alfa{p} .* v{p}(i,:);
    end
end

function F = godunov_flux(q, dx, dt, gamma, BC, X, X_c)
g1 = gamma - 1;
[RHO, U, P] = q_to_w(q, g1);
if min(P)< 0
    disp('!');
end
rho_ghost = [RHO(1), RHO, RHO(end)];
p_ghost =   [  P(1), P,     P(end)];

if strcmpi(BC{1,1}, 'reflective') || strcmpi(BC{1,1}, 'wall')
    u_ghost =   [ -U(1) + 2*BC{2,1}, U];
else
    u_ghost =   [  U(1),             U];
end
if strcmpi(BC{1,2}, 'reflective') || strcmpi(BC{1,2}, 'wall')
    u_ghost =   [ u_ghost, -U(end) + 2*BC{2,2} ];
else
    u_ghost =   [ u_ghost,  U(end)             ];
end
Nx_c = size(q,2);
F = zeros(3,Nx_c+1); U_debug = zeros(1,Nx_c+1);
for i = 1:Nx_c + 1
    inner_Nx = 9;  i_mid = 5;
    x = linspace(-dx, dx, inner_Nx);
    % w :: rho, u, p
    w_l = [rho_ghost(i);   u_ghost(i);   p_ghost(i)  ];
    w_r = [rho_ghost(i+1); u_ghost(i+1); p_ghost(i+1)];
    if w_l(1) == w_r(1) && w_l(2) == w_r(2) && w_l(3) == w_r(3)
        rho = rho_ghost(i); u = u_ghost(i); p = p_ghost(i);
    else
        w = RP_Euler_Solver_ZZ(w_l, w_r, gamma, x, dt/4);
        rho = w(1, i_mid); u = w(2, i_mid); p = w(3, i_mid);
    end
    U_debug(i) = u;
    F(1,i) = rho * u;
    F(2,i) = rho * u^2 + p;
    F(3,i) = (p/g1 + 0.5*rho*u^2 + p) * u;
end

function [rho, u, p] = q_to_w(q, g1)
% w :: [rho; u; p;]
% q :: [rho; rho*u; rho*E;]
rho =   q(1,:);
u   =   q(2,:) ./ q(1,:);
p   = ( q(3,:) - 0.5 * rho .* u.^2 ) * g1;

function [q1, q2, q3] = w_to_q(w, g1)
% q :: [rho; rho*u; rho*E;]
% w :: [rho; u; p;]
q1 = w(1,:);
q2 = w(2,:) .* w(1,:);
q3 = w(3,:)/g1 + 0.5 * w(1,:) .* w(2,:).^2;

%% monitor within main solver
% % initialize monitor plot
% if monitor_flag == 1
% hFig = figure('Color',[1,1,1]);
% title_str = {'q1', 'F_2', 'q2', 'u', 'q3', 'p'};
% for p = 1:6; 
%     ax(p)=subplot(3,2,p); title(ax(p), title_str{p}); hold(ax(p), 'on');
% end
% for i_var = 1:3
%     HP_Q(i_var) =    plot(ax(2*i_var-1), X_c,    q1(1,:));
%     HP_Q_rp(i_var) = plot(ax(2*i_var-1), X_RP+1, O_RP);
%     HP_W(i_var) =    plot(ax(2*i_var),   X_c,    q1(1,:));
%     HP_W_rp(i_var) = plot(ax(2*i_var),   X_RP+1, O_RP);
% end
% end
% 
% if monitor_flag == 2; 
%     subplot(2,1,1); h_flux = plot(X_c, X_c, 'linewidth', 1.0);
%     subplot(2,1,2); h_q = plot(X_c, X_c, 'linewidth', 1.0);
% end; drawnow;

%% monitor within Godunov flux
% if monitor_flag == 2;
%     set(h_flux, 'YData', diff(F(2,:))); 
%     set(h_q, 'YData', q2(it+1,:));
% end; drawnow;
%     
% if monitor_flag == 1    
%     % monitor \vec{q}
%     set(HP_Q(1), 'YData', q1(it+1,:)); 
%     set(HP_Q(2), 'YData', q2(it+1,:)); 
%     set(HP_Q(3), 'YData', q3(it+1,:));
%     
% %     % monitor W
%     [rho, u, p] = q_to_w([q1(it+1,:);q2(it+1,:);q3(it+1,:);], g1);
%     set(HP_W(1), 'YData', rho); 
%     set(HP_W(2), 'YData', u); 
%     set(HP_W(3), 'YData', p);
%     
%     set(hFig, 'Name', ['t = ' num2str(it*dt*1000) '(ms)']); 
%     
%     [W_RP, ~, ~] = ...
%         RP_Euler_Solver_ZZ(W_L_rp, W_R_rp, gamma, X_RP, (it+1)*dt);
%     %     % monitor W
%     set(HP_W_rp(1), 'YData', W_RP(1,:)); 
%     set(HP_W_rp(2), 'YData', W_RP(2,:)); 
%     set(HP_W_rp(3), 'YData', W_RP(3,:)); 
% 
%     % monitor q
%     [q1_rp, q2_rp, q3_rp] = w_to_q(W_RP(1:3,:), g1);
%     set(HP_Q_rp(1), 'YData', q1_rp); 
%     set(HP_Q_rp(2), 'YData', q2_rp); 
%     set(HP_Q_rp(3), 'YData', q3_rp);
%         
%     set(HP_W(1),    'YData', diff(F(2,:)));
%     set(HP_W_rp(1), 'YData', O_RP);
%     for p = 1:6; xlim(ax(p),[0,4]); end; drawnow;
% end