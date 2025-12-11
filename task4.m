%Saker vi vet
rho = 1;
c = 1;
beta = 0;
r_star = 0.1;
x_l = -1;
x_r = 1;
t_star = 1.8;
alpha = 0.05;
L_domain = x_r - x_l;

%antal punkter
%m = [101, 201, 401, 601, 801];



%Gaussian profiles:
theta_1 = @(x, t) exp(-((x-t)/r_star).^2);
theta_2 = @(x, t) -exp(-((x+t)/r_star).^2);


%initial data
m = 101;

domain_width = 2;   %Detta är x_r - x_l
h = domain_width/(m - 1);

x = linspace(-1, 1, m);
p0 = theta_1(x,0) - theta_2(x,0);
v0 = theta_1(x,0) + theta_2(x,0);

%u0 = [p0; v0];
u0 = [p0(:); v0(:)]; % Kolumnvektor (2m x 1)


[H, HI, Dp, Dm, ~, ~] = SBP7_Upwind(m, h);
H_bar = kron(eye(2), H);  %Kan strunta i C eftersom den blir ehnhetsmatris (2m x 2m) eftersom roh = c =1
D_x = [zeros(m), Dp;
    Dm, zeros(m)];     

% Dirichlet
e_1 = zeros(m,1);
e_1(1) = 1;
e_m = zeros(m, 1);
e_m(m) = 1;
e1 = [1, 0];    % Detta är e^(1)
e2 = [0, 1];    % e^(2)

Ll_d = kron(e1, e_1.');     %Pressure är 0 dirichlet
Lr_d = kron(e1, e_m.');
L = [Ll_d;
    Lr_d];
P = eye(2*m) - H_bar\L.' / (L/H_bar * L.')*L;
M_d = -P*(D_x)*P;



[t, u] = RK4(M_d, u0, [0, 1.8], 0.05*h);
p = u(1:m, end);
v = u(m+1:2*m, end);


%% Plot vid t=1.8
plot(x, p, 'r-', 'LineWidth', 1.5);
hold on
plot(x, v, 'b--', 'LineWidth', 1.5);
hold off

legend('p', 'v')
title('Numerical solution at t=1.8')



%% Plot vid t=0
plot(x, p0, 'r-', 'LineWidth', 1.5);
hold on
plot(x, v0, 'b--', 'LineWidth', 1.5);
hold off

legend('p', 'v')
title('Numerical solution at t=0')

%% Error SBP7
m_values = [101, 201, 401, 601, 801];
err_norms = zeros(size(m_values));

for i = 1:length(m_values)
    
    m= m_values(i);    
    x = linspace(-1, 1, m);
    h = domain_width/(m - 1);
    
    %analytic solution vid t_star
    p_analytic = theta_2(x, 2-t_star) - theta_1(x, 2-t_star);
    v_analytic = theta_1(x, 2-t_star) + theta_2(x, 2-t_star);
    u_analytic = [p_analytic(:); v_analytic(:)];

    % 1. Återskapa initialdata u0 för aktuellt m
    p0 = theta_1(x,0) - theta_2(x,0);
    v0 = theta_1(x,0) + theta_2(x,0);
    u0_current = [p0(:); v0(:)]; % Kolumnvektor (2m x 1)

    %UPWIND SBP
    [H, HI, Dp, Dm, ~, ~] = SBP7_Upwind(m, h);
    H_bar = kron(eye(2), H);
    D_x = [zeros(m), Dp;
        Dm, zeros(m)];     
    
    % Dirichlet
    e_1 = zeros(m,1);
    e_1(1) = 1;
    e_m = zeros(m, 1);
    e_m(m) = 1;
    e1 = [1, 0];    % Detta är e^(1)
    e2 = [0, 1];    % e^(2)
    
    Ll_d = kron(e1, e_1.');
    Lr_d = kron(e1, e_m.');
    L = [Ll_d;
        Lr_d];

    P = eye(2*m) - H_bar\L.' / (L/H_bar * L.')*L;
    M_d = -P*(D_x)*P;
    
    [~, u_matrix] = RK4(M_d, u0_current, [0, t_star], alpha*h);
    u_m = u_matrix(:, end);

    %Beräkna error grejerna
    err = u_analytic - u_m;
    err_norm = sqrt(h) * norm(err);
    %Choppa in i listan
    err_norms_SBP7(i) = err_norm;
end

%% q för SBP7
q_values_SBP7 = zeros(size(m_values));

for i = 2:length(m_values)
    m_prev = m_values(i-1);
    m_now = m_values(i);
    err_prev = err_norms_SBP7(i-1);
    err_now = err_norms_SBP7(i);

    q = log10(err_prev/err_now) / log10(m_now/m_prev);
    q_values_SBP7(i) = q;
end

%% Testar min analytic och det blir jättebra

p_anal_plot = u_analytic(1:m, end);
v_anal_plot = u_analytic(m+1:2*m, end);


plot(x, p_anal_plot, 'r-', 'LineWidth', 1.5);
hold on
plot(x, v_anal_plot, 'b--', 'LineWidth', 1.5);
hold off

legend('p', 'v')
title('Analytic')

%error för m punkter
%err_m = u_t_star - u_m;
%err_m_norm = sqrt(h) * norm(err_m);

%Convergence rate (q):
%q = log10(err_m_norm/err_n_norm) / log10(n/m);


%% Error SBP 6

m_values = [101, 201, 401, 601, 801];
err_norms = zeros(size(m_values));

for i = 1:length(m_values)
    
    m= m_values(i);    
    x = linspace(-1, 1, m);
    h = domain_width/(m - 1);
    
    %analytic solution vid t_star
    p_analytic = theta_2(x, 2-t_star) - theta_1(x, 2-t_star);
    v_analytic = theta_1(x, 2-t_star) + theta_2(x, 2-t_star);
    u_analytic = [p_analytic(:); v_analytic(:)];

    % 1. Återskapa initialdata u0 för aktuellt m
    p0 = theta_1(x,0) - theta_2(x,0);
    v0 = theta_1(x,0) + theta_2(x,0);
    u0_current = [p0(:); v0(:)]; % Kolumnvektor (2m x 1)

    %UPWIND SBP
    [H, HI, D1, ~, ~, ~] = SBP6(m, h);

    Dp = D1;
    Dm = D1; 

    H_bar = kron(eye(2), H);
    D_x = [zeros(m), Dp;
        Dm, zeros(m)];     
    
    % Dirichlet
    e_1 = zeros(m,1);
    e_1(1) = 1;
    e_m = zeros(m, 1);
    e_m(m) = 1;
    e1 = [1, 0];    % Detta är e^(1)
    e2 = [0, 1];    % e^(2)
    
    Ll_d = kron(e1, e_1.');
    Lr_d = kron(e1, e_m.');
    L = [Ll_d;
        Lr_d];

    P = eye(2*m) - H_bar\L.' / (L/H_bar * L.')*L;
    M_d = -P*(D_x)*P;
    
    [~, u_matrix] = RK4(M_d, u0_current, [0, t_star], alpha*h);
    u_m = u_matrix(:, end);

    %Beräkna error grejerna
    err = u_analytic - u_m;
    err_norm = sqrt(h) * norm(err);
    %Choppa in i listan
    err_norms_SBP6(i) = err_norm;
end   

%% q för SBP6

q_values_SBP6 = zeros(size(m_values));

for i = 2:length(m_values)
    m_prev = m_values(i-1);
    m_now = m_values(i);
    err_prev = err_norms_SBP6(i-1);
    err_now = err_norms_SBP6(i);

    q = log10(err_prev/err_now) / log10(m_now/m_prev);
    q_values_SBP6(i) = q;
end