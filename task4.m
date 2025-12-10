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
m = [101, 201, 401, 601, 801];



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
u0 = [p0 v0];


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

Ll_d = kron(e1, e_1.');
Lr_d = kron(e1, e_m.');
L = [Ll_d;
    Lr_d];
P = eye(2*m) - H_bar\L.' / (L/H_bar * L.')*L;
M_d = -P*(D_x)*P;



[t, u] = RK4(M, u0, [0, 1.8], 0.05*h);
p = u(1:m, end);
v = u(m+1:2*m, end);


plot(x, p, 'r-', 'LineWidth', 1.5);
hold on
plot(x, v, 'b--', 'LineWidth', 1.5);
hold off

legend('p', 'v')
title('RK4 grejer')



%%
plot(x, p0, 'r-', 'LineWidth', 1.5);
hold on
plot(x, v0, 'b--', 'LineWidth', 1.5);
hold off

legend('p', 'v')
title('Initialvärden')

%% Tabellen

%analytic solution
p_analytic = @(x, t) theta_2(x, 2-t) - theta_1(x, 2-t);
v_analytic = @(x, t) theta_1(x, 2-t) + theta_2(x, 2-t);
u_analytic = @(x, t) [p_analytic(x,t); v_analytic(x,t)];

%m_values = [101, 201, 401, 601, 801];
m_values = [101];
err_norms = zeros(size(m_values));

for i = 1:length(m_values)
    m= m_values(i);    
    
    x = linspace(-1, 1, m);
    h = domain_width/(m - 1);

    % 1. Återskapa initialdata u0 för aktuellt m
    p0 = theta_1(x,0) - theta_2(x,0);
    v0 = theta_1(x,0) + theta_2(x,0);
    u0_current = [p0(:); v0(:)]; % Kolumnvektor (2m x 1)

    %UPWIND SBP
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
    
    Ll_d = kron(e1, e_1');
    Lr_d = kron(e1, e_m');
    L = [Ll_d;
        Lr_d];

    P = eye(2*m) - H_bar\L.' / (L/H_bar * L.')*L;
    M_d = -P*(D_x)*P;
    
    [~, u_matrix] = RK4(M_d, u0_current, [0, t_star], alpha*h);
    u_m = u_matrix(:, end);

    u_anal_matrix = u_analytic(x, t_star); %2 x m matrix
    u_anal = u_anal_matrix(:); % 2m x 1 

    %Beräkna error grejerna
    err = u_anal - u_m;
    err_norm = sqrt(h) * norm(err);
    %Choppa in i listan
    err_norms(i) = err_norm;
end   

%% Testar min analytic och det blir skit

p_anal_plot = u_anal(1:m, end);
v_anal_plot = u_anal(m+1:2*m, end);


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






