%Saker vi vet
rho = 1;
c = 1;
beta = 0;
r_star = 0.1;
x_l = -1;
x_r = 1;
t_star = 1.8;
alpha = 0.05;
L = x_r - x_l;

%antal punkter
m = [101, 201, 401, 601, 801];


%error för m punkter
%err_m = u_t_star - u_m;
%err_m_norm = sqrt(h) * norm(err_m);

%error för n punkter
%err_n = u_t_star - u_n;
%err_n_norm = sqrt(h) * norm(err_n);


%Convergence rate (q):
%q = log10(err_m_norm/err_n_norm) / log10(n/m);


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




e_1 = zeros(m,1);
e_1(1) = 1;
e_m = zeros(m, 1);
e_m(m) = 1;
e1 = [1, 0];    % Detta är e^(1)
e2 = [0, 1];    % e^(2)

[H, HI, Dp, Dm, ~, ~] = SBP7_Upwind(m, h);

H_bar = kron(eye(2), H);  %Kan strunta i C eftersom den blir ehnhetsmatris (2m x 2m) eftersom roh = c =1

D_x = [zeros(m), Dp;
    Dm, zeros(m)];     

% Dirichlet

Ll_d = kron(e2, e_1');
Lr_d = kron(e2, e_m');
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




%analytic solution
p_analytic = @(x, t) theta_2(x, L-t) - theta_1(x, L-t);
v_analytic = @(x, t) theta_1(x, L-t) - theta_2(x, L-t);
u_analytic = @(x, t) [p_analytic; v_analytic];



%% Gemini med ett m-värde

% --- SBP-Projektionsmetoden: Körning för ett Enskilt Nätsteg (m) ---
clear; clc;

% =================================================================
% ➡️ PARAMETER SOM DU SKA ÄNDRA MELLAN KÖRNINGARNA
m = 101; 
% Testa m = 101, spara resultatet. 
% Byt till m = 201, spara resultatet, osv.
% =================================================================

% Saker vi vet (Problemparametrar)
rho = 1;        % Från text: rho=1
c = 1;          % Från text: c=1
beta = 0;       % Från text: beta=0
r_star = 0.1;   % Bredden på Gaussianerna
x_l = -1;       % Vänster rand
x_r = 1;        % Höger rand
T = 1.8;        % Sluttid
CFL = 0.05;
L = x_r - x_l; % Domänens längd (L=2)

% --- Analytiska Funktioner ---
% Gaussian profiles
theta_1 = @(x, t) exp(-((x-t)/r_star).^2);
theta_2 = @(x, t) -exp(-((x+t)/r_star).^2);

% Analytisk lösning u(x,t)
p_analytic = @(x, t) theta_2(x, L-t) - theta_1(x, L-t);
v_analytic = @(x, t) theta_1(x, L-t) - theta_2(x, L-t);
%analytic_solution = @(x, t) [p_analytic; v_analytic];
analytic_solution = @(x, t) theta_2(x - c*t, 0) - theta_1(x - c*t, 0);

% --- Nät och Tidstegsberäkning ---
h = L / (m - 1); 
x_grid = linspace(x_l, x_r, m)';
k_t = CFL * h / c; 
num_timesteps = ceil(T / k_t);
k_t = T / num_timesteps; % Justera k_t

fprintf('## SBP-Projektionsmetod: Enkel Körning ##\n');
fprintf('Nätstorlek m=%d (h=%.4f), dt=%.6f, N_t=%d\n', m, h, k_t, num_timesteps);
fprintf('-----------------------------------------------------\n');

% Initialisera lagring
Err_6th = 0;
Err_7th = 0;
u_exact = analytic_solution(x_grid, T);

% --- Lösare Logik (Körs för 6:e och 7:e ordningen) ---
orders = {'6th', '7th'};

for order_idx = 1:length(orders)
    order = orders{order_idx};
    
    % Steg 1: Skapa SBP-Operatorer (Anropar externa .m-filer)
    if strcmp(order, '6th')
        [D, H, HI] = SBP6_CD(m, h);
    else % '7th'
        [H, HI, Dp, Dm, ~, ~] = SBP7_Upwind(m, h);
    end
    
    % Steg 2: Skapa Projektionsmatriser och Randvärdesfunktion
    e_1 = zeros(m, 1); e_1(1) = 1;
    e_m = zeros(m, 1); e_m(m) = 1;
    E = e_1 * e_1' + e_m * e_m'; % Projektionsmatrisen E
    
    % Hämta Dirichlet-värden från den analytiska lösningen
    g_0 = @(t) analytic_solution(x_l, t);
    g_L = @(t) analytic_solution(x_r, t);
    
    % Randvärdesvektor g_vec(t)
    g_vec = @(t) g_0(t) * e_1 + g_L(t) * e_m;
    
    % Steg 3: Definiera Högerledet F(u, t)
    F = @(u, t) -c * D * u - invH * E * (u - g_vec(t));

    % Steg 4: Initialisera Lösningen u(x, 0)
    u_num = analytic_solution(x_grid, 0); 
    
    % Steg 5: Tidstegning med Klassisk RK4
    t = 0;
    for step = 1:num_timesteps
        k1 = F(u_num, t);
        k2 = F(u_num + k_t/2 * k1, t + k_t/2);
        k3 = F(u_num + k_t/2 * k2, t + k_t/2);
        k4 = F(u_num + k_t * k3, t + k_t);
        u_num = u_num + k_t/6 * (k1 + 2*k2 + 2*k3 + k4);
        t = t + k_t;
    end

    % Steg 6: Beräkna Fel och Skriv ut
    err = u_num - u_exact;
    Err_h = sqrt(h) * norm(err, 2);
    
    if strcmp(order, '6th')
        Err_6th = Err_h;
        fprintf('6:e ordningen (SBP-CD): ||err(m=%d)||_h = %.8e\n', m, Err_6th);
    else
        Err_7th = Err_h;
        fprintf('7:e ordningen (SBP-UW): ||err(m=%d)||_h = %.8e\n', m, Err_7th);
    end
    
    % (Valfritt: Spara felet till en extern fil/logg om du vill)

end % Slut på order_idx loop

% --- HUR MAN BERÄKNAR q MANUELLT ---
fprintf('\n-----------------------------------------------------\n');
fprintf('För att beräkna konvergenshastigheten (q) måste du:\n');
fprintf('1. Spara felet (Err_h) från denna körning (m_k).\n');
fprintf('2. Ändra m till nästa värde (m_{k+1}, t.ex. 201) och kör igen.\n');
fprintf('3. Använd formeln: q = log10(Err_k / Err_{k+1}) / log10(m_{k+1} / m_k)\n');
