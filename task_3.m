m = 51;     %Detta ska även testas för 51, 101

%   Givna kosntanter
rho = 1;
c = 1;
beta = 0;

%   Definera matriser
A = [0, 1;
    1, 0];
C = [1/(rho*c^2), 0;
    0, rho];
D = [beta, 0;
    0, 0];
e_1 = zeros(m,1);
e_1(1) = 1;
e_m = zeros(m, 1);
e_m(m) = 1;
e1 = [1, 0];    % Detta är e^(1)
e2 = [0, 1];    % e^(2)

domain_width = 1;   %Detta är x_r - x_l
h = domain_width/(m - 1);

[H, HI, Dp, Dm, ~, ~] = SBP7_Upwind(m, h);

H_bar = kron(eye(2), H);  %Kan strunta i C eftersom den blir ehnhetsmatris (2m x 2m) eftersom rho = c =1


%%

D_x = [zeros(m), Dp;
       Dm, zeros(m)];

% Definera L för characteristic
Ll_c = kron(e1 + e2, e_1');   %e_1' är conjugate transpose
Lr_c = kron(e1 - e2, e_m');
L = [Ll_c;
    Lr_c];
P = eye(2*m) - H_bar\L.' / (L/H_bar * L.')*L;
M_c = -P*(D_x)*P;


% Dirichlet
Ll_d = kron(e2, e_1');
Lr_d = kron(e2, e_m');
L = [Ll_d;
    Lr_d];
P = eye(2*m) - H_bar\L.' / (L/H_bar * L.')*L;
M_d = -P*(D_x)*P;

%%
% RK4 Stabilitetsplot

% --- Inställningar ---
N = 400;              % Upplösning
x = linspace(-4, 4, N);
y = linspace(-4, 4, N);
[X, Y] = meshgrid(x, y);
Z = X + 1i*Y;

% --- RK4 stabilitetsfunktion ---
R = 1 + Z + (Z.^2)/2 + (Z.^3)/6 + (Z.^4)/24;
S = abs(R);

% --- Plot Characteristic ---
figure; hold on;
contour(X, Y, S, [1 1], 'b', 'LineWidth', 1.5); % Stabilitetsgränsen
scatter(real(h*eig(M_c)), imag(h*eig(M_c)), 5, 'filled', 'r') % Characteristic egenvärden
xlabel('Re(z)');
ylabel('Im(z)');
title('Stability domain for RK4 with eigenvalues of M (Characteristic BC)');
axis equal
grid on

%%
% --- Plot Dirichlet ---
figure; hold on;
contour(X, Y, S, [1 1], 'b', 'LineWidth', 1.5); % Stabilitetsgränsen
scatter(real(h*eig(M_d)), imag(h*eig(M_d)), 5, 'filled', 'r') % Dirichlet egenvärden
xlabel('Re(z)');
ylabel('Im(z)');
title('Stability domain for RK4 with eigenvalues of M (Dirichlet BC)');
axis equal
grid on

%% Hitta CFL enligt RK4 stabilitetsekvationen
% För characteristic
k = 0.01; %delta t
eigs = eig(M_d);
z = k*eigs;

R = 1 + z + (z.^2)/2 + (z.^3)/6 + (z.^4)/24;
while all(abs(R) <= 1)
    k = k + 0.000001;
    z = k*eigs;
    R = 1 + z + (z.^2)/2 + (z.^3)/6 + (z.^4)/24;
end
k = k - 0.000001
CFL_char = k/h
