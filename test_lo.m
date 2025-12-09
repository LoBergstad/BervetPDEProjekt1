m = 51;     %Detta ska även testas för 51, 101

%   Givna kosntanter
roh = 1;
c = 1;
beta = 0;

%   Definera matriser
A = [0, 1;
    1, 0];
C = [1/(roh*c^2), 0;
    0, roh];
D = [beta, 0;
    0, 0];
e_1 = zeros(m,1);
e_1(1) = 1;
e_m = zeros(m, 1);
e_m(m) = 1;
e1 = [1, 0];    % Detta är e^(1)
e2 = [0, 1];    % e^(2)
% Definera L för characteristic
Ll_c = kron(e1 + e2, e_1');   %e_1' är conjugate transpose
Lr_c = kron(e1 - e2, e_m');
L_c = [Ll_c;
    Lr_c];
% V_t = -PC^-1(D_x + D)PV = MV
% M = -PC^-1(D_x + D)P

domain_width = 1;   %Detta är x_r - x_l
h = domain_width/(m - 1);

[H, HI, Dp, Dm, ~, ~] = SBP7_Upwind(m, h);

H_bar = kron(eye(2), H);  %Kan strunta i C eftersom den blir ehnhetsmatris (2m x 2m) eftersom roh = c =1

%Definera P för characteristic
P_c = eye(2*m) - inv(H_bar)*transpose(L_c)*inv((L_c*inv(H_bar)*transpose(L_c)))*L_c;

%%

% Sätt ihop till M
D_x = [zeros(m), Dp;
    Dm, zeros(m)];

M = -P_c*(D_x)*P_c;     % Struntar i C invers eftersom det blir enhetsmatrisen, D-matrisen försvinner eftersom beta=0  

% Plotta egenvärden för Characteristic
M_eig_char = eig(M);
h_Meig_char = h*M_eig_char;
figure;
scatter(real(h_Meig_char), imag(h_Meig_char), 80, 'filled')
grid on
xlabel('Re')
ylabel('Im')
title('Eigenvalues of M using characteristic boundary conditions')
xlim([-0.5 0.2])
ylim([-2.5 2.5])

%%
% Dirichlet
Ll_d = kron(e2, e_1');
Lr_d = kron(e2, e_m');
L_d = [Ll_d;
    Lr_d];
%Definera P för dirichlet
P_d = eye(2*m) - inv(H_bar)*transpose(L_d)*inv((L_d*inv(H_bar)*transpose(L_d)))*L_d;

% Sätt ihop till M
D_x = [zeros(m), Dp;
    Dm, zeros(m)];
M = -P_d*(D_x)*P_d;     % Struntar i C invers eftersom det blir enhetsmatrisen, D-matrisen försvinner eftersom beta=0  

M_eig_dir = eig(M);
h_Meig_dir = h*M_eig_dir;
figure;
scatter(real(h_Meig_dir), imag(h_Meig_dir), 80, 'filled')
grid on
xlabel('Re')
ylabel('Im')
title('Eigenvalues of M using Dirichlet boundary conditions')
xlim([-0.5 0.2])
ylim([-2.5 2.5])

%% Plotta RK4 stability domain och egenvärden för de båda randvillkoren
% Characteristic

% --- Inställningar ---
N = 400;              % Upplösning
x = linspace(-4, 4, N);
y = linspace(-4, 4, N);
[X, Y] = meshgrid(x, y);
Z = X + 1i*Y;

% --- RK4 stabilitetsfunktion ---
R = 1 + Z + (Z.^2)/2 + (Z.^3)/6 + (Z.^4)/24;
S = abs(R);
% --- Plot ---
figure; hold on;
contour(X, Y, S, [1 1], 'b', 'LineWidth', 1.5); % Stabilitetsgränsen
scatter(real(h_Meig_char), imag(h_Meig_char), 80, 'filled', 'r') % Characteristic egenvärden
xlabel('Re(z)');
ylabel('Im(z)');
title('Stability domain for RK4 with eigenvalues of M (Characteristic BC)');
axis equal
grid on

%%
% --- Plot ---
figure; hold on;
contour(X, Y, S, [1 1], 'b', 'LineWidth', 1.5); % Stabilitetsgränsen
scatter(real(h_Meig_dir), imag(h_Meig_dir), 80, 'filled', 'r') % Dirichlet egenvärden
xlabel('Re(z)');
ylabel('Im(z)');
title('Stability domain for RK4 with eigenvalues of M (Dirichlet BC)');
axis equal
grid on

%% Hitta CFL m.h.a 2.78
%        % För characteristic, utan att 
%        lambda_max_char = max(abs(M_eig_char));
%        % Ska vara mindre än
%        R = 2.78;
%        k_char = R/lambda_max_char;  % = delta t
%        
%        
%        CFL_char = k_char/h
%        
%        % För dirichlet
%        lambda_max_dir = max(abs(M_eig_dir));
%        k_dir = R/lambda_max_dir;
%        CFL_dir = k_dir/h

%% Hitta CFL enligt RK4 stabilitetsekvationen
% För characteristic

lambda_char = M_eig_char;
k_char = 0.01;
z1 = k_char*lambda_char;
R_char = 1 + z1 + (z1.^2)/2 + (z1.^3)/6 + (z1.^4)/24;
while all(abs(R_char) <= 1)
    k_char = k_char + 0.000001;
    z1 = k_char*lambda_char;
    R_char = 1 + z1 + (z1.^2)/2 + (z1.^3)/6 + (z1.^4)/24;
end
k_char - 0.000001;
CFL_char = k_char/h

% För dirichlet
lambda_dir = M_eig_dir;
k_dir = 0.01;
z2 = k_dir * lambda_dir;
R_dir = 1 + z2 + (z2.^2)/2 + (z2.^3)/6 + (z2.^4)/24;
while all(abs(R_dir) <= 1)
    k_dir = k_dir + 0.000001;
    z2 = k_dir * lambda_dir;
    R_dir = 1 + z2 + (z2.^2)/2 + (z2.^3)/6 + (z2.^4)/24;
end
k_dir = k_dir - 0.000001;
CFL_dir = k_dir/h

