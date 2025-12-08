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
Ll = kron(e1 + e2, e_1');   %e_1' är conjugate transpose
Lr = kron(e1 - e2, e_m');
L = [Ll;
    Lr];
% V_t = -PC^-1(D_x + D)PV = MV
% M = -PC^-1(D_x + D)P

% Nedan högst oklart
domain_width = 1;   %Detta är x_r - x_l
h = domain_width/(m - 1);
%H = eye(m);
%H(1, 1) = 1/2;
%H(m, m) = 1/2;
%H = h*H;

[H, HI, Dp, Dm, ~, ~] = SBP7_Upwind(m, h);

H_bar = kron(eye(2), H);  %Kan strunta i C eftersom den blir ehnhetsmatris (2m x 2m) eftersom roh = c =1

%Definera P
P = eye(2*m) - inv(H_bar)*transpose(L)*inv((L*inv(H_bar)*transpose(L)))*L;

%%

% Sätt ihop till M
D_x = [zeros(m), Dp;
    Dm, zeros(m)];
M = -P*(D_x)*P;     % Struntar i C invers eftersom det blir enhetsmatrisen, D-matrisen försvinner eftersom beta=0  

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
%Definerar här om L till Dirichlet
Ll = kron(e2, e_1');
Lr = kron(e2, e_m');
L = [Ll;
    Lr];
%Definera P
P = eye(2*m) - inv(H_bar)*transpose(L)*inv((L*inv(H_bar)*transpose(L)))*L;

% Sätt ihop till M
D_x = [zeros(m), Dp;
    Dm, zeros(m)];
M = -P*(D_x)*P;     % Struntar i C invers eftersom det blir enhetsmatrisen, D-matrisen försvinner eftersom beta=0  

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

