m = 2;     %Detta ska även testas för 51, 101

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
H = eye(m);
H(1, 1) = 1/2;
H(m, m) = 1/2;
H = h*H;

H_bar = kron(eye(2), H)*C;  %Detta blir fel. Vet ej hur man ska få till [H, 0; 0, H]


