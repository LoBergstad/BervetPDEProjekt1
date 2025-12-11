%Task 2 Ass 2

% Startvärden

order = 7;
m = 1001;
t_end = 5;

rho_1 = 1;
c_1 = 1;
beta_1 = 0;

rho_2 = 2;
c_2 = 2;
beta_2 = 10;

r_s = 0.1;
x_l = -4;
x_r = 6;

e1 = [1, 0];
e2 = [0, 1];
e_1 = [1; zeros(m-1,1)];
e_m = [zeros(m-1, 1); 1];

theta_1 = @(x, t) exp(-((x-t)/r_s).^2);
theta_2 = @(x, t) -exp(-((x+t)/r_s).^2);

alpha = 0.05; %Räkna ut kanske?

h = (x_r-x_l)/(m - 1);
x = linspace(x_l, x_r, m);

%% Simulera vågen

% Hitta C utifrån väggposition
x_wall = 1;
gap = 0.6;
thickness = 0.2;

x_a1 = x_wall;
x_a2 = x_wall+thickness;
x_b1 = x_wall+thickness+gap;
x_b2 = x_wall+2*thickness+gap;

wall = ((x >= x_a1 & x <= x_a2) | (x >= x_b1 & x <= x_b2));
is_wall = [wall, wall]';

C_1 = [(1/(rho_1*c_1^2))*ones(m, 1); rho_1*ones(m, 1)];
C_2 = [(1/(rho_2*c_2^2))*ones(m, 1); rho_2*ones(m, 1)];

C = diag(is_wall .* C_2 + (~is_wall) .* C_1);

D_1 = [beta_1*ones(m, 1); zeros(m, 1)];
D_2 = [beta_2*ones(m, 1); zeros(m, 1)];

D = diag(is_wall .* D_2 + (~is_wall) .* D_1);


% Derivator % Integrator
%Kolla för ordning 7
%[H, ~, Dp, Dm, ~, ~] = SBP7_Upwind(m, h);
%Kolla för ordning 6
[H, ~, D1, ~, ~, ~] = SBP6(m, h);
Dp = D1;
Dm = D1;

D_x = [zeros(m), Dp; Dm, zeros(m)]; 
H_bar = kron(eye(2), H)*C;


% Characteristic BC
L = [kron(e1 + e2, e_1');   kron(e1 - e2, e_m')];
P = eye(2*m) - H_bar\L.' / (L/H_bar * L.')*L;
M = -P/C*(D_x+D)*P;

% Initial Values
p0 = theta_1(x,0);% - theta_2(x,0);
v0 = theta_1(x,0);% + theta_2(x,0);
u0 = [p0 v0];

[t, u] = RK4(M, u0, [0, t_end], alpha*h);


%% Plotta
ymin = -1;
ymax = 1.5;

figure;
ylim([ymin ymax]);

% Rita väggen
patch([x_a1 x_a2 x_a2 x_a1], [ymin ymin ymax ymax], [0.8 0.8 0.8], ...
      'EdgeColor', 'none', 'HandleVisibility', 'off');
hold on;

patch([x_b1 x_b2 x_b2 x_b1], [ymin ymin ymax ymax], [0.8 0.8 0.8], ...
      'EdgeColor', 'none', 'HandleVisibility', 'off');

% Rita Initialvärden

p = plot(x, u(1:m, 1), 'r-', 'LineWidth', 1.5);
hold on
v = plot(x, u(m+1:2*m, 1), 'b--', 'LineWidth', 1.5);
hold off

tit = title(sprintf('t = %.1f', t(1)));
legend('p', 'v');

drawnow;

% Animera

frameTime = 0.005;
n = 5;
for k = 1:length(t)/n
    p.YData = u(1:m, k*n);          
    v.YData = u(m+1:2*m, k*n);
    tit.String = sprintf('t = %.2f', t(k*n));
    drawnow limitrate;
    pause(frameTime);
end

