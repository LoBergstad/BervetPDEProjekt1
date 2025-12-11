%Task 1 Ass 2
rho_1 = 1;
c_1 = 1;
Z_1 = rho_1*c_1;

rho_2 = 2;
c_2 = 2;
Z_2 = rho_2*c_2;

T = 2*Z_2/(Z_1+Z_2);
R = (Z_2-Z_1)/(Z_1+Z_2);

r_s = 0.1;
x_l = -2;
x_r = 4;
x_0 = 1;

e1 = [1, 0];
e2 = [0, 1];

theta_1 = @(x, t) exp(-((x-t)/r_s).^2);
theta_2 = @(x, t) -exp(-((x+t)/r_s).^2);

order = 7;
m = 401;
t_end = 2.5;


alpha = 0.05; %Räkna ut kanske?

e_1 = [1; zeros(m-1,1)];
e_m = [zeros(m-1, 1); 1];

h = (x_r-x_l)/(m - 1);
x = linspace(x_l, x_r, m);

if order == 6
    [H, ~, D1, ~, ~, ~] = SBP6(m, h);
    D_x = [zeros(m), D1; D1, zeros(m)];  
else
    [H, ~, Dp, Dm, ~, ~] = SBP7_Upwind(m, h);
    D_x = [zeros(m), Dp; Dm, zeros(m)];  
end

m_half = floor(m/2);
C = diag([ ...
    (1/(rho_1*c_1^2)) * ones(m_half+1, 1);
    (1/(rho_2*c_2^2)) * ones(m_half, 1);
    rho_1 * ones(m_half+1, 1);
    rho_2 * ones(m_half, 1)
]);
%C = [1/(rho_1*c_1*c_1)*ones() 1/(rho_2*c_2*c_2)*ones(m), 0;
%    0, rho_1*ones(m) rho_2*ones(m)];


H_bar = kron(eye(2), H)*C;


% Characteristic BC
L = [kron(e1 + e2, e_1'); 
    kron(e1 - e2, e_m')];
P = eye(2*m) - H_bar\L.' / (L/H_bar * L.')*L;
M = -P/C*(D_x)*P;

% Initial Values
p0 = theta_1(x,0);% - theta_2(x,0);
v0 = theta_1(x,0);% + theta_2(x,0);
u0 = [p0 v0];

[t, u] = RK4(M, u0, [0, t_end], alpha*h);

%% ANIMATE SOLUTION


figure;

patch([1 4 4 1], [-1 -1 2 2], [0.8 0.8 0.8], ...
      'EdgeColor', 'none', 'HandleVisibility', 'off');     % light grey, no border
hold on;

p = u(1:m, 1);
v = u(m+1:2*m, 1);

p_plot = plot(x, p, 'r-', 'LineWidth', 1.5);
hold on
v_plot = plot(x, v, 'b--', 'LineWidth', 1.5);
hold off

ylim([-1 2]);

tit = title(sprintf('t = %.1f', t(1)));
legend('p', 'v');
drawnow;

frameTime = 0.005;
n = 5;
for k = 2:length(t)/n
    p_plot.YData = u(1:m, k*n);          
    v_plot.YData = u(m+1:2*m, k*n);
    tit.String = sprintf('t = %.2f', t(k*n));
    drawnow limitrate;
    pause(frameTime);
end

