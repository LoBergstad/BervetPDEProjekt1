%please snälla se detta på github!!
%Saker vi vet
rho = 1;
c = 1;
beta = 0;
r_star = 0.1;
x_l = -1;
x_r = 1;
t = 1.8;
CFL = 0.05;

L = x_r -x_l;

%antal punkter
m = 101; %[101, 201, 401, 801]
n = m-1; %Kolla om det är sant

u_t_star = 1;
u_m = 1;
u_n = 1;

%error för m punkter
err_m = u_t_star - u_m;
err_m_norm = sqrt(h) * norm(err_m);

%error för n punkter
err_n = u_t_star - u_n;
err_n_norm = sqrt(h) * norm(err_n);


%Convergence rate (q):
q = log10(err_m_norm/err_n_norm) / log10(n/m);


%Gaussian profiles:
theta_1 = @(x, t) exp(-((x-t)/r_star).^2);
theta_2 = @(x, t) -exp(-((x+t)/r_star).^2);



