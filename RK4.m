function [t, W] = rk4_linear(M, W0, tspan, k)
%RK4_LINEAR  Löser systemet W_t = M*W med Runge–Kutta 4
%
%   [t, W] = rk4_linear(M, W0, tspan, k)
%
%   M      - systemmatris (t.ex. från SBP-discretisering)
%   W0     - begynnelsevärde (vektor)
%   tspan  - [t0, t_slut]
%   k      - tidssteg
%
%   Returnerar:
%   t  - tidsvektor
%   W  - lösning, varje rad motsvarar ett tidssteg

t0 = tspan(1);
t_end = tspan(2);
t = (t0:k:t_end).';
n_steps = length(t);

n_vars = length(W0);
W = zeros(n_steps, n_vars);
W(1,:) = W0(:).';

for i = 1:n_steps-1
    Wi = W(i,:).';
    k1 = M * Wi;
    k2 = M * (Wi + 0.5*k*k1);
    k3 = M * (Wi + 0.5*k*k2);
    k4 = M * (Wi + k*k3);
    W(i+1,:) = (Wi + (k/6)*(k1 + 2*k2 + 2*k3 + k4)).';
end
end