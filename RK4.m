function [t, W] = rk4_linear(M, W0, tspan, k)
%RK4_LINEAR  Löser systemet W_t = M*W med Runge–Kutta 4
%
%   [t, W] = rk4_linear(M, W0, tspan, k)
%
%   M      - systemmatris (t.ex. från SBP-discretisering)
%   W0     - initialtillstånd, vektor [p; v] av längd 2m
%   tspan  - [t0, t_slut]
%   k      - tidssteg
%
%   Returnerar:
%   t  - tidsvektor (kolumnvektor)
%   W  - lösning, varje kolumn motsvarar W vid ett tidssteg:
%         W(:,i) = [p_i; v_i]

% --- Initiering ---
t0 = tspan(1);
t_end = tspan(2);
t = (t0:k:t_end).';
n_steps = length(t);

n_vars = length(W0);
W = zeros(n_vars, n_steps);   % varje kolumn = tillståndsvektor
W(:,1) = W0(:);

% --- RK4-loop ---
for i = 1:n_steps-1
    Wi = W(:,i);

    k1 = M * Wi;
    k2 = M * (Wi + 0.5*k*k1);
    k3 = M * (Wi + 0.5*k*k2);
    k4 = M * (Wi + k*k3);

    W(:,i+1) = Wi + (k/6)*(k1 + 2*k2 + 2*k3 + k4);
end
end
