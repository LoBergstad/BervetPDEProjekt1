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

for order = [6, 7]
    for m = [201, 401]
        for t_end = [1.5, 2.5]
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
            L = [kron(e1 + c1*e2, e_1'); 
                kron(e1 - c1*e2, e_m')];
            P = eye(2*m) - H_bar\L.' / (L/H_bar * L.')*L;
            M = -P/C*(D_x)*P;
            
            % Initial Values
            p0 = theta_1(x,0) - theta_2(x,0);
            v0 = theta_1(x,0) + theta_2(x,0);
            u0 = [p0 v0];

            [t, u] = RK4(M, u0, [0, t_end], alpha*h);

            p = u(1:m, end);
            v = u(m+1:2*m, end);
            
            figure;

            plot(x, p, 'r-', 'LineWidth', 1.5);
            hold on
            plot(x, v, 'b--', 'LineWidth', 1.5);
            hold off

            legend('p', 'v');
            %sprintf('%.0f th order SBP operator, m = %.0f, at t = %.1f', order, m, t_end)
            title(sprintf('%.0f th order SBP operator, m = %.0f, at t = %.1f', order, m, t_end));
            
            drawnow
        end
    end
end
