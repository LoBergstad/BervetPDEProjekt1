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
        for t = [1.5, 2.5]
            disp(order)
            disp(m)
            disp(t)
            alpha = 0.05; %Räkna ut kanske?
            
            e_1 = [1; zeros(m-1,1)];
            e_m = [zeros(m-1, 1); 1];

            h = (x_r-x_l)/(m - 1);
            x = linspace(x_l, x_r, m);

            if order == 6
                %[H, HI, Dp, Dm, ~, ~] = SBP7_Upwind(m, h);
                [H, HI, Dp, Dm, ~, ~] = SBP6_CD(m, h);
            else
                [H, HI, Dp, Dm, ~, ~] = SBP7_Upwind(m, h);
            end

            H_bar = kron(eye(2), H);
            D_x = [zeros(m), Dp; Dm, zeros(m)];    
            
            % Characteristic BC
            L = [kron(e1 + e2, e_1'); 
                kron(e1 - e2, e_m')];
            P = eye(2*m) - H_bar\L.' / (L/H_bar * L.')*L;
            M = -P*(D_x)*P;
            
            % Initial Values
            p0 = theta_1(x,0) - theta_2(x,0);
            v0 = theta_1(x,0) + theta_2(x,0);
            u0 = [p0 v0];

            [t, u] = RK4(M, u0, [0, t], alpha*h);
            p = u(1:m, end);
            v = u(m+1:2*m, end);
            
            figure;

            plot(x, p, 'r-', 'LineWidth', 1.5);
            hold on
            plot(x, v, 'b--', 'LineWidth', 1.5);
            hold off

            legend('p', 'v');
            title(string(sprintf('%d th order SBP operator, m = %d, at t = %.1f', order, m, t)));
            
            drawnow

        end
    end
end
