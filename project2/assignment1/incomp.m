% Simulates incompressible Navier-Stokes equations in 2D on a backward
% facing step problem using second order central SBP finite differences.
% The computational domain consists of three blocks arranged as follows:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %                    %
%         %                    %
%    1    %         2          %
%         %                    %
%         %                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %                    %
%         %                    %
%         %         3          %
%         %                    %
%         %                    %
%         %%%%%%%%%%%%%%%%%%%%%%
% 
% Input parameters:
% m       - number of grid points on the inflow boundary (integer > 2)
% Re      - Reynolds number (positive float)
% Tend    - End time (positive float)
% Lout    - Length of outflow pipe (positive float)
% animate - if animating solution or not (boolean)

function incomp(m,Re,Tend,Lout,animate)
tic
% Add the paths to various functions
addpath('utils')
addpath('plot_utils')

% Numbers in operators below these values are set to zero, for efficiency
nnz_tolerance_rhs = 1e-10;
nnz_tolerance_ppe = 1e-8;

% Spatial domain
Lin = 1;
H1 = 0.5;
H2 = 0.5;

% Spatial limits in each block
xl = {-Lin,0,0};
xr = {0,Lout,Lout};
yl = {0,0,-H2};
yr = {H1,H1,0};

% Grid step size
h = H1/(m-1);

% Build spatial operators
fprintf("Computing spatial operators...")
ops = get_ops(h,xl,xr,yl,yr);
fprintf(" done!\n")

% Temporal domain
k = 0.02;
dt_try = k*Re*ops.ops_loc{1}.SBPx.h.^2;
[dt,mt] = alignedTimestep(dt_try,Tend);
tvec = linspace(0,Tend,mt+1);

% Inflow boundary data
tmax = 3;
slope = 1;
Vmax = 1.5;
R = (ops.ops_loc{1}.SBPy.x(1) + ops.ops_loc{1}.SBPy.x(end))/2;
r = ops.ops_loc{1}.SBPy.x - R;
t_scaling = @(t) ((t < tmax).*exp(-slope*(t-tmax).*(t-tmax)) + (t >= tmax));
gin_data = @(t) Vmax*(1 - r.^2/R.^2)*t_scaling(t);
gin_data_t = @(t) Vmax*(1 - r.^2/R.^2)*((t < tmax).*-slope*exp(-slope*(t - tmax)^2)*(2*t - 2*tmax));

% Function solving the PPE
fprintf("Computing PPE solver...")
solve_ppe = get_ppe_solver(ops,Re,gin_data_t,nnz_tolerance_ppe);
fprintf(" done!\n")

% RHS function in ODE and projection post process function
fprintf("Computing RHS function...")
[rhs,apply_P_full_div] = get_rhs(ops,Re,gin_data,gin_data_t,nnz_tolerance_rhs);
fprintf(" done!\n")

% Initial data
w = zeros(2*ops.mtot,1);
p = solve_ppe(w,0);

fprintf("Time elapsed in setup: %.2f s\n",toc)

if animate
    srfs = setup_plot(w,p,ops,[xl{1},xr{end},yl{end},yr{1}]);
end

% Time stepping loop
fprintf("---- Running simulation ----\n")
tic
for tidx = 1:mt
    t = tvec(tidx);

    % RK4 iteration
    k1 = rhs(w,t,p);
    k2 = rhs(w + 0.5*dt*k1, t + 0.5*dt,p);
    k3 = rhs(w + 0.5*dt*k2, t + 0.5*dt,p);
    k4 = rhs(w + dt*k3, t + dt,p);
    w = w + dt/6*(k1 + 2*k2 + 2*k3 + k4);

    % Filter divergence every 20th time step
    if mod(tidx,20) == 0 || tidx == mt
        w = apply_P_full_div(w,gin_data(t+dt));
    end

    % Solve PPE
    p = solve_ppe(w,t+dt);

    if mod(tidx,round(mt/10)) == 0
        fprintf("%d%% done, t = %.2f, elapsed time: %.2f s\n",round(100*tidx/mt),t+dt,toc);
    end

    % Update plot every 0.05 time units
    if animate && (mod(tidx,round(0.05/dt)) == 0 || tidx == mt)
        w = apply_P_full_div(w,gin_data(t+dt));
        p = solve_ppe(w,t+dt);
        srfs = update_plots(w,p,ops,srfs);
    end
end

if animate
    plot_streamlines(w,ops,[xl{1},xr{end},yl{end},yr{1}]);
end

% ----------
% u3 = w(2*ops.ops_loc{1}.mtot + 2*ops.ops_loc{2}.mtot+1:2*ops.ops_loc{1}.mtot + 2*ops.ops_loc{2}.mtot + ops.ops_loc{3}.mtot);
% U3 = reshape(u3,[ops.ops_loc{3}.my,ops.ops_loc{3}.mx]);
% x3 = ops.ops_loc{3}.SBPx.x;
% y3 = ops.ops_loc{3}.SBPy.x;
% 
% Uncomment the above four lines and add your own code here to automatically 
% compute the reattachment length (non-mandatory). All you should need to 
% compute the reattachment length is the x-velocity in block 3, which is 
% given by the matrix U3, and the vector x3 containing the x-values in block 3.
% ----------

