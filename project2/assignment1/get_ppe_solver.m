% Generates function that solves the pressure Poisson equation (PPE)
function solve_ppe = get_ppe_solver(ops,Re,gin_data_t,nnz_tolerance)

ops_loc = ops.ops_loc;

% Boundary conditions operators for the pressure, Neumann everywhere except on outflow (Dirichlet)
Lbc1 = [ops_loc{1}.dw';ops_loc{1}.ds';ops_loc{1}.dn'];
Lbc2 = [ops_loc{2}.ee';ops_loc{2}.dn'];
Lbc3 = [ops_loc{3}.dw';ops_loc{3}.ee';ops_loc{3}.ds'];

% Interface condition operators, pure projection
Lint12 = [
    ops_loc{1}.ee',-ops_loc{2}.ew';
    ops_loc{1}.de',-ops_loc{2}.dw'];
Lint23 = [
    ops_loc{2}.es',-ops_loc{3}.en';
    ops_loc{2}.ds',-ops_loc{3}.dn'];

% Assemble full operator
L = [Lbc1,sparse(size(Lbc1,1),ops_loc{2}.mtot + ops_loc{3}.mtot);
    sparse(size(Lbc2,1),ops_loc{1}.mtot),Lbc2,sparse(size(Lbc2,1),ops_loc{3}.mtot);
    sparse(size(Lbc3,1),ops_loc{1}.mtot + ops_loc{2}.mtot),Lbc3;
    Lint12,sparse(size(Lint12,1),ops_loc{3}.mtot);
    sparse(size(Lint23,1),ops_loc{1}.mtot),Lint23];

% Compute projection operator
Lplus = ops.HI*L'*sparse(pinv(full(L*ops.HI*L'),1e-2));
Lplus = set2zero(Lplus,nnz_tolerance,0);
P = speye(ops.mtot) - Lplus*L;
P = set2zero(P,nnz_tolerance,0);

% Build Poisson LHS matrix
sigma = -1/ops.h;
DL = ops.Dxx + ops.Dyy;
A = ops.H*P*DL*P + sigma*ops.H*Lplus*L;

% Compute Cholesky decomposition
[R_chol,~,P_chol] = chol(-A);
R_chol = set2zero(R_chol,nnz_tolerance,0);

% Build temporary constant matrices
zeros_int = zeros(size(Lint12,1)+size(Lint23,1),1);
tmp1 = ops.H*P;
tmp2 = -ops.H*P*DL*Lplus + sigma*ops.H*Lplus;
tmp1 = set2zero(tmp1,nnz_tolerance,0);
tmp2 = set2zero(tmp2,nnz_tolerance,0);

    % Define function for solving the PPE
    function p = solve_ppe_(w,t)

        % Split solution vector
        u1 = w(1:ops_loc{1}.mtot);
        v1 = w(ops_loc{1}.mtot+1:2*ops_loc{1}.mtot);
        u2 = w(2*ops_loc{1}.mtot+1:2*ops_loc{1}.mtot + ops_loc{2}.mtot);
        v2 = w(2*ops_loc{1}.mtot + ops_loc{2}.mtot+1:2*ops_loc{1}.mtot + 2*ops_loc{2}.mtot);
        u3 = w(2*ops_loc{1}.mtot + 2*ops_loc{2}.mtot+1:2*ops_loc{1}.mtot + 2*ops_loc{2}.mtot + ops_loc{3}.mtot);
        v3 = w(2*ops_loc{1}.mtot + 2*ops_loc{2}.mtot + ops_loc{3}.mtot+1:end);

        u = [u1;u2;u3];
        v = [v1;v2;v3];

        % Compute boundary conditions for the pressure
        ue2 = ops_loc{2}.ee'*u2;
        Lambda_5_e_2 = (ue2 - sqrt(ue2.^2 + 8))/2;

        ue3 = ops_loc{3}.ee'*u3;
        Lambda_5_e_3 = (ue3 - sqrt(ue3.^2 + 8))/2;

        u1_sq = u1.*u1;
        v1_sq = v1.*v1;
        u1v1 = u1.*v1;
        ux1 = ops_loc{1}.Dx*u1;
        uy1 = ops_loc{1}.Dy*u1;
        vx1 = ops_loc{1}.Dx*v1;
        vy1 = ops_loc{1}.Dy*v1;

        v2_sq = v2.*v2;
        u2v2 = u2.*v2;
        ux2 = ops_loc{2}.Dx*u2;
        vx2 = ops_loc{2}.Dx*v2;

        u3_sq = u3.*u3;
        v3_sq = v3.*v3;
        u3v3 = u3.*v3;
        ux3 = ops_loc{3}.Dx*u3;
        uy3 = ops_loc{3}.Dy*u3;
        vx3 = ops_loc{3}.Dx*v3;
        vy3 = ops_loc{3}.Dy*v3;

        gw1 = -gin_data_t(t) + ops_loc{1}.ew'*(...
            - 0.5*(ops_loc{1}.Dx*u1_sq) ...
            - 0.5*(v1.*uy1 + ops_loc{1}.Dy*u1v1 - vy1.*u1) ...
            + 1/Re*ops_loc{1}.DL*u1);
        gs1 = ops_loc{1}.es'*(...
            - 0.5*(u1.*vx1 + ops_loc{1}.Dx*u1v1 - ux1.*v1) ...
            - 0.5*(ops_loc{1}.Dy*v1_sq) ...
            + 1/Re*ops_loc{1}.DL*v1);
        gn1 = ops_loc{1}.en'*(...
            - 0.5*(u1.*vx1 + ops_loc{1}.Dx*u1v1 - ux1.*v1) ...
            - 0.5*(ops_loc{1}.Dy*v1_sq) ...
            + 1/Re*ops_loc{1}.DL*v1);
        
        ge2 = -Lambda_5_e_2.*ue2 + 1/Re*ops_loc{2}.de'*u2;
        gn2 = ops_loc{2}.en'*(...
            - 0.5*(u2.*vx2 + ops_loc{2}.Dx*u2v2 - ux2.*v2) ...
            - 0.5*(ops_loc{2}.Dy*v2_sq) ...
            + 1/Re*ops_loc{2}.DL*v2);

        gw3 = ops_loc{3}.ew'*(...
            - 0.5*(ops_loc{3}.Dx*u3_sq) ...
            - 0.5*(v3.*uy3 + ops_loc{3}.Dy*u3v3 - vy3.*u3) ...
            + 1/Re*ops_loc{3}.DL*u3);
        ge3 = -Lambda_5_e_3.*ue3 + 1/Re*ops_loc{3}.de'*u3;
        gs3 = ops_loc{3}.es'*(...
            - 0.5*(u3.*vx3 + ops_loc{3}.Dx*u3v3 - ux3.*v3) ...
            - 0.5*(ops_loc{3}.Dy*v3_sq) ...
            + 1/Re*ops_loc{3}.DL*v3);

        % Assemble Poisson equation
        gbc = [gw1;gs1;gn1;ge2;gn2;gw3;ge3;gs3];
        g = [gbc;zeros_int];
        f = -(ops.Dx*u).^2 - 2*(ops.Dy*u).*(ops.Dx*v) - (ops.Dy*v).^2;
        b = tmp1*f + tmp2*g;

        % Solve PPE
        p = -P_chol*(R_chol\(R_chol'\(P_chol'*b)));
    end
solve_ppe = @(w,t) solve_ppe_(w,t);
end