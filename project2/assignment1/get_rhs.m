% Generates RHS-function and post processing projection filtering
function [rhs,apply_P_full_div] = get_rhs(ops,Re,gin_data,gin_data_t,nnz_tolerance)

ops_loc = ops.ops_loc;
mtot = ops.mtot;

% Boundary divergence operators
Ldiv1 = [
    ops_loc{1}.dw',ops_loc{1}.SBPy.D1*ops_loc{1}.ew';
    ops_loc{1}.de',ops_loc{1}.SBPy.D1*ops_loc{1}.ee';
    ops_loc{1}.SBPx.D1*ops_loc{1}.es',ops_loc{1}.ds';
    ops_loc{1}.SBPx.D1*ops_loc{1}.en',ops_loc{1}.dn'];

Ldiv2 = [
    ops_loc{2}.dw',ops_loc{2}.SBPy.D1*ops_loc{2}.ew';
    ops_loc{2}.de',ops_loc{2}.SBPy.D1*ops_loc{2}.ee';
    ops_loc{2}.SBPx.D1*ops_loc{2}.es',ops_loc{2}.ds';
    ops_loc{2}.SBPx.D1*ops_loc{2}.en',ops_loc{2}.dn'];

Ldiv3 = [
    ops_loc{3}.dw',ops_loc{3}.SBPy.D1*ops_loc{3}.ew';
    ops_loc{3}.de',ops_loc{3}.SBPy.D1*ops_loc{3}.ee';
    ops_loc{3}.SBPx.D1*ops_loc{3}.es',ops_loc{3}.ds';
    ops_loc{3}.SBPx.D1*ops_loc{3}.en',ops_loc{3}.dn'];

% Whole domain divergence operators
Ldiv_full1 = [ops_loc{1}.Dx,ops_loc{1}.Dy];
Ldiv_full2 = [ops_loc{2}.Dx,ops_loc{2}.Dy];
Ldiv_full3 = [ops_loc{3}.Dx,ops_loc{3}.Dy];

% Boundary condition operators, Dirichlet everywhere except outflow
Lbc1 = [ops_loc{1}.ew';ops_loc{1}.es';ops_loc{1}.en'];
Lbc1 = [Lbc1,sparse(size(Lbc1,1),ops_loc{1}.mtot);
    sparse(size(Lbc1,1),ops_loc{1}.mtot),Lbc1];

Lbc2 = [ops_loc{2}.en'];
Lbc2 = [Lbc2,sparse(size(Lbc2,1),ops_loc{2}.mtot);
    sparse(size(Lbc2,1),ops_loc{2}.mtot),Lbc2];

Lbc3 = [ops_loc{3}.ew';ops_loc{3}.es'];
Lbc3 = [Lbc3,sparse(size(Lbc3,1),ops_loc{3}.mtot);
    sparse(size(Lbc3,1),ops_loc{3}.mtot),Lbc3];

% Interface conditions operators, pure projection
Lint12 = [
    ops_loc{1}.ee',sparse(ops_loc{1}.my,ops_loc{1}.mtot),-ops_loc{2}.ew',sparse(ops_loc{2}.my,ops_loc{2}.mtot+2*ops_loc{3}.mtot);
    sparse(ops_loc{1}.my,ops_loc{1}.mtot),ops_loc{1}.ee',sparse(ops_loc{2}.my,ops_loc{2}.mtot),-ops_loc{2}.ew',sparse(ops_loc{2}.my,2*ops_loc{3}.mtot);
    ops_loc{1}.de',sparse(ops_loc{1}.my,ops_loc{1}.mtot),-ops_loc{2}.dw',sparse(ops_loc{2}.my,ops_loc{2}.mtot+2*ops_loc{3}.mtot);
    sparse(ops_loc{1}.my,ops_loc{1}.mtot),ops_loc{1}.de',sparse(ops_loc{2}.my,ops_loc{2}.mtot),-ops_loc{2}.dw',sparse(ops_loc{2}.my,2*ops_loc{3}.mtot);];

Lint23 = [
    sparse(ops_loc{2}.mx,2*ops_loc{1}.mtot),-ops_loc{2}.es',sparse(ops_loc{2}.mx,ops_loc{2}.mtot),ops_loc{3}.en',sparse(ops_loc{2}.mx,ops_loc{3}.mtot);
    sparse(ops_loc{2}.mx,2*ops_loc{1}.mtot),sparse(ops_loc{2}.mx,ops_loc{2}.mtot),-ops_loc{2}.es',sparse(ops_loc{2}.mx,ops_loc{3}.mtot),ops_loc{3}.en';
    sparse(ops_loc{2}.mx,2*ops_loc{1}.mtot),-ops_loc{2}.ds',sparse(ops_loc{2}.mx,ops_loc{2}.mtot),ops_loc{3}.dn',sparse(ops_loc{2}.mx,ops_loc{3}.mtot);
    sparse(ops_loc{2}.mx,2*ops_loc{1}.mtot),sparse(ops_loc{2}.mx,ops_loc{2}.mtot),-ops_loc{2}.ds',sparse(ops_loc{2}.mx,ops_loc{3}.mtot),ops_loc{3}.dn'];

% Assemble full operator with boundary divergence
L = [Lbc1,sparse(size(Lbc1,1),2*(ops_loc{2}.mtot + ops_loc{3}.mtot));
    sparse(size(Lbc2,1),2*ops_loc{1}.mtot),Lbc2,sparse(size(Lbc2,1),2*ops_loc{3}.mtot);
    sparse(size(Lbc3,1),2*(ops_loc{1}.mtot + ops_loc{2}.mtot)),Lbc3;
    Ldiv1,sparse(size(Ldiv1,1),2*(ops_loc{2}.mtot + ops_loc{3}.mtot));
    sparse(size(Ldiv2,1),2*ops_loc{1}.mtot),Ldiv2,sparse(size(Ldiv2,1),2*ops_loc{3}.mtot);
    sparse(size(Ldiv3,1),2*(ops_loc{1}.mtot + ops_loc{2}.mtot)),Ldiv3;
    Lint12;
    Lint23];

% Assemble full operator with global divergence
L_div_full = [Lbc1,sparse(size(Lbc1,1),2*(ops_loc{2}.mtot + ops_loc{3}.mtot));
    sparse(size(Lbc2,1),2*ops_loc{1}.mtot),Lbc2,sparse(size(Lbc2,1),2*ops_loc{3}.mtot);
    sparse(size(Lbc3,1),2*(ops_loc{1}.mtot + ops_loc{2}.mtot)),Lbc3;
    Ldiv_full1,sparse(size(Ldiv_full1,1),2*(ops_loc{2}.mtot + ops_loc{3}.mtot));
    sparse(size(Ldiv_full2,1),2*ops_loc{1}.mtot),Ldiv_full2,sparse(size(Ldiv_full2,1),2*ops_loc{3}.mtot);
    sparse(size(Ldiv_full3,1),2*(ops_loc{1}.mtot + ops_loc{2}.mtot)),Ldiv_full3;
    Lint12;
    Lint23];

% Remove linearly dependent rows
L = remove_deps(L,ops.H_2comp,1e-6);

% Compute projection operators
HI_2comp = inv(ops.H_2comp);
Lplus = HI_2comp*L'*inv(L*HI_2comp*L');
Lplus = set2zero(Lplus,nnz_tolerance,0);
P = speye(2*mtot) - Lplus*L;
P = set2zero(P,nnz_tolerance,0);

tmp = L_div_full*HI_2comp*L_div_full';
Lplus_div_full = HI_2comp*L_div_full'*sparse(pinv(full(tmp),1e-2));
Lplus_div_full = set2zero(Lplus_div_full,nnz_tolerance,0);
P_full_div = speye(2*mtot) - Lplus_div_full*L_div_full;
P_full_div = set2zero(P_full_div,nnz_tolerance,0);

Lplus_small = Lplus(:,1:ops.ops_loc{1}.my);
Lplus_div_full_small = Lplus_div_full(:,1:ops.ops_loc{1}.my);

apply_P = @(w,gin) P*w + Lplus_small*gin;
apply_P_full_div = @(w,gin) P_full_div*w + Lplus_div_full_small*gin;

    % Define RHS function of ODE with pressure already computed
    function [wt,p] = rhs_(w,t,p)
        % Inner projection
        w = apply_P(w,gin_data(t));
        
        % Split solution vector
        u1 = w(1:ops_loc{1}.mtot);
        v1 = w(ops_loc{1}.mtot+1:2*ops_loc{1}.mtot);
        u2 = w(2*ops_loc{1}.mtot+1:2*ops_loc{1}.mtot + ops_loc{2}.mtot);
        v2 = w(2*ops_loc{1}.mtot + ops_loc{2}.mtot+1:2*ops_loc{1}.mtot + 2*ops_loc{2}.mtot);
        u3 = w(2*ops_loc{1}.mtot + 2*ops_loc{2}.mtot+1:2*ops_loc{1}.mtot + 2*ops_loc{2}.mtot + ops_loc{3}.mtot);
        v3 = w(2*ops_loc{1}.mtot + 2*ops_loc{2}.mtot + ops_loc{3}.mtot+1:end);

        % Outflow SATs (homogeneous characteristic bc)
        ue2 = ops_loc{2}.ee'*u2;
        ve2 = ops_loc{2}.ee'*v2;
        Lambda_5_e_2 = (ue2 - sqrt(ue2.^2 + 8))/2;
        Lambda_2_e_2 = (ue2 - sqrt(ue2.^2 + 4))/2;

        ue3 = ops_loc{3}.ee'*u3;
        ve3 = ops_loc{3}.ee'*v3;
        Lambda_5_e_3 = (ue3 - sqrt(ue3.^2 + 8))/2;
        Lambda_2_e_3 = (ue3 - sqrt(ue3.^2 + 4))/2;

        p1 = p(1:ops_loc{1}.mtot);
        p2 = p(ops_loc{1}.mtot+1:ops_loc{1}.mtot+ops_loc{2}.mtot);
        p3 = p(ops_loc{1}.mtot+ops_loc{2}.mtot+1:end);

        SATu2 = -ops_loc{2}.HI*ops_loc{2}.ee*ops_loc{2}.He*...
            (-Lambda_5_e_2.*ue2 + 1/Re*ops_loc{2}.de'*u2 - ops_loc{2}.ee'*p2);
        SATv2 = -ops_loc{2}.HI*ops_loc{2}.ee*ops_loc{2}.He*...
            (-Lambda_2_e_2.*ve2 + 1/Re*ops_loc{2}.de'*v2);
        SATu3 = -ops_loc{3}.HI*ops_loc{3}.ee*ops_loc{3}.He*...
            (-Lambda_5_e_3.*ue3 + 1/Re*ops_loc{3}.de'*u3 - ops_loc{3}.ee'*p3);
        SATv3 = -ops_loc{3}.HI*ops_loc{3}.ee*ops_loc{3}.He*...
            (-Lambda_2_e_3.*ve3 + 1/Re*ops_loc{3}.de'*v3);

        % Evaluate momentum equations
        u1_sq = u1.*u1;
        v1_sq = v1.*v1;
        u1v1 = u1.*v1;
        ux1 = ops_loc{1}.Dx*u1;
        uy1 = ops_loc{1}.Dy*u1;
        vx1 = ops_loc{1}.Dx*v1;
        vy1 = ops_loc{1}.Dy*v1;

        u2_sq = u2.*u2;
        v2_sq = v2.*v2;
        u2v2 = u2.*v2;
        ux2 = ops_loc{2}.Dx*u2;
        uy2 = ops_loc{2}.Dy*u2;
        vx2 = ops_loc{2}.Dx*v2;
        vy2 = ops_loc{2}.Dy*v2;

        u3_sq = u3.*u3;
        v3_sq = v3.*v3;
        u3v3 = u3.*v3;
        ux3 = ops_loc{3}.Dx*u3;
        uy3 = ops_loc{3}.Dy*u3;
        vx3 = ops_loc{3}.Dx*v3;
        vy3 = ops_loc{3}.Dy*v3;

        u1t = -0.5*(ops_loc{1}.Dx*u1_sq + v1.*uy1 + ops_loc{1}.Dy*u1v1 - vy1.*u1) ...
            - ops_loc{1}.Dx*p1 + 1/Re*ops_loc{1}.DL*u1;
        v1t = -0.5*(u1.*vx1 + ops_loc{1}.Dx*u1v1 - ux1.*v1 + ops_loc{1}.Dy*v1_sq) ...
            - ops_loc{1}.Dy*p1 + 1/Re*ops_loc{1}.DL*v1;
        u2t = -0.5*(ops_loc{2}.Dx*u2_sq + v2.*uy2 + ops_loc{2}.Dy*u2v2 - vy2.*u2) ...
            - ops_loc{2}.Dx*p2 + 1/Re*ops_loc{2}.DL*u2 + SATu2;
        v2t = -0.5*(u2.*vx2 + ops_loc{2}.Dx*u2v2 - ux2.*v2 + ops_loc{2}.Dy*v2_sq) ...
            - ops_loc{2}.Dy*p2 + 1/Re*ops_loc{2}.DL*v2 + SATv2;
        u3t = -0.5*(ops_loc{3}.Dx*u3_sq + v3.*uy3 + ops_loc{3}.Dy*u3v3 - vy3.*u3) ...
            - ops_loc{3}.Dx*p3 + 1/Re*ops_loc{3}.DL*u3 + SATu3;
        v3t = -0.5*(u3.*vx3 + ops_loc{3}.Dx*u3v3 - ux3.*v3 + ops_loc{3}.Dy*v3_sq) ...
            - ops_loc{3}.Dy*p3 + 1/Re*ops_loc{3}.DL*v3 + SATv3;

        wt = [u1t;v1t;u2t;v2t;u3t;v3t];

        % Outer projection
        wt = apply_P(wt,gin_data_t(t));
    end

rhs = @(w,t,pold) rhs_(w,t,pold);
end