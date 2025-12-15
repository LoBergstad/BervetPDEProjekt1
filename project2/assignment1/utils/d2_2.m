function SBP = d2_2(m,lim)
x_l = lim{1};
x_r = lim{2};
L = x_r-x_l;
h = L/(m-1);

BP = 1;
if(m<2*BP)
    error(['Operator requires at least ' num2str(2*BP) ' grid points']);
end

e_1=sparse(m,1);e_1(1)=1;
e_m=sparse(m,1);e_m(m)=1;

H=(speye(m,m));H(1,1)=0.5;H(m,m)=0.5;
H=h*H;
HI=inv(H);

diags   = -1:1;
stencil = [-1/2 0 1/2];
D1 = stripeMatrix(stencil, diags, m);

D1(1,1) = -1;
D1(1,2) = 1;

D1(m,m-1) = -1;
D1(m,m)   = 1;

D1 = D1/h;

Q = H*D1 + 1/2*(e_1*e_1') - 1/2*(e_m*e_m');

diags   = -1:1;
stencil = [1 -2 1];
D2 = stripeMatrix(stencil, diags, m);

D2(1,1) = 1;
D2(1,2) = -2;
D2(1,3) = 1;
D2(m,m-2) = 1;
D2(m,m-1) = -2;
D2(m,m)   = 1;
D2 = D2/h^2;

S_U = [-3/2, 2, -1/2]/h;
S_1 = sparse(1,m);
S_1(1:3) = S_U;
S_m = sparse(1,m);
S_m(m-2:m) = fliplr(-S_U);

M = -H*D2-e_1*S_1+e_m*S_m;
d1_l = S_1';
d1_r = S_m';

SBP.H = H;
SBP.HI = HI;
SBP.D1 = D1;
SBP.D2 = D2;
SBP.e_l = e_1;
SBP.e_r = e_m;
SBP.d1_l = d1_l;
SBP.d1_r = d1_r;
SBP.x = linspace(x_l,x_r,m)';
SBP.h = L/(m-1);
SBP.m = m;
end