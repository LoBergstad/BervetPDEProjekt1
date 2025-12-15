% Constructs LHS matrix and RHS vector for Poisson system on a 2D square
% box.
function [A,b] = get_system(m)
% Add the paths to various functions
addpath('utils')

% Spatial limits
xl = 0;
xr = 1;
yl = 0;
yr = 1;

% 2D SBP operators
SBPx = d2_2(m,{xl,xr});
SBPy = d2_2(m,{yl,yr});
Imx = speye(m);
Imy = speye(m);
Dxx = kr(SBPx.D2,Imy);
Dyy = kr(Imx,SBPy.D2);
DL = Dxx + Dyy;
H = kr(SBPx.H,SBPy.H);
HI = inv(H);

dw = kr(SBPx.d1_l,Imy);
ee = kr(SBPx.e_r,Imy);
ds = kr(Imx,SBPy.d1_l);
dn = kr(Imx,SBPy.d1_r);

L = [dw';ee';ds';dn'];
L = remove_deps(L,H,1e-14);

Lplus = HI*L'*inv(L*HI*L');
P = speye(m*m) - Lplus*L;

% Build Poisson LHS matrix
sigma = -1/SBPx.h;
A = H*P*DL*P + sigma*H*Lplus*L;
A = -A;

% Build RHS vector
b = ones(m*m,1);
