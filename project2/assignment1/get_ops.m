% Construct two-dimensional SBP operators (local and global)
function ops = get_ops(h,xl,xr,yl,yr)

% Increase the outflow pipe length for outflow BC region
xr{2} = xr{2} + 0;
xr{3} = xr{3} + 0;
for bidx = 1:length(xl)
    mx{bidx} = ceil(1 + (xr{bidx} - xl{bidx})/h);
    my{bidx} = ceil(1 + (yr{bidx} - yl{bidx})/h);
    mtot{bidx} = mx{bidx}*my{bidx};
end

nBlocks = length(mx);

% Build local operators
for bidx = 1:nBlocks
    SBPx = d2_2(mx{bidx},{xl{bidx},xr{bidx}});
    SBPy = d2_2(my{bidx},{yl{bidx},yr{bidx}});
    SBPx.DI = -SBPx.H\(transpose(SBPx.D2*SBPx.h^2)*SBPx.D2*SBPx.h^2);
    SBPy.DI = -SBPy.H\(transpose(SBPy.D2*SBPy.h^2)*SBPy.D2*SBPy.h^2);

    [X,Y] = meshgrid(SBPx.x,SBPy.x);
    Imx = speye(mx{bidx});
    Imy = speye(my{bidx});

    Dx = kr(SBPx.D1,Imy);
    Dy = kr(Imx,SBPy.D1);
    DIx = kr(SBPx.DI,Imy);
    DIy = kr(Imx,SBPy.DI);
    Dxx = kr(SBPx.D2,Imy);
    Dyy = kr(Imx,SBPy.D2);
    DL = Dxx + Dyy;
    H = kr(SBPx.H,SBPy.H);

    ew = kr(SBPx.e_l,Imy);
    ee = kr(SBPx.e_r,Imy);
    es = kr(Imx,SBPy.e_l);
    en = kr(Imx,SBPy.e_r);
    dw = kr(SBPx.d1_l,Imy);
    de = kr(SBPx.d1_r,Imy);
    ds = kr(Imx,SBPy.d1_l);
    dn = kr(Imx,SBPy.d1_r);
    Hw = SBPy.H;
    He = SBPy.H;
    Hs = SBPx.H;
    Hn = SBPx.H;

    ops_loc{bidx}.SBPx = SBPx;
    ops_loc{bidx}.SBPy = SBPy;
    ops_loc{bidx}.Dx = Dx;
    ops_loc{bidx}.Dy = Dy;
    ops_loc{bidx}.DIx = DIx;
    ops_loc{bidx}.DIy = DIy;
    ops_loc{bidx}.Dxx = Dxx;
    ops_loc{bidx}.Dyy = Dyy;
    ops_loc{bidx}.DL = DL;
    ops_loc{bidx}.H = H;
    ops_loc{bidx}.HI = inv(H);
    ops_loc{bidx}.ew = ew;
    ops_loc{bidx}.ee = ee;
    ops_loc{bidx}.es = es;
    ops_loc{bidx}.en = en;
    ops_loc{bidx}.dw = dw;
    ops_loc{bidx}.de = de;
    ops_loc{bidx}.ds = ds;
    ops_loc{bidx}.dn = dn;
    ops_loc{bidx}.Hw = Hw;
    ops_loc{bidx}.He = He;
    ops_loc{bidx}.Hs = Hs;
    ops_loc{bidx}.Hn = Hn;
    ops_loc{bidx}.mx = mx{bidx};
    ops_loc{bidx}.my = my{bidx};
    ops_loc{bidx}.mtot = mx{bidx}*my{bidx};
    ops_loc{bidx}.X = X;
    ops_loc{bidx}.Y = Y;
end

% Build global operators
H = cell(nBlocks, nBlocks);
Dx = cell(nBlocks, nBlocks);
Dy = cell(nBlocks, nBlocks);
DIx = cell(nBlocks, nBlocks);
DIy = cell(nBlocks, nBlocks);
Dxx = cell(nBlocks, nBlocks);
Dyy = cell(nBlocks, nBlocks);
ops.mtot = 0;
for bidx = 1:nBlocks
    H{bidx,bidx} = ops_loc{bidx}.H;
    H_2comp{2*(bidx-1)+1,2*(bidx-1)+1} = ops_loc{bidx}.H;
    H_2comp{2*(bidx-1)+2,2*(bidx-1)+2} = ops_loc{bidx}.H;
    Dx{bidx,bidx} = ops_loc{bidx}.Dx;
    Dy{bidx,bidx} = ops_loc{bidx}.Dy;
    Dxx{bidx,bidx} = ops_loc{bidx}.Dxx;
    Dyy{bidx,bidx} = ops_loc{bidx}.Dyy;
    ops.mtot = ops.mtot + ops_loc{bidx}.mtot;
end
ops.H = toMatrix(H);
ops.H_2comp = toMatrix(H_2comp);
ops.HI = inv(ops.H);
ops.Dx = toMatrix(Dx);
ops.Dy = toMatrix(Dy);
ops.Dxx = toMatrix(Dxx);
ops.Dyy = toMatrix(Dyy);
ops.ops_loc = ops_loc;
ops.h = h;