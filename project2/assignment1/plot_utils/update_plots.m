function srfs = update_plots(w,p,ops,srfs)
u1 = w(1:ops.ops_loc{1}.mtot);
v1 = w(ops.ops_loc{1}.mtot+1:2*ops.ops_loc{1}.mtot);
u2 = w(2*ops.ops_loc{1}.mtot+1:2*ops.ops_loc{1}.mtot + ops.ops_loc{2}.mtot);
v2 = w(2*ops.ops_loc{1}.mtot + ops.ops_loc{2}.mtot+1:2*ops.ops_loc{1}.mtot + 2*ops.ops_loc{2}.mtot);
u3 = w(2*ops.ops_loc{1}.mtot + 2*ops.ops_loc{2}.mtot+1:2*ops.ops_loc{1}.mtot + 2*ops.ops_loc{2}.mtot + ops.ops_loc{3}.mtot);
v3 = w(2*ops.ops_loc{1}.mtot + 2*ops.ops_loc{2}.mtot + ops.ops_loc{3}.mtot+1:end);
p1 = p(1:ops.ops_loc{1}.mtot);
p2 = p(ops.ops_loc{1}.mtot+1:ops.ops_loc{1}.mtot+ops.ops_loc{2}.mtot);
p3 = p(ops.ops_loc{1}.mtot+ops.ops_loc{2}.mtot+1:end);

U1 = reshape(u1,[ops.ops_loc{1}.my,ops.ops_loc{1}.mx]);
V1 = reshape(v1,[ops.ops_loc{1}.my,ops.ops_loc{1}.mx]);

U2 = reshape(u2,[ops.ops_loc{2}.my,ops.ops_loc{2}.mx]);
V2 = reshape(v2,[ops.ops_loc{2}.my,ops.ops_loc{2}.mx]);

U3 = reshape(u3,[ops.ops_loc{3}.my,ops.ops_loc{3}.mx]);
V3 = reshape(v3,[ops.ops_loc{3}.my,ops.ops_loc{3}.mx]);

P1 = reshape(p1,[ops.ops_loc{1}.my,ops.ops_loc{1}.mx]);
P2 = reshape(p2,[ops.ops_loc{2}.my,ops.ops_loc{2}.mx]);
P3 = reshape(p3,[ops.ops_loc{3}.my,ops.ops_loc{3}.mx]);

srfs.srf1_u.CData = U1;
srfs.srf2_u.CData = U2;
srfs.srf3_u.CData = U3;
srfs.srf1_v.CData = V1;
srfs.srf2_v.CData = V2;
srfs.srf3_v.CData = V3;
srfs.srf1_p.CData = P1;
srfs.srf2_p.CData = P2;
srfs.srf3_p.CData = P3;
srfs.srf1_u.ZData = U1;
srfs.srf2_u.ZData = U2;
srfs.srf3_u.ZData = U3;
srfs.srf1_v.ZData = V1;
srfs.srf2_v.ZData = V2;
srfs.srf3_v.ZData = V3;
srfs.srf1_p.ZData = P1;
srfs.srf2_p.ZData = P2;
srfs.srf3_p.ZData = P3;

drawnow