function srfs = setup_plot(w,p,ops,lims)
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

    figure
    subplot(3,1,1)
    hold on
    srf1_u = surf(ops.ops_loc{1}.X,ops.ops_loc{1}.Y,U1,U1,'Edgecolor','none');
    srf2_u = surf(ops.ops_loc{2}.X,ops.ops_loc{2}.Y,U2,U2,'Edgecolor','none');
    srf3_u = surf(ops.ops_loc{3}.X,ops.ops_loc{3}.Y,U3,U3,'Edgecolor','none');
    view(2)
    axis(lims)
    xlabel('x')
    ylabel('y')
    title('u (x-velocity)')
    colormap(linspecer)
    if isMATLABReleaseOlderThan('R2022a') % change to caxis if matlab version < R2022a 
        caxis([-0.1,1.5])
    else
        clim([-0.1,1.5])
    end
    colorbar

    subplot(3,1,2)
    hold on
    srf1_v = surf(ops.ops_loc{1}.X,ops.ops_loc{1}.Y,V1,V1,'Edgecolor','none');
    srf2_v = surf(ops.ops_loc{2}.X,ops.ops_loc{2}.Y,V2,V2,'Edgecolor','none');
    srf3_v = surf(ops.ops_loc{3}.X,ops.ops_loc{3}.Y,V3,V3,'Edgecolor','none');
    view(2)
    axis(lims)
    xlabel('x')
    ylabel('y')
    title('v (y-velocity)')
    colormap(linspecer)
    if isMATLABReleaseOlderThan('R2022a') % change to caxis if matlab version < R2022a 
        caxis([-0.25,0.1])
    else
        clim([-0.25,0.1])
    end
    colorbar

    subplot(3,1,3)
    hold on
    srf1_p = surf(ops.ops_loc{1}.X,ops.ops_loc{1}.Y,P1,P1,'Edgecolor','none');
    srf2_p = surf(ops.ops_loc{2}.X,ops.ops_loc{2}.Y,P2,P2,'Edgecolor','none');
    srf3_p = surf(ops.ops_loc{3}.X,ops.ops_loc{3}.Y,P3,P3,'Edgecolor','none');
    view(2)
    axis(lims)
    xlabel('x')
    ylabel('y')
    title('p (pressure)')
    colormap(linspecer)
    if isMATLABReleaseOlderThan('R2022a') % change to caxis if matlab version < R2022a 
        caxis([0.1,1])
    else
        clim([0.1,1])
    end
    colorbar

    srfs.srf1_u = srf1_u;
    srfs.srf1_v = srf1_v;
    srfs.srf1_p = srf1_p;
    srfs.srf2_u = srf2_u;
    srfs.srf2_v = srf2_v;
    srfs.srf2_p = srf2_p;
    srfs.srf3_u = srf3_u;
    srfs.srf3_v = srf3_v;
    srfs.srf3_p = srf3_p;