function plot_streamlines(w,ops,lims)

u1 = w(1:ops.ops_loc{1}.mtot);
v1 = w(ops.ops_loc{1}.mtot+1:2*ops.ops_loc{1}.mtot);
u2 = w(2*ops.ops_loc{1}.mtot+1:2*ops.ops_loc{1}.mtot + ops.ops_loc{2}.mtot);
v2 = w(2*ops.ops_loc{1}.mtot + ops.ops_loc{2}.mtot+1:2*ops.ops_loc{1}.mtot + 2*ops.ops_loc{2}.mtot);
u3 = w(2*ops.ops_loc{1}.mtot + 2*ops.ops_loc{2}.mtot+1:2*ops.ops_loc{1}.mtot + 2*ops.ops_loc{2}.mtot + ops.ops_loc{3}.mtot);
v3 = w(2*ops.ops_loc{1}.mtot + 2*ops.ops_loc{2}.mtot + ops.ops_loc{3}.mtot+1:end);

U1 = reshape(u1,[ops.ops_loc{1}.my,ops.ops_loc{1}.mx]);
V1 = reshape(v1,[ops.ops_loc{1}.my,ops.ops_loc{1}.mx]);

U2 = reshape(u2,[ops.ops_loc{2}.my,ops.ops_loc{2}.mx]);
V2 = reshape(v2,[ops.ops_loc{2}.my,ops.ops_loc{2}.mx]);

U3 = reshape(u3,[ops.ops_loc{3}.my,ops.ops_loc{3}.mx]);
V3 = reshape(v3,[ops.ops_loc{3}.my,ops.ops_loc{3}.mx]);

X1 = ops.ops_loc{1}.X;
X2 = ops.ops_loc{2}.X;
X3 = ops.ops_loc{3}.X;
Y1 = ops.ops_loc{1}.Y;
Y2 = ops.ops_loc{2}.Y;
Y3 = ops.ops_loc{3}.Y;

U4 = nan*zeros(ops.ops_loc{3}.my,ops.ops_loc{1}.mx);
V4 = nan*zeros(ops.ops_loc{3}.my,ops.ops_loc{1}.mx);
[X4,Y4] = meshgrid(linspace(min(min(X1)),max(max(X1)),ops.ops_loc{1}.mx),linspace(min(min(Y3)),max(max(Y3)),ops.ops_loc{3}.my));

U4 = U4(1:end-1,1:end-1);
V4 = V4(1:end-1,1:end-1);
X4 = X4(1:end-1,1:end-1);
Y4 = Y4(1:end-1,1:end-1);

U1 = U1(:,1:end-1);
V1 = V1(:,1:end-1);
X1 = X1(:,1:end-1);
Y1 = Y1(:,1:end-1);

U3 = U3(1:end-1,:);
V3 = V3(1:end-1,:);
X3 = X3(1:end-1,:);
Y3 = Y3(1:end-1,:);

U = [U4,U3;U1,U2];
V = [V4,V3;V1,V2];
X = [X4,X3;X1,X2];
Y = [Y4,Y3;Y1,Y2];

U = U(:,X(1,:) <= lims(2));
V = V(:,X(1,:) <= lims(2));
X = X(:,X(1,:) <= lims(2));
Y = Y(:,X(1,:) <= lims(2));

figure
hold on
streamlines.even_stream_arrow(X,Y,U,V,2,3);
xlabel('x')
ylabel('y')
title('Streamlines of velocity field')
axis(lims)
drawnow