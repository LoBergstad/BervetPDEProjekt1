% Removes linearly dependent rows of Lold
function [Lnew,ids_remove] = remove_deps(Lold,H,TOL)
% clear
% Lold = rand(5,5);

% Norms smaller than TOL are assumed 0
% TOL = 1e-6;

% Gram-Schmidt
L = Lold';

inprod = @(u,v) u'*H*v;
norm = @(u) sqrt(inprod(u,u));

N = size(L,2);
Lort = sparse(size(L,1),size(L,2));
Lort(:,1) = L(:,1)/norm(L(:,1));
ids_remove = []; % indices of removed rows

% Gram-Schmidt
for i = 2:N
    tmp = inprod(L(:,i),Lort(:,1:i-1))';
    r = Lort(:,1:i-1)*tmp;
    v = L(:,i) - r;

    if norm(v) < TOL
        % collect indices of linearly dependent vectors in L
        ids_remove = [ids_remove;i];
    else
        Lort(:,i) = v/norm(v);
    end
end

Lnew = Lold;
Lnew(ids_remove,:) = [];