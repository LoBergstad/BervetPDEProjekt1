function A = set2zero(A,tol,verbose)
if verbose
    nnz_before = nnz(A);
end
[i,j,vals] = find(A);
vals(abs(vals) < tol) = 0;
A = sparse(i,j,vals,size(A,1),size(A,2));

if verbose
    nnz_after = nnz(A);
    fprintf("nnz before: %.2e, nnz after: %.2e\n",nnz_before/numel(A),nnz_after/numel(A));
end