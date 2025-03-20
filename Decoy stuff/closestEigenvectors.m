%{
FUNCTION NAME: closestEigenvectors
 This function computes the maximum distance between
 the eigenvectors of two matrices that are close in 1-norm.

 Parameters: 

 * rho - The known finite dimensional state.

 * ep - The maximum distance between the known and unknown state in
       1-norm.

* n - Ordering the eigenvalues from largest to smallest, the eigenvector
      corresponding to the nth eigenvalue is the one whose distance we wish to
      find. For eg.- n=1 means to find the max dist between the first
      eigenvector of rho to the corresponding eigenvector in the matrix
      at most ep away from rho (in 1-norm).

Return:
 
* vec - The eigenvector corresponding to the nth eigenvalue of rho.

* eigvalue - The nth eigenvalue of rho.

* dist - The maxumum 1-norm of the difference of the two eigenvectors.
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vec,eigvalue,dist] = closestEigenvectors(rho, ep, n)
    
    [V,D] = eig(rho);
    [~,ind] = sort(diag(D),'descend');
    D = D(ind,ind);
    V = V(:,ind);
    e = diag(D);
    vec = V(:,n);
    eigvalue = e(n);
    if n == 1
        delta = e(1)-e(2)-ep;
    elseif n == length(rho)
        delta = e(n-1)-e(n)-ep;
    else
        delta = min(e(n-1)-e(n)-ep, e(n)-e(n+1)-ep);
    end
    
    if delta<0 || ep>delta
        dist = 2;
    else
        dist = 2*ep/delta;
    end
    if dist > 0 && dist <eps
        dist = eps
    end

end
