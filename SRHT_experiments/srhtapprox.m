function [residspecnorm1, residspecnorm2, ...
          forwardspecnorm1, forwardspecnorm2, ...
          residfrobnorm1, residfrobnorm2, ...
          forwardfrobnorm1, forwardfrobnorm2] = ...
                               srhtapprox(A, Ak, k, r, numiters)
%function [residspecnorm1, residspecnorm2, ...
%          forwardspecnorm1, forwardspecnorm2, ...
%          residfrobnorm1, residfrobnorm2, ...
%          forwardfrobnorm1, forwardfrobnorm2] = ...
%                               srhtapprox(A, Ak, k, r, numiters)
%
%   Forms two SRHT approximations of the matrix A numiters times,
%   each using r column samples, and returns the spectral and frobenius 
%   norm errors in vectors of length numiters. Ak is the optimal rank k
%   approximation to A.
%
%   Approximation 1 is YY^\dagger A, which may have rank > k,
%   Approximation 2 is QX_opt, which has rank at most k
%

n = size(A, 2);
if (ceil(log2(n)) - log2(n)) 
    error('The number of columns in the matrix should be a power of 2')
end

H = walsh(n);
parfor iter = 1:numiters
    d = 1 - 2*(rand(1,n)>1/2);
    cols = randperm(n);
    cols = cols(1:r);

    DH = bsxfun(@times, d', H);
    Y = 1/sqrt(r)*A*DH(:,cols);
    [Q,~] = qr(Y,0);
    QtA = Q'*A;
    
    approx1 = Q*QtA;
    
    [U,S,V] = svd(QtA);
    approx2 = Q*(U(:,1:k)*S(1:k, 1:k)*V(:, 1:k)');
    
    residspecnorm1(iter) = normest(A - approx1);
    residspecnorm2(iter) = normest(A - approx2);
    residfrobnorm1(iter) = norm(A - approx1, 'fro');
    residfrobnorm2(iter) = norm(A - approx2, 'fro');
    
    forwardspecnorm1(iter) = normest(Ak - approx1);
    forwardspecnorm2(iter) = normest(Ak - approx2);
    forwardfrobnorm1(iter) = norm(Ak - approx1, 'fro');
    forwardfrobnorm2(iter) = norm(Ak - approx2, 'fro');
end


end
