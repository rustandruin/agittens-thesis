function [err, timing] = rrqrNystromExtension(A, l)
    tic;
        [V,D] = eig(A);
        [~,~,p,~] = rrqr(V*sqrt(D)*V');
        C = A(:, p(1:l));
        W = C(p(1:l), :);
    timing = toc;
    
    err = normest(A - C*pinv(W)*C'); 
end

