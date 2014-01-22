function [err, timing] = regularizedWNystromExtension(A, rho, indices)
    tic;
    C = A(:, indices);
    W = C(indices,:);
    
    s = sort(eig(W));
    if s(1) < rho
        W = W+rho*eye(size(W));
    end
    Wrho = pinv(W,0);
    
    timing = toc;
    
    err = norm(A - C*Wrho*C');
end

