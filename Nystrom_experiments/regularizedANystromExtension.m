function [err, timing] = regularizedANystromExtension(A, rho, indices)
    tic;
    E = eye(size(A));
    C = A(:, indices) + rho*E(:, indices);
    W = C(indices,:);
    Wrho = pinv(W,0);
    timing = toc;
    
    err = norm(A - C*Wrho*C');
end

