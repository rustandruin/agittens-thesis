function [err, timing] = truncatedWNystromExtension(A, rho, indices)
    tic;
        C = A(:, indices);
        W = C(indices,:);
    
        Wrho = pinv(W, rho);
    timing = toc;

    err = norm(A - C*Wrho*C');
end

