function [err, timing] = simpleNystromExtension(A, indices)
    tic;
    C = A(:, indices);
    W = C(indices,:);
    Y = pinv(W)*C';
    
    err = norm(A - C*Y); 
    timing = toc;
end

