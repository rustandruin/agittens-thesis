% test if the multiplicative factor in the error bound is optimal in n and
% l
% Export to 4 by 6 pdf, then crop

n = 1000;

alpha = 1;
M = eye(n+1);
M = M(:, 2:(n+1)) + alpha*repmat([1; zeros(n,1)], 1, n);
A = M'*M;

k = 10;
dl = 20;
llist = dl:dl:(n-dl);
err = zeros(length(llist), 1);

for lidx = 1:length(llist)
    p = randperm(n);
    p = p(1:llist(lidx));
    [err(lidx), ~] = simpleNystromExtension(A, p);
end

optkerr = svds(A, k+1);
optkerr = optkerr(end);

myplot = @semilogy;

myplot(llist, n./llist, '-', 'Color', [.7 .7 .7], 'LineWidth', 2);
xlab = xlabel('$\ell$ (column samples)');
ylabel('error ratio');
set(xlab, 'Interpreter', 'Latex');
hold on;
myplot(llist, err/optkerr, 'k+', 'MarkerSize', 8, 'LineWidth', 2)
lh = legend('$\quad n/\ell$', '$\|\mathbf{A} - \mathbf{C}\mathbf{W}^\dagger\mathbf{C}^T\|_2/\lambda_{11}(\mathbf{A})$');
set(lh, 'Interpreter', 'Latex');
hold off;

printcf('optimalityverification.pdf', 12, 6, 4);
