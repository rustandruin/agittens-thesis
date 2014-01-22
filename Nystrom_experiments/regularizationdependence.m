% see how the error behaves as the regularization rho is varied
clear
clc

matlabpool open 4;

n = 500;
k = 20;
l = 200; % fixed number of column samples

residual = 10^(-10);
Sigma = diag([logspace(log10(1), log10(residual), 2*k), residual*ones(1, n-2*k)]);
opterr = Sigma(k+1, k+1);

numrhos = 15;
%rholist = logspace(log10(opterr), log10(1), numrhos);
rholist = logspace(-15,-5,numrhos);

tau = k/n; % maximally incoherent

[U1, truetau] = mtxGenMethod1(n,k,tau); % make U1 dense
[R,~] = eig(randn(n-k, n-k));
U2 = orth(eye(n) - U1*U1')*R; % the random rotation prevents U2 from being axis aligned
A = [U1 U2]*Sigma*[U1 U2]';
%[U, ~] = eig(randn(n));
%U = 1/sqrt(n)*fft(eye(n));
%A = U*Sigma*U';

numreps = 60;

for rhoidx = 1:numrhos
    
    rho = rholist(rhoidx);
    
    parfor repeatidx = 1:numreps
        fprintf('On rho %d/%d rep %d/%d\n', rhoidx, numrhos, repeatidx, numreps);
        indices = randperm(n);
        indices = indices(1:l);
        
        [simpleerr(rhoidx, repeatidx), ~] = ...
            simpleNystromExtension(A, indices);
        
        [truncatedWerr(rhoidx, repeatidx), ~] = ...
           truncatedWNystromExtension(A, rho, indices);
   
        [regularizedAerr(rhoidx, repeatidx), ~] = ...
          regularizedANystromExtension(A, rho, indices);
        
        [regularizedWerr(rhoidx, repeatidx), ~] = ...
            regularizedWNystromExtension(A, rho, indices);
    end

    meansimpleerr(rhoidx) = max(simpleerr(rhoidx, :));
    meantruncatedWerr(rhoidx) = max(truncatedWerr(rhoidx, :));
    meanregularizedAerr(rhoidx) = max(regularizedAerr(rhoidx, :));
    meanregularizedWerr(rhoidx) = max(regularizedWerr(rhoidx, :));
end

matlabpool close

% Visualization

ms = 8;
lw = 2;

style1 = '+-';
style2 = 'd-';
style3 = 'o-';
style4 = 'v-';

color1 = [.1 .1 .1];
color2 = [.5 .5 .5];
color3 = [.3 .3 .3];
color4 = [.4 .4 .4];

loglog(rholist, meansimpleerr/opterr, style1, 'Color', color1, 'MarkerFaceColor', color1, 'MarkerSize', ms, 'LineWidth', lw);
hold on;
loglog(rholist, meanregularizedAerr/opterr, style3, 'Color', color3, 'MarkerFaceColor', color3, 'MarkerSize', ms, 'LineWidth', lw);
loglog(rholist, meantruncatedWerr/opterr, style2, 'Color', color2, 'MarkerFaceColor', color2, 'MarkerSize', ms, 'LineWidth', lw);
%loglog(rholist, meanregularizedWerr/opterr, style4, 'Color', color4, 'MarkerFaceColor', color4, 'MarkerSize', ms, 'LineWidth', lw);
hold off;
xlab = xlabel('$\rho$ (regularization parameter)');
set(xlab, 'Interpreter', 'Latex');
yl = ylabel('$\|\mathbf{A} - \mathbf{C}\mathbf{W}^\dagger\mathbf{C}^T\|_2/\lambda_{11}(\mathbf{A})$');
set(yl, 'Interpreter', 'Latex');
legend('Nystrom', 'Algorithm 6.1', 'Algorithm 6.2', 'Location', 'NorthWest');
%legend('simple Nystrom', '[CD]', 'Algorithm 1', 'Algorithm 2', 'Location', 'NorthWest');

printcf('regularizationdependence.pdf', 12, 8, 6);