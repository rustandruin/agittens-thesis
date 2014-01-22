% Compare the performance of SRHT-based low-rank approximation to that 
% guaranteed by the bounds in "Improved matrix algorithms via the
% Subsampled Randomized Hadamard Transform", and for reference, against
% the performance of Gaussian-based low-rank approximations

%% Clear
clear; clc; close all

%% Set up test matrix

numrepeats = 30; % number of times to rerun the low-rank approximation for each value of r

p = 10; 
n = 2^p; % the number of columns in the test matrix (power of 2s only)

alpha = .1;
M = 10;
A = alpha*eye(n+1);
A = A(:,2:end);
e1 = [1; zeros(n,1)];
A = M*(A + repmat(e1,1,n));

[U,S,V] = svd(A);
s = diag(S)';
fnorms = sqrt(cumsum(s.^2));

rho = n; % rank of A
opt_fnorm_errs = [fliplr(sqrt(cumsum(fliplr(s(2:end).^2)))) 0];
opt_snorm_errs = [s(2:end),0];

%% Experiment 1: Calculate approx errors as function of number of col samples, 
% for both SRHT and Gaussian methods.

kmax = floor(n/(2*log(n))); % choose so 2*kmax*log(n) < n
numpts = 30;
kgrid = ceil(linspace(2,kmax,numpts));
rgrid = ceil(2*kgrid*log(n));

% do the most time consuming experiments first so I can determine if I need
% to drop the size of the experiments
for idx = length(rgrid):-1:1
    fprintf('On trial %d of %d\n', length(rgrid)-idx+1, length(rgrid)); 
    r = rgrid(idx);
    k = kgrid(idx);
    
    Ak = U(:, 1:k)*S(1:k, 1:k)*V(:, 1:k)';
    
    tic
    [srht_resid_spec1(:, idx), srht_resid_spec2(:, idx), ...
     srht_forward_spec1(:, idx), srht_forward_spec2(:, idx), ...
     srht_resid_frob1(:, idx), srht_resid_frob2(:, idx), ...
     srht_forward_frob1(:, idx), srht_forward_frob2(:, idx)] = ...
        srhtapprox(A, Ak, k, r, numrepeats);
    
    [gauss_resid_spec1(:, idx), gauss_resid_spec2(:, idx), ...
     gauss_forward_spec1(:, idx), gauss_forward_spec2(:, idx), ...
     gauss_resid_frob1(:, idx), gauss_resid_frob2(:, idx), ...
     gauss_forward_frob1(:, idx), gauss_forward_frob2(:, idx)] = ...
        gaussianapprox(A, Ak, k, r, numrepeats);
    
    fprintf('... took %d seconds\n', toc);
end

%% Compute means and optimal errors for plotting

opt_resid_spec_errs = opt_snorm_errs(kgrid);
opt_forward_spec_errs = s(kgrid);
opt_resid_frob_errs = opt_fnorm_errs(kgrid);
opt_forward_frob_errs = fnorms(kgrid);

meansrht_resid_spec1 = mean(srht_resid_spec1);
meansrht_resid_spec2 = mean(srht_resid_spec2);
meansrht_resid_frob1 = mean(srht_resid_frob1);
meansrht_resid_frob2 = mean(srht_resid_frob2);

meangauss_resid_spec1 = mean(gauss_resid_spec1);
meangauss_resid_spec2 = mean(gauss_resid_spec2);
meangauss_resid_frob1 = mean(gauss_resid_frob1);
meangauss_resid_frob2 = mean(gauss_resid_frob2);

meansrht_forward_spec1 = mean(srht_forward_spec1);
meansrht_forward_spec2 = mean(srht_forward_spec2);
meansrht_forward_frob1 = mean(srht_forward_frob1);
meansrht_forward_frob2 = mean(srht_forward_frob2);

meangauss_forward_spec1 = mean(gauss_forward_spec1);
meangauss_forward_spec2 = mean(gauss_forward_spec2);
meangauss_forward_frob1 = mean(gauss_forward_frob1);
meangauss_forward_frob2 = mean(gauss_forward_frob2);

%% Visualize errors

myplot = @plot; % easily facilitate changing plot types

lwbottom = 2;
lwtop = 2;
mksize = 10;

srs1_style = 'ko';
srs2_style = 'ms';
grs1_style = 'r+';
grs2_style = 'b.';

sfs1_style = 'ko';
sfs2_style = 'ms';
gfs1_style = 'r+';
gfs2_style = 'b.';

srf1_style = 'ko';
srf2_style = 'ms';
grf1_style = 'r+';
grf2_style = 'b.';

sff1_style = 'ko';
sff2_style = 'ms';
gff1_style = 'r+';
gff2_style = 'b.';

% Plot relative residual spectral errors for both
% gaussian and srht methods. 

figure(1);
myplot(kgrid, meangauss_resid_spec1./opt_resid_spec_errs, grs1_style, 'LineWidth', lwbottom, 'MarkerSize', mksize);
hold on;
myplot(kgrid, meansrht_resid_spec1./opt_resid_spec_errs, srs1_style, 'LineWidth', lwtop, 'MarkerSize', mksize);

%myplot(kgrid, meangauss_resid_spec2./opt_resid_spec_errs, grs2_style, 'LineWidth', lwbottom, 'MarkerSize', mksize);
%myplot(kgrid, meansrht_resid_spec2./opt_resid_spec_errs, srs2_style, 'LineWidth', lwtop, 'MarkerSize', mksize);

hold off;
xlabel('k (target rank)');
yh = ylabel('relative error');
th = title('Relative residual spectral errors for $\mathbf{A}$');
set(th, 'interpreter', 'latex');
lh = legend('$\|\mathbf{A} - \mathbf{P}_{\mathbf{A}\mathbf{S}} \mathbf{A} \|_2 / \|\mathbf{A} - \mathbf{A}_k\|_2$, Gaussian', ...
            '$\|\mathbf{A} - \mathbf{P}_{\mathbf{A}\mathbf{S}} \mathbf{A} \|_2/ \|\mathbf{A} - \mathbf{A}_k\|_2$, SRHT');            
%lh = legend('$\|\mathbf{A} - \mathbf{P}_{\mathbf{A}\mathbf{S}} \mathbf{A} \|_2 / \|\mathbf{A} - \mathbf{A}_k\|_2$, Gaussian', ...
%            '$\|\mathbf{A} - \mathbf{P}_{\mathbf{A}\mathbf{S}} \mathbf{A} \|_2/ \|\mathbf{A} - \mathbf{A}_k\|_2$, SRHT', ...
%            '$\|\mathbf{A} - \mathbf{\Pi}_{\mathbf{A} \mathbf{S}, k}^{\mathrm{F}}(\mathbf{A})\|_2/\|\mathbf{A} - \mathbf{A}_k\|_2$, Gaussian', ...
%            '$\|\mathbf{A} - \mathbf{\Pi}_{\mathbf{A} \mathbf{S}, k}^{\mathrm{F}}(\mathbf{A})\|_2/\|\mathbf{A} - \mathbf{A}_k\|_2$, SRHT')
set(lh, 'interpreter', 'latex');

printcf('experimentA-residual-spectral.pdf', 16, 7, 6);

% Plot corresponding Frobenius relative errors

figure(2);
myplot(kgrid, meangauss_resid_frob1./opt_resid_frob_errs, grf1_style, 'LineWidth', lwbottom, 'MarkerSize', mksize);
hold on;
myplot(kgrid, meansrht_resid_frob1./opt_resid_frob_errs, srf1_style, 'LineWidth', lwtop, 'MarkerSize', mksize);

%myplot(kgrid, meangauss_resid_frob2./opt_resid_frob_errs, grf2_style, 'LineWidth', lwbottom, 'MarkerSize', mksize);
%myplot(kgrid, meansrht_resid_frob2./opt_resid_frob_errs, srf2_style, 'LineWidth', lwtop, 'MarkerSize', mksize);

hold off;
xlabel('k (target rank)');
yh = ylabel('relative error');
th = title('Relative residual Frobenius errors for $\mathbf{A}$');
set(th, 'interpreter', 'latex');
lh = legend('$\|\mathbf{A} - \mathbf{P}_{\mathbf{A}\mathbf{S}} \mathbf{A} \|_F / \|\mathbf{A} - \mathbf{A}_k\|_F$, Gaussian', ...
             '$\|\mathbf{A} - \mathbf{P}_{\mathbf{A}\mathbf{S}} \mathbf{A}\|_F / \|\mathbf{A} - \mathbf{A}_k\|_F$, SRHT', ...
             'Location', 'SouthWest');
% lh = legend('$\|\mathbf{A} - \mathbf{P}_{\mathbf{A}\mathbf{S}} \mathbf{A} \|_F / \|\mathbf{A} - \mathbf{A}_k\|_F$, Gaussian', ...
%             '$\|\mathbf{A} - \mathbf{P}_{\mathbf{A}\mathbf{S}} \mathbf{A} \|_F / \|\mathbf{A} - \mathbf{A}_k\|_F$, SRHT', ...
%             '$\|\mathbf{A} - \mathbf{\Pi}_{\mathbf{A} \mathbf{S}, k}^{\mathrm{F}}(\mathbf{A})\|_F / \|\mathbf{A} - \mathbf{A}_k\|_F$, Gaussian', ...
%             '$\|\mathbf{A} - \mathbf{\Pi}_{\mathbf{A} \mathbf{S}, k}^{\mathrm{F}}(\mathbf{A})\|_F / \|\mathbf{A} - \mathbf{A}_k\|_F$, SRHT', 'Location', 'SouthWest');
set(lh, 'interpreter', 'latex');

printcf('experimentA-residual-frobenius.pdf', 16, 7, 6);

% Plot forward spectral errors

figure(3);
myplot(kgrid, meangauss_forward_spec1./opt_forward_spec_errs, gfs1_style, 'LineWidth', lwbottom, 'MarkerSize', mksize);
hold on;
myplot(kgrid, meansrht_forward_spec1./opt_forward_spec_errs, sfs1_style, 'LineWidth', lwtop, 'MarkerSize', mksize);

myplot(kgrid, meangauss_forward_spec2./opt_forward_spec_errs, gfs2_style, 'LineWidth', lwbottom, 'MarkerSize', mksize);
myplot(kgrid, meansrht_forward_spec2./opt_forward_spec_errs, sfs2_style, 'LineWidth', lwtop, 'MarkerSize', mksize);

hold off;
xlabel('k (target rank)');
yh = ylabel('relative error');
th = title('Relative forward spectral errors for $\mathbf{A}$');
set(th, 'interpreter', 'latex');
lh = legend('$\|\mathbf{A}_k - \mathbf{P}_{\mathbf{A}\mathbf{S}} \mathbf{A} \|_2 / \|\mathbf{A}_k\|_2$, Gaussian', ...
            '$\|\mathbf{A}_k - \mathbf{P}_{\mathbf{A}\mathbf{S}} \mathbf{A} \|_2/ \|\mathbf{A}_k\|_2$, SRHT', ...
            '$\|\mathbf{A}_k - \mathbf{\Pi}_{\mathbf{A} \mathbf{S}, k}^{\mathrm{F}}(\mathbf{A})\|_2/\|\mathbf{A}_k\|_2$, Gaussian', ...
            '$\|\mathbf{A}_k - \mathbf{\Pi}_{\mathbf{A} \mathbf{S}, k}^{\mathrm{F}}(\mathbf{A})\|_2/\|\mathbf{A}_k\|_2$, SRHT');
set(lh, 'interpreter', 'latex');

printcf('experimentA-forward-spectral.pdf', 16, 7, 6);

% Plot forward Frobenius errors

figure(4);
myplot(kgrid, meangauss_forward_frob1./opt_forward_frob_errs, gff1_style, 'LineWidth', lwbottom, 'MarkerSize', mksize);
hold on;
myplot(kgrid, meansrht_forward_frob1./opt_forward_frob_errs, sff1_style, 'LineWidth', lwtop, 'MarkerSize', mksize);

myplot(kgrid, meangauss_forward_frob2./opt_forward_frob_errs, gff2_style, 'LineWidth', lwbottom, 'MarkerSize', mksize);
myplot(kgrid, meansrht_forward_frob2./opt_forward_frob_errs, sff2_style, 'LineWidth', lwtop, 'MarkerSize', mksize);

hold off;
xlabel('k (target rank)');
yh = ylabel('relative error');
th = title('Relative forward Frobenius errors for $\mathbf{A}$');
set(th, 'interpreter', 'latex');
lh = legend('$\|\mathbf{A}_k - \mathbf{P}_{\mathbf{A}\mathbf{S}} \mathbf{A} \|_F / \|\mathbf{A}_k\|_F$, Gaussian', ...
            '$\|\mathbf{A}_k - \mathbf{P}_{\mathbf{A}\mathbf{S}} \mathbf{A} \|_F/ \|\mathbf{A}_k\|_F$, SRHT', ...
            '$\|\mathbf{A}_k - \mathbf{\Pi}_{\mathbf{A} \mathbf{S}, k}^{\mathrm{F}}(\mathbf{A})\|_F/\|\mathbf{A}_k\|_F$, Gaussian', ...
            '$\|\mathbf{A}_k - \mathbf{\Pi}_{\mathbf{A} \mathbf{S}, k}^{\mathrm{F}}(\mathbf{A})\|_F/\|\mathbf{A}_k\|_F$, SRHT', 'Location', 'NorthWest');
set(lh, 'interpreter', 'latex');


printcf('experimentA-forward-frobenius.pdf', 16, 7, 6);

%% Compare the empirically necessary number of samples to that predicted by our bounds
%% for relative error frobenius bound

% The lower bound on r necessary for the guarantees of the theorem to hold,
% ignoring constants
sufficientr = @(k, vareps, delta) ...
    ceil(vareps.^(-1) * (sqrt(k) + sqrt(log(n/delta))).^2 .* log(k/delta));

kgrid2 = 1:20;
eps1 = 1/2;
eps2 = 1/3;
eps3 = 1/6;
delta1 = 1/2;

rest_eps1 = sufficientr(kgrid2, eps1, delta1);
rest_eps2 = sufficientr(kgrid2, eps2, delta1);
rest_eps3 = sufficientr(kgrid2, eps3, delta1);

% calculate the empirically necessary r to achieve (1+eps) frobenius error 
% guarantees

% for approx 1
[remp_eps1_rf1, maxk_eps1_rf1] = invertDecreasingCurve(srht_resid_frob1, opt_fnorm_errs*(1+eps1), delta1);
[remp_eps1_ff1, maxk_eps1_ff1] = invertDecreasingCurve(srht_forward_frob1, opt_fnorm_errs*(1+eps1), delta1);
remp_eps1_rf1 = rgrid(remp_eps1_rf1);
remp_eps1_ff1 = rgrid(remp_eps1_ff1);

[remp_eps2_rf1, maxk_eps2_rf1] = invertDecreasingCurve(srht_resid_frob1, opt_fnorm_errs*(1+eps2), delta1);
[remp_eps2_ff1, maxk_eps2_ff1] = invertDecreasingCurve(srht_forward_frob1, opt_fnorm_errs*(1+eps2), delta1);
remp_eps2_rf1 = rgrid(remp_eps2_rf1);
remp_eps2_ff1 = rgrid(remp_eps2_ff1);

[remp_eps3_rf1, maxk_eps3_rf1] = invertDecreasingCurve(srht_resid_frob1, opt_fnorm_errs*(1+eps3), delta1);
[remp_eps3_ff1, maxk_eps3_ff1] = invertDecreasingCurve(srht_forward_frob1, opt_fnorm_errs*(1+eps3), delta1);
remp_eps3_rf1 = rgrid(remp_eps3_rf1);
remp_eps3_ff1 = rgrid(remp_eps3_ff1);

% for approx 2
[remp_eps1_rf2, maxk_eps1_rf2] = invertDecreasingCurve(srht_resid_frob2, opt_fnorm_errs*(1+eps1), delta1);
[remp_eps1_ff2, maxk_eps1_ff2] = invertDecreasingCurve(srht_forward_frob2, opt_fnorm_errs*(1+eps1), delta1);
remp_eps1_rf2 = rgrid(remp_eps1_rf2);
remp_eps1_ff2 = rgrid(remp_eps1_ff2);

[remp_eps2_rf2, maxk_eps2_rf2] = invertDecreasingCurve(srht_resid_frob2, opt_fnorm_errs*(1+eps2), delta1);
[remp_eps2_ff2, maxk_eps2_ff2] = invertDecreasingCurve(srht_forward_frob2, opt_fnorm_errs*(1+eps2), delta1);
remp_eps2_rf2 = rgrid(remp_eps2_rf2);
remp_eps2_ff2 = rgrid(remp_eps2_ff2);

[remp_eps3_rf2, maxk_eps3_rf2] = invertDecreasingCurve(srht_resid_frob2, opt_fnorm_errs*(1+eps3), delta1);
[remp_eps3_ff2, maxk_eps3_ff2] = invertDecreasingCurve(srht_forward_frob2, opt_fnorm_errs*(1+eps3), delta1);
remp_eps3_rf2 = rgrid(remp_eps3_rf2);
remp_eps3_ff2 = rgrid(remp_eps3_ff2);

% this is the highest k for which there exists a relative error
% approximation for both residual and forward frobenius errors, for both
% approximation methods, at all values of eps
maxk = max([maxk_eps1_rf1 maxk_eps1_ff1 maxk_eps2_rf1 maxk_eps2_ff1 ...
            maxk_eps3_rf1 maxk_eps3_ff1 maxk_eps1_rf2 maxk_eps1_ff2 ...
            maxk_eps2_rf2 maxk_eps2_ff2 maxk_eps3_rf2 maxk_eps3_ff2]);
        
%% Visualize the sufficient bound on r for (1+eps) Frobenius approximations
% to hold versus the empirically sufficient bound on r 

myplot = @plot; % conveniently change type of plot

% the data points for the YY^\pinv A approx are too dense, so subsample
% them
downsamplerate = 7;
krange1 = downsample(1:maxk_eps1_rf1, downsamplerate);
remp_eps1_rf1_downsampled = downsample(remp_eps1_rf1, downsamplerate);
krange2 = downsample(1:maxk_eps2_rf1, downsamplerate);
remp_eps2_rf1_downsampled = downsample(remp_eps2_rf1, downsamplerate);
krange3 = downsample(1:maxk_eps3_rf1, downsamplerate);
remp_eps3_rf1_downsampled = downsample(remp_eps3_rf1, downsamplerate);

myplot(krange1, remp_eps1_rf1_downsampled, 'k.', 'LineWidth', 1);
hold on;
myplot(krange2, remp_eps2_rf1_downsampled, 'k*', 'LineWidth', 1);
myplot(krange3, remp_eps3_rf1_downsampled, 'k+', 'LineWidth', 1);

myplot(1:maxk_eps1_rf2, remp_eps1_rf2, 'd', 'LineWidth', 1);
myplot(1:maxk_eps2_rf2, remp_eps2_rf2, 's', 'LineWidth', 1);
myplot(1:maxk_eps3_rf2, remp_eps3_rf2, 'o', 'LineWidth', 1);
hold off;

xlabel('k (target rank)')
yh = ylabel('$\ell$ (column samples)');
th = title({'Number of samples needed to achieve $\|\mathbf{A} - \mathbf{L}\|_F \leq (1+\varepsilon) \|\mathbf{A} - \mathbf{A}_k\|_F$ with probability at least 1/2'});
set(th, 'interpreter', 'latex');
set(yh, 'interpreter', 'latex');
lh = legend('$\mathbf{L} = \mathbf{P}_{\mathbf{A}\mathbf{S}} \mathbf{A}, \varepsilon = 1/2$', ...
            '$\mathbf{L} = \mathbf{P}_{\mathbf{A}\mathbf{S}} \mathbf{A}, \varepsilon = 1/3$', ...
            '$\mathbf{L} = \mathbf{P}_{\mathbf{A}\mathbf{S}} \mathbf{A}, \varepsilon = 1/6$', ...
            '$\mathbf{L} = \mathbf{\Pi}_{\mathbf{A} \mathbf{S}, k}^{\mathrm{F}}(\mathbf{A}), \varepsilon = 1/2$', ...
            '$\mathbf{L} = \mathbf{\Pi}_{\mathbf{A} \mathbf{S}, k}^{\mathrm{F}}(\mathbf{A}), \varepsilon = 1/3$', ...
            '$\mathbf{L} = \mathbf{\Pi}_{\mathbf{A} \mathbf{S}, k}^{\mathrm{F}}(\mathbf{A}), \varepsilon = 1/6$', ...
            'Location', 'SouthEast');
set(lh, 'interpreter', 'latex');

printcf('experimentA-empirically-necessary-r.pdf', 16, 12, 6);

%% For some values of k, determine the lower bound on r necessary for the 
% theorem to apply, calculate the upper bound on the spectral errors, and
% compare to the empirical spectral errors

eps4 = 1/2;
delta2 = 1/2;
kgrid3 = 1:20;

rgrid3 = sufficientr(kgrid3, eps4, delta2);

% find an upper bound on the spectral errors
predictedspecerr = @(r, optspec, optfroerr, rho, delta) ...
     (1 + sqrt(log(n/delta)*log(rho/delta)./r)).*optspec + ...
     sqrt(log(rho/delta)./r).*optfroerr;

serr_pred = predictedspecerr(rgrid3, opt_snorm_errs(kgrid3), ...
                             opt_fnorm_errs(kgrid3), rho, delta2);

% calculate the actual spectral errors
for idx = length(rgrid3):-1:1
    fprintf('On trial %d of %d\n', length(rgrid3)-idx+1, length(rgrid3)); 
    r = rgrid3(idx);
    k = kgrid3(idx);
    
    Ak = U(:, 1:k)*S(1:k, 1:k)*V(:, 1:k)';
    
    tic
    [emp_sr1(:, idx), emp_sr2(:, idx), emp_sf1(:, idx), emp_sf2(:, idx),...
       ~, ~, ~, ~] = srhtapprox(A, Ak, k, r, numrepeats);
    fprintf('... took %d seconds\n', toc);
end

%% Visualize the actual spectral errors vs the guaranteed upper bound

myplot = @plot; % conveniently change type of plot

% residuals
figure();

myplot(kgrid3, serr_pred./opt_snorm_errs(kgrid3), 'ko-', 'LineWidth', 2);
hold on;
myplot(kgrid3, mean(emp_sr1)./opt_snorm_errs(kgrid3), 'r+-', 'LineWidth', 2);
myplot(kgrid3, mean(emp_sr2)./opt_snorm_errs(kgrid3), 'bs-', 'LineWidth', 2);
xlabel('k (target rank)')
ylabel('relative error')
th = title('Guaranteed vs observed spectral errors for $\mathbf{A}$');
set(th, 'interpreter', 'latex');
lh = legend('spectral error bound / $\|\mathbf{A} - \mathbf{A_k}\|_2$', ...
            '$\|\mathbf{A} - \mathbf{P}_{\mathbf{A}\mathbf{S}} \mathbf{A}\|_2 / \|\mathbf{A} - \mathbf{A_k}\|_2$',...
            '$\|\mathbf{A} - \mathbf{\Pi}_{\mathbf{A} \mathbf{S}, k}^{\mathrm{F}}(\mathbf{A})\|_2 / \|\mathbf{A} - \mathbf{A}_k\|_2$');
set(lh, 'interpreter', 'latex');

printcf('experimentA-actual-versus-predicted-spectral-error.pdf', 16, 6, 6);

%% Save data so the visualization parts can be run without the data generation parts

save('experimentAdata.mat');