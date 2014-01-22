% verify that coherence affects the error bounds: as coherence decreases,
% the error incurred by a fixed number of samples decreases
%
clear
clc

matlabpool open 4; % use parallelization if possible

showplots = true; % flag of whether to show plots
savedatasets = true; % flag of whether to save the generated datasets

datasetnames = {'truncatedWNystromExtensionDatasetSparseU1', ...
    'regularizedANystromExtensionDatasetSparseU1', ...
    'regularizedWNystromExtensionDatasetSparseU1', ...
    'truncatedWNystromExtensionDatasetDenseU1', ...
    'regularizedANystromExtensionDatasetDenseU1', ...
    'regularizedWNystromExtensionDatasetDenseU1'};
eigenvectortypes = ...
    {'sparse', 'sparse', 'sparse', 'dense', 'dense', 'dense'};
extensionfunctionlist = {@truncatedWNystromExtension, ...
    @regularizedANystromExtension, ...
    @regularizedWNystromExtension, ...
    @truncatedWNystromExtension, ...
    @regularizedANystromExtension, ...
    @regularizedWNystromExtension};

n = 500; % the size of the PSD matrix
k = 10; % the target rank we will compare approximation errors to
opterr = 10^(-15); % residual: the value for all the singular values above 2*k
Sigma = diag([logspace(1, -3, 2*k), opterr*ones(1,n-2*k)]);
opterr = Sigma(k+1,k+1);
rho = opterr; % optimal choice for regularization parameter

dl = 10; % increment size for l
llist = k:dl:(n-dl);
numreps = 60; % number of times to repeat each experiment

numcoherences = 4; % the number of values of coherences to use
% tau is the maximum diagonal entry in the projection
taulist = linspace(k/n, 1, numcoherences);

for datasetnum = 1:length(datasetnames)
    datasetname = datasetnames{datasetnum};
    eigenvectortype = eigenvectortypes{datasetnum};
    extensionfunction = extensionfunctionlist{datasetnum};
    
    for coherenceidx = 1:numcoherences
        % To generate the eigenvectors with specified sparsity type:
        % mtxGenMethod3 gives sparse eigenvectors
        % mtxGenMethod1 gives dense eigenvectors but hangs when tau is too high
        % mtxGenMethod2 gives dense eigenvectors (use instead of Method1 when
        % coherence is too high, but prefer Method1 as it gives denser results)
        
        switch eigenvectortype
            case 'sparse'
                [U1, truetau] = mtxGenMethod3(n,k,taulist(coherenceidx));
            case 'dense'
                if taulist(coherenceidx) < .99
                    [U1, truetau] = mtxGenMethod1(n,k,taulist(coherenceidx));
                else
                    [U1, truetau] = mtxGenMethod2(n,k,taulist(coherenceidx));
                end
        end
        
        % generate a U2 so that [U1 U2] is an onb, then form A
        % the random rotation prevents U2 from being too regular (e.g. when U1
        % has structure, don't want it to be passed through to U2)
        [R,~] = eig(randn(n-k, n-k));
        U2 = orth(eye(n) - U1*U1')*R;
        A = [U1 U2]*Sigma*[U1 U2]';
        
        % Generate the appropriate Nystrom extension and store statistics
        for lidx = 1:length(llist)
            fprintf(['Experiment ''%s'' (%d/%d), ', ...
                'coherence %d/%d, col samples %d/%d\n'], ...
                datasetname, datasetnum, length(datasetnames), ...
                coherenceidx, numcoherences, lidx, length(llist));
            
            % use parfor if available
            parfor repeatidx = 1:numreps
                %indices = randperm(n, llist(lidx));
                indices = randperm(n);
                indices = indices(1:llist(lidx));
                
                [err(coherenceidx, lidx, repeatidx), ...
                    timing(coherenceidx, lidx, repeatidx)] = ...
                    extensionfunction(A, rho, indices);
            end
            
            meanerr(coherenceidx, lidx) = mean(err(coherenceidx, lidx, :));
            stderr(coherenceidx, lidx) = std(err(coherenceidx, lidx, :));
            meantiming(coherenceidx, lidx) = mean(timing(coherenceidx, lidx, :));
        end
    end
    
    if showplots
        
        % Plot the results for this dataset: the mean error, the standard
        % deviation, and the timing
        % remember that mu = n/k * tau! The legenda need to be manually changed if
        % the parameters of the experiment do.
        
        myplot = @semilogy;
        figure(2*datasetnum-1);
        subplot(1,2,1);
        
        myplot(llist, meanerr(1,:)/opterr, 'r+-', 'LineWidth', 2);
        hold on;
        myplot(llist, meanerr(2,:)/opterr, 'go-', 'LineWidth', 2);
        myplot(llist, meanerr(3,:)/opterr, 'bd-', 'LineWidth', 2);
        myplot(llist, meanerr(4,:)/opterr, 'kv-', 'LineWidth', 2);
        title(datasetname);
        xh = xlabel('$\ell$ (columns samples)');
        yh = ylabel('$\|\mathbf{A} - \mathbf{C}\mathbf{W}^\dagger\mathbf{C}^\mathrm{T}\|_2 / \|\mathbf{A} - \mathbf{A}_k\|_2$');
        leg = legend('$$\mu= 1$$', '$\mu = 17$', '$\mu = 34$', '$\mu = 50$');
        set(xh, 'Interpreter', 'Latex');
        set(yh, 'Interpreter', 'Latex');
        set(leg, 'Interpreter', 'Latex');
        
        hold off;
        
        subplot(1,2,2);
        
        myplot(llist, stderr(1,:)/opterr, 'g+-', 'LineWidth', 2);
        hold on;
        myplot(llist, stderr(2,:)/opterr, 'ro-', 'LineWidth', 2);
        myplot(llist, stderr(3,:)/opterr, 'bd-', 'LineWidth', 2);
        myplot(llist, stderr(4,:)/opterr, 'kv-', 'LineWidth', 2);
        title(datasetname);
        xlabel('number of columns sampled');
        ylabel('std of relative spectral norm error');
        leg = legend('$$\mu= 1$$', '$\mu = 17$', '$\mu = 34$', '$\mu = 50$');
        set(leg, 'Interpreter', 'Latex')
        
        hold off;
        
        figure(2*datasetnum);
        myplot(llist, meantiming(1,:), 'g+-', 'LineWidth', 2);
        xlabel('number of columns sampled');
        ylabel('time to compute extension (s)');
        title('datasetname');
    end
    
    if savedatasets
        % Save the generated dataset for later visualization
        save([datasetname '.mat']);
    end
    
end

matlabpool close;
