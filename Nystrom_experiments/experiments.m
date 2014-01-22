%% Experiments with PSD matrices which have specified coherence:
% - is the coherence a reasonable indicator of whether the Nystrom
%   procedure is feasible?
% - how many measurements are needed as a function of the coherence, for
%   the Nystrom procedure to give good results?
% - how do the different methods compare?

%% Clear workspace (because will save all variables later) and add the
%% matrix generation function to the path
clear;
cd kappaSQ3.20.12
path(path, pwd)
cd ..

%% Experiment parameters: the size of the matrix, the spectrum, the number
%% of trials to run, the coherence, the regularization parameter, etc.
n = 1000; % size of the matrix
p = -.8; % exponent in the power law decay
C = 1000; % maximum eigenvalue
spectrum = C*[1:n].^p; % list of the eigenvalues
Sigma = diag(spectrum); % diagonal matrix of eigenvalues

rho = 10^-1; % parameter for truncation and Tikhonov regularization
k = 20; % target rank of the approximation
maxsamples = n; % maximum number of columns to use in Nystrom approxs
% number of sampling parameters between k and maxsamples to use to generate
% plots
numells = 10; 
% number of coherences between 1 and n/k to use to generate plots
numcoherences = 3; 
% list of sampling parameters ranging from k to maxsamples
elllist = ceil([k:(maxsamples-k)/(numells-1):maxsamples]); 
% list of taus ranging from k/n to 1 (note that tau = k/n*mu is the maximum
% diagonal entry in the projection onto the target rank-k space, not the 
% coherence (that is mu))
taulist = [k/n:(1-k/n)/(numcoherences-1):1]; 
numtrials = 20; % number of trials to run per sample

% for quick runs
numcoherences = 1;
taulist = k/n;
numtrials = 5;

%% Run the trials

matlabpool open

% keep time measurement so I'll have a rough idea of what feasible
% parameter settings are
starttime = cputime;

% the relative spectral error as a function of the coherence, l, the number
% of columns sampled, and the current trial
simplenysrelerr = zeros(numtrials, numcoherences, numells); 
regularizedWrelerr = zeros(numtrials, numcoherences, numells);
regularizedArelerr = zeros(numtrials, numcoherences, numells);
truncatedWrelerr = zeros(numtrials, numcoherences, numells);
rrqrrelerr = zeros(numtrials, numcoherences, numells);

% the cost of all these computations
simplenyscost = zeros(numtrials, numcoherences, numells);
regularizedWcost = zeros(numtrials, numcoherences, numells);
regularizedAcost = zeros(numtrials, numcoherences, numells);
truncatedWcost = zeros(numtrials, numcoherences, numells);
rrqrcost = zeros(numtrials, numcoherences, numells);

parfor trialno = 1:numtrials
    
    cursimpletrialerror = zeros(numcoherences, numells);
    curregularizedWtrialerror = zeros(numcoherences, numells);
    curregularizedAtrialerror = zeros(numcoherences, numells);
    curtruncatedWtrialerror = zeros(numcoherences, numells);
    currrqrtrialerror = zeros(numcoherences, numells);
    
    cursimpletrialcost = zeros(numcoherences, numells);
    curregularizedWtrialcost = zeros(numcoherences, numells);
    curregularizedAtrialcost = zeros(numcoherences, numells);
    curtruncatedWtrialcost = zeros(numcoherences, numells);
    currrqrtrialcost = zeros(numcoherences, numells);
    
    for coherenceno = 1:numcoherences
        for ellno = 1:numells
            
            % generate a matrix with the given spectrum and top k-dim 
            % invariant subspace with the specified coherence
            % mtxGenMethod3 makes as many columns zero as is possible, so 
            % that e.g. tau=1 ensures U1 consists of columns from the 
            % identity matrix
            % NOTE: if only set coherence to a certain value by making one 
            % col have that norm, as opposed to making as many cols as 
            % possible have that norm and the remaining be zeroes, then you
            % can use much fewer samples to get low error, so maybe max
            % leverage is too restrictive: look for bounds using average
            % leverage?
            [U1, truetau] = mtxGenMethod3(n,k,taulist(coherenceno)); 
            U2 = orth(eye(n) - U1*U1');
            A = [U1 U2]*Sigma*[U1 U2]';
            
            % REMOVE when needed!
            % set regularization parameter to minimize theoretical bounds 
            rho = spectrum(k+1)*elllist(ellno)/(2*n);
            
            % compute the relative spectral norm errors of the various
            % extensions
            l = elllist(ellno);
            indices = randperm(n);
            indices = indices(1:l);
            
            % simple Nystrom extension
            [err, timing] = simpleNystromExtension(A, indices);
            cursimpletrialerror(coherenceno, ellno) = err/spectrum(k+1);
            cursimpletrialcost(coherenceno, ellno) = timing;
            
            % truncated W Nystrom extension
            [err, timing] = truncatedWNystromExtension(A, rho, indices);
            curtruncatedWtrialerror(coherenceno, ellno) = err/spectrum(k+1);
            curtruncatedWtrialcost(coherenceno, ellno) = timing;
            
            % regularized coupling matrix extension
            [err, timing] = regularizedWNystromExtension(A, rho, indices);
            curregularizedWtrialerror(coherenceno, ellno) = err/spectrum(k+1);
            curregularizedWtrialcost(coherenceno, ellno) = timing;
            
            % regularized A Nystrom extension
            [err, timing] = regularizedANystromExtension(A, rho, indices);
            curregularizedAtrialerror(coherenceno, ellno) = err/spectrum(k+1);
            curregularizedAtrialcost(coherenceno, ellno) = timing;
            
            % RRQR Nystrom extension: select the columns from A
            % corresponding to the most informative columns of A^{1/2}
            [err, timing] = rrqrNystromExtension(A, l);
            currrqrtrialerror(coherenceno, ellno) = err/spectrum(k+1);
            currrqrtrialcost(coherenceno, ellno) = timing;
            
            % give a status update
            fprintf(['completed trial %d/%d, coherence %d/%d, ' ...
                'l-value %d/%d \n'], trialno, numtrials, coherenceno, ...
                numcoherences, ellno, numells);
        end
    end
    
    simplenysrelerr(trialno,:,:) = cursimpletrialerror;
    regularizedWrelerr(trialno,:,:) = curregularizedWtrialerror;
    regularizedArelerr(trialno,:,:) = curregularizedAtrialerror;
    truncatedWrelerr(trialno,:,:) = curtruncatedWtrialerror;
    rrqrrelerr(trialno,:,:) = currrqrtrialerror;
    
    simplenyscost(trialno,:,:) = cursimpletrialcost;
    regularizedWcost(trialno,:,:) = curregularizedWtrialcost;
    regularizedAcost(trialno,:,:) = curregularizedAtrialcost;
    truncatedWcost(trialno,:,:) = curtruncatedWtrialcost;
    rrqrcost(trialno,:,:) = currrqrtrialcost;
end

time_elapsed = cputime - starttime;

matlabpool close

%% Save data
save(strcat('varying-coherence-theoretical-rho',datestr(now, 'dd-mm-yyyy-HH-MM')));
 
%% Visualization
% generate side-by-side boxplots displaying the performance of each
% extension as a function of the number of samples
coherencenum = 1;

simplenystrom_error = ...
    squeeze(shiftdim(simplenysrelerr(:,coherencenum,:),1))';
regularizedW_error = ...
    squeeze(shiftdim(regularizedWrelerr(:,coherencenum,:),1))';
regularizedA_error = ...
    squeeze(shiftdim(regularizedArelerr(:,coherencenum,:),1))';
truncatedW_error = ...
    squeeze(shiftdim(truncatedWrelerr(:,coherencenum,:),1))';
rrqr_error = ...
    squeeze(shiftdim(rrqrrelerr(:,coherencenum,:),1))';

% boxplotCsub can't handle group specifications as matrices, so everything
% has to be vectorized
elllist_vec = repmat(elllist, numtrials, 1);
elllist_vec = elllist_vec(:);

l= [200:n];
plot(l, 2*n./l, 'k', 'LineWidth', 2);
hold on
plot(elllist, mean(rrqr_error), 'y-*')
plot(elllist, mean(simplenystrom_error), 'r-*')
plot(elllist, mean(regularizedW_error), 'g-*')
plot(elllist, mean(regularizedA_error), 'm-*')
plot(elllist, mean(truncatedW_error), 'b-*')
hold off

figure();
boxplotCsub(rrqr_error(:), elllist_vec, 0, ['+'], 1, 1.5, 'r', ...
    false, 0.5, true, [1 5] );
boxplotCsub(simplenystrom_error(:), elllist_vec, 0, ['+'], 1, 1.5, 'g', ...
    false, 0.5, true, [2 5] );
boxplotCsub(regularizedW_error(:), elllist_vec, 0, ['+'], 1, 1.5, 'm', ...
    false, 0.5, true, [3 5] );
boxplotCsub(regularizedA_error(:), elllist_vec, 0, ['+'], 1, 1.5, 'b', ...
    false, 0.5, true, [4 5] );
boxplotCsub(truncatedW_error(:), elllist_vec, 0, ['+'], 1, 1.5, 'y', ...
    false, 0.5, true, [5 5] );

ylabel('Relative spectral norm error')
xlabel('Number of columns sampled', 'interpreter', 'latex')
