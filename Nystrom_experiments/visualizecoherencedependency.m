%% parameters
ms = 8; % markersize
lw = 2; % linewidth
style1 = '+-';
style2 = 'o-';
style3 = 'd-';
style4 = 'v-';
col1 = [.2 .2 .2];
col2 = [.4 .4 .4];
col3 = [.6 .6 .6];
col4 = [.7 .7 .7];

% col1 = 'green';
% col2 = 'blue';
% col3 = 'black';
% col4 = 'red';

gap = [0.05 0.09]; % vert, horiz gap between adjacent axes
vertmargin = [0.08 0.025]; % bottom, top margins
horizmargin = [0.05 0.02]; % left, right margin

fontsize = 12;
figwidth = 12;
figheight = 5;

%% RegularizedW algorithm

load 'regularizedWNystromExtensionDatasetDenseU1'
densemeanerr = meanerr;
load 'regularizedWNystromExtensionDatasetSparseU1'
sparsemeanerr = meanerr;

myplot = @semilogy;

figure();
subtightplot(1,2,1,gap,vertmargin,horizmargin);
myplot(llist, densemeanerr(1,:)/opterr, style1, 'Color', col1, 'MarkerFaceColor', col1, 'MarkerSize', ms, 'LineWidth', lw);
hold on;
myplot(llist, densemeanerr(2,:)/opterr, style2, 'Color', col2, 'MarkerFaceColor', col2, 'MarkerSize', ms, 'LineWidth', lw);
myplot(llist, densemeanerr(3,:)/opterr, style3, 'Color', col3, 'MarkerFaceColor', col3, 'MarkerSize', ms, 'LineWidth', lw);
myplot(llist, densemeanerr(4,:)/opterr, style4, 'Color', col4, 'MarkerFaceColor', col4, 'MarkerSize', ms, 'LineWidth', lw);
xlab = xlabel('$\ell$ (column samples)');
set(xlab, 'Interpreter', 'Latex')
yl = ylabel('$\|\mathbf{A} - \mathbf{C}\mathbf{W}^\dagger\mathbf{C}^T\|_2/\lambda_{11}(\mathbf{A})$');
set(yl, 'Interpreter', 'Latex');
% remember that mu = n/k * tau
leg = legend('$$\mu= 1$$', '$\mu = 17$', '$\mu = 34$', '$\mu = 50$');
set(leg, 'Interpreter', 'Latex')
hold off;

subtightplot(1,2,2,gap,vertmargin,horizmargin);
myplot(llist, sparsemeanerr(1,:)/opterr, style1, 'Color', col1, 'MarkerFaceColor', col1, 'MarkerSize', ms, 'LineWidth', lw);
hold on;
myplot(llist, sparsemeanerr(2,:)/opterr, style2, 'Color', col2, 'MarkerFaceColor', col2, 'MarkerSize', ms, 'LineWidth', lw);
myplot(llist, sparsemeanerr(3,:)/opterr, style3, 'Color', col3, 'MarkerFaceColor', col3, 'MarkerSize', ms, 'LineWidth', lw);
myplot(llist, sparsemeanerr(4,:)/opterr, style4, 'Color', col4, 'MarkerFaceColor', col4, 'MarkerSize', ms, 'LineWidth', lw);
xlab = xlabel('$\ell$ (column samples)');
set(xlab, 'Interpreter', 'Latex')
yl = ylabel('$\|\mathbf{A} - \mathbf{C}\mathbf{W}^\dagger\mathbf{C}^T\|_2/\lambda_{11}(\mathbf{A})$');
set(yl, 'Interpreter', 'Latex');
% remember that mu = n/k * tau
leg = legend('$$\mu= 1$$', '$\mu = 17$', '$\mu = 34$', '$\mu = 50$');
set(leg, 'Interpreter', 'Latex')
hold off;

printcf('RegularizedWDenseVsSparse.pdf', fontsize, figwidth, figheight);

%% RegularizedA algorithm

load 'regularizedANystromExtensionDatasetDenseU1'
densemeanerr = meanerr;
load 'regularizedANystromExtensionDatasetSparseU1'
sparsemeanerr = meanerr;

myplot = @semilogy;

% for talk
figure();
myplot(llist, densemeanerr(1,:)/opterr, style1, 'Color', col1, 'MarkerFaceColor', col1, 'MarkerSize', ms, 'LineWidth', lw);
hold on;
myplot(llist, densemeanerr(2,:)/opterr, style2, 'Color', col2, 'MarkerFaceColor', col2, 'MarkerSize', ms, 'LineWidth', lw);
myplot(llist, densemeanerr(3,:)/opterr, style3, 'Color', col3, 'MarkerFaceColor', col3, 'MarkerSize', ms, 'LineWidth', lw);
myplot(llist, densemeanerr(4,:)/opterr, style4, 'Color', col4, 'MarkerFaceColor', col4, 'MarkerSize', ms, 'LineWidth', lw);
xlab = xlabel('$\ell$ (column samples)');
set(xlab, 'Interpreter', 'Latex')
yl = ylabel('$\|\mathbf{A} - \mathbf{C}\mathbf{W}^\dagger\mathbf{C}^T\|_2/\lambda_{11}(\mathbf{A})$');
set(yl, 'Interpreter', 'Latex');
% remember that mu = n/k * tau
leg = legend('$$\mu= 1$$', '$\mu = 17$', '$\mu = 34$', '$\mu = 50$');
set(leg, 'Interpreter', 'Latex')
hold off;

printcf('nystromcoherencedependence.pdf', 16, 8, 6);

figure();
subtightplot(1,2,1,gap,vertmargin,horizmargin);
myplot(llist, densemeanerr(1,:)/opterr, style1, 'Color', col1, 'MarkerFaceColor', col1, 'MarkerSize', ms, 'LineWidth', lw);
hold on;
myplot(llist, densemeanerr(2,:)/opterr, style2, 'Color', col2, 'MarkerFaceColor', col2, 'MarkerSize', ms, 'LineWidth', lw);
myplot(llist, densemeanerr(3,:)/opterr, style3, 'Color', col3, 'MarkerFaceColor', col3, 'MarkerSize', ms, 'LineWidth', lw);
myplot(llist, densemeanerr(4,:)/opterr, style4, 'Color', col4, 'MarkerFaceColor', col4, 'MarkerSize', ms, 'LineWidth', lw);
xlab = xlabel('$\ell$ (column samples)');
set(xlab, 'Interpreter', 'Latex')
yl = ylabel('$\|\mathbf{A} - \mathbf{C}\mathbf{W}^\dagger\mathbf{C}^T\|_2/\lambda_{11}(\mathbf{A})$');
set(yl, 'Interpreter', 'Latex');
% remember that mu = n/k * tau
leg = legend('$$\mu= 1$$', '$\mu = 17$', '$\mu = 34$', '$\mu = 50$');
set(leg, 'Interpreter', 'Latex')
hold off;

subtightplot(1,2,2,gap,vertmargin,horizmargin);
myplot(llist, sparsemeanerr(1,:)/opterr, style1, 'Color', col1, 'MarkerFaceColor', col1, 'MarkerSize', ms, 'LineWidth', lw);
hold on;
myplot(llist, sparsemeanerr(2,:)/opterr, style2, 'Color', col2, 'MarkerFaceColor', col2, 'MarkerSize', ms, 'LineWidth', lw);
myplot(llist, sparsemeanerr(3,:)/opterr, style3, 'Color', col3, 'MarkerFaceColor', col3, 'MarkerSize', ms, 'LineWidth', lw);
myplot(llist, sparsemeanerr(4,:)/opterr, style4, 'Color', col4, 'MarkerFaceColor', col4, 'MarkerSize', ms, 'LineWidth', lw);
xlab = xlabel('$\ell$ (column samples)');
set(xlab, 'Interpreter', 'Latex')
yl = ylabel('$\|\mathbf{A} - \mathbf{C}\mathbf{W}^\dagger\mathbf{C}^T\|_2/\lambda_{11}(\mathbf{A})$');
set(yl, 'Interpreter', 'Latex');
% remember that mu = n/k * tau
leg = legend('$$\mu= 1$$', '$\mu = 17$', '$\mu = 34$', '$\mu = 50$');
set(leg, 'Interpreter', 'Latex')
hold off;

printcf('RegularizedADenseVsSparse.pdf', fontsize, figwidth, figheight);

%% Chiu--Demanet truncatedW algorithm

load 'truncatedWNystromExtensionDatasetDenseU1'
densemeanerr = meanerr;
load 'truncatedWNystromExtensionDatasetSparseU1'
sparsemeanerr = meanerr;

myplot = @semilogy;

figure();
subtightplot(1,2,1,gap,vertmargin,horizmargin+[.005,0]);
myplot(llist, densemeanerr(1,:)/opterr, style1, 'Color', col1, 'MarkerFaceColor', col1, 'MarkerSize', ms, 'LineWidth', lw);
hold on;
myplot(llist, densemeanerr(2,:)/opterr, style2, 'Color', col2, 'MarkerFaceColor', col2, 'MarkerSize', ms, 'LineWidth', lw);
myplot(llist, densemeanerr(3,:)/opterr, style3, 'Color', col3, 'MarkerFaceColor', col3, 'MarkerSize', ms, 'LineWidth', lw);
myplot(llist, densemeanerr(4,:)/opterr, style4, 'Color', col4, 'MarkerFaceColor', col4, 'MarkerSize', ms, 'LineWidth', lw);
xlab = xlabel('$\ell$ (column samples)');
set(xlab, 'Interpreter', 'Latex')
yl = ylabel('$\|\mathbf{A} - \mathbf{C}\mathbf{W}^\dagger\mathbf{C}^T\|_2/\lambda_{11}(\mathbf{A})$');
set(yl, 'Interpreter', 'Latex');
% remember that mu = n/k * tau
leg = legend('$$\mu= 1$$', '$\mu = 17$', '$\mu = 34$', '$\mu = 50$');
set(leg, 'Interpreter', 'Latex')
hold off;

subtightplot(1,2,2,gap,vertmargin,horizmargin);
myplot(llist, sparsemeanerr(1,:)/opterr, style1, 'Color', col1, 'MarkerFaceColor', col1, 'MarkerSize', ms, 'LineWidth', lw);
hold on;
myplot(llist, sparsemeanerr(2,:)/opterr, style2, 'Color', col2, 'MarkerFaceColor', col2, 'MarkerSize', ms, 'LineWidth', lw);
myplot(llist, sparsemeanerr(3,:)/opterr, style3, 'Color', col3, 'MarkerFaceColor', col3, 'MarkerSize', ms, 'LineWidth', lw);
myplot(llist, sparsemeanerr(4,:)/opterr, style4, 'Color', col4, 'MarkerFaceColor', col4, 'MarkerSize', ms, 'LineWidth', lw);
xlab = xlabel('$\ell$ (column samples)');
set(xlab, 'Interpreter', 'Latex')
yl = ylabel('$\|\mathbf{A} - \mathbf{C}\mathbf{W}^\dagger\mathbf{C}^T\|_2/\lambda_{11}(\mathbf{A})$');
set(yl, 'Interpreter', 'Latex');
% remember that mu = n/k * tau
leg = legend('$$\mu= 1$$', '$\mu = 17$', '$\mu = 34$', '$\mu = 50$');
set(leg, 'Interpreter', 'Latex')
hold off;

printcf('TruncatedWDenseVsSparse.pdf', fontsize, figwidth, figheight);

%% remove whitespace from produced plots

% run truncateplots.pl to crop off whitespace of ALL pdfs in the directory
[status, result] = system('./truncateplots.pl');
clc
disp(result)