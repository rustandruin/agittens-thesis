numptstobeplotted = 10;
indicestobeplotted = ceil(linspace(1, length(kgrid), numptstobeplotted));
kvalstobeplotted = kgrid(indicestobeplotted);

data1 = srht_forward_spec1(:, indicestobeplotted);
data2 = srht_forward_spec2(:, indicestobeplotted);
data3 = gauss_forward_spec1(:, indicestobeplotted);
data4 = gauss_forward_spec2(:, indicestobeplotted);
data = [data1 data2 data3 data4];

klabels = repmat(kvalstobeplotted, 1, 4);

grouplabels = [repmat({'S1'}, 1, numptstobeplotted), ...
               repmat({'S2'}, 1, numptstobeplotted), ...
               repmat({'G1'}, 1, numptstobeplotted), ...
               repmat({'G2'}, 1, numptstobeplotted)];

boxplot(data, {klabels, grouplabels}, ... % 'plotstyle', 'compact', 
        'color', repmat('rgbk', 1, numptstobeplotted), ... %'factorgap', [5 2], ...
        'labelverbosity', 'minor');
xlabel('k (target rank)')
yh = ylabel('$\|\mathbf{A} - \tilde{\mathbf{A}}\|_2/\|\mathbf{A} - \mathbf{A}_k\|_2$');
set(yh, 'Interpreter', 'latex');

printcf('experimentA-forward-spectral-boxplot.pdf', 13, 14, 7);