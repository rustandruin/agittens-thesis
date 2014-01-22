function printcf(outputfname, fontsize, width, height)
% printcf(outputfname, fontsize, width, height)
%
% Prints the current figure w/ the given dimensions and all text at the
% specified fontsize to a pdf file named in outputfname
%

% handle subplots appropriately
%set(findobj(get(gca, 'Parent'), 'Type', 'axes'), 'FontSize', fontsize, 'FontWeight', 'Bold')
set(findobj(get(gca, 'Parent'), 'Type', 'axes'), 'FontSize', fontsize, 'FontName', 'Times')
titlehs = get(findobj(get(get(gca, 'Parent'), 'Children'), 'Type','axes'), 'Title');
if length(titlehs) > 1
    for hidx = 1:length(titlehs)
        set(titlehs{hidx}, 'FontSize', fontsize, 'FontName', 'Times');
       % set(titlehs{hidx}, 'FontSize', fontsize, 'FontWeight', 'Bold');
    end
else
    set(titlehs, 'FontSize', fontsize, 'FontName', 'Times');
   % set(titlehs, 'FontSize', fontsize, 'FontWeight', 'Bold');
end
set(findall(gcf,'type','text'),'FontSize',fontsize, 'FontName', 'Times');
set(gcf, 'PaperUnits', 'inches')
set(gcf, 'PaperSize', [width height])
set(gcf, 'PaperPosition', [0 0 width height])
print(gcf, '-dpdf', outputfname);
end