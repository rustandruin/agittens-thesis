eval([plotstr2(1:end-1),',''LineWidth'',2);']);%PLOT DATA
if isequal(xAxisType,1)
    axis([min(muRange) max(muRange) 0 1])
elseif isequal(xAxisType,2)
    axis([min(cRange) max(cRange) 0 1])
elseif isequal(xAxisType,3)
    axis([min(cRange) max(cRange) 0 1])
end
title(titlestr2(1:end-2))
ylabel('%')
xlabel(xlbstr)