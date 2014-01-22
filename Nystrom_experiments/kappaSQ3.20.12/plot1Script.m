if loadPrev==2
if cohErr>0
    %When c is the x-axis we want to have all of the matrices have the same
    %coherenc  This is not always possible, so we will print the error
    cohErrStatement=['(+/- < ',num2str(cohErr),'%)'];
else
    cohErrStatement='';
end
if delta>0 && ~isequal(deltaStr(1),',')%make sure we don't do it twice
    deltaStr=[', \delta=',num2str(delta),deltaStr];
end
end
eval([plotstr(1:end-1),',''LineWidth'',2);']);%PLOT DATA
eval([lgndstr(1:end-1),');']);%ADD LEGEND
xlabel(xlbstr)%ADD XLABEL
ylabel(ylbstr)%ADD YLABEL
if isequal(xAxisType,1)
    title([titlestr(1:end-2), '\newline m=',num2str(m),', n=',num2str(n),mtxgstr,', c=',num2str(cVals),titleAdd,deltaStr])%ADD TITLE for xAxisType 1
    axis([min(muRange) max(muRange) yPlotmin max(11,yPlotmax)])
elseif isequal(xAxisType,2)
    title([titlestr(1:end-2), '\newline m=',num2str(m),', n=',num2str(n),mtxgstr,', \mu=',num2str(muVals),cohErrStatement,deltaStr])%ADD TITLE for xAxisType 2
    axis([min(cRange) max(cRange) yPlotmin max(11,yPlotmax)])
elseif isequal(xAxisType,3)
    title([titlestr(1:end-2), '\newline m=',num2str(m),', n=',num2str(n),mtxgstr,', c=',num2str(cVals),', \mu=',num2str(muVals),cohErrStatement,deltaStr])%ADD TITLE for xAxisType 3
    axis([min(cRange) max(cRange) yPlotmin max(11,yPlotmax)])
end