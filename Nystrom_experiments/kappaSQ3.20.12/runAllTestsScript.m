for i=1:length(plotTypes)
    global wtBr;
    wtBr=waitbar(0,['Working on (',plotTypes(i),')']);
    if isequal(plotTypes(i),'1')%MAIN AREA FOR RUNNING MODULES
        [xAxis1,data1, xfailureP1,failureP1, cohErrTmp] = dataCollector(m,n,rns,cRange, muRange, xAxisType , mtxM, @exactlyCSample);
        lgndstr=[lgndstr,'''',plt1LgndStr,''','];%Exactly(c) Experiment
        if coloryn==1
            plotstr=[plotstr,'xAxis1,data1,''bx'','];
        else
            plotstr=[plotstr,'xAxis1,data1,''kx'','];
        end
        titlestr=[titlestr, title1,', '];
        yPlotmax=max(yPlotmax,max(data1));
        
        titlestr2=' Failure %, ';
        if coloryn==1
            plotstr2=[plotstr2, 'xfailureP1, failureP1,''bx-'','];
        else
            plotstr2=[plotstr2, 'xfailureP1, failureP1,''kx-'','];
        end
        cohErr=max(cohErr,cohErrTmp);
    elseif isequal(plotTypes(i),'2')
        [xAxis2,data2, xfailureP2,failureP2,cohErrTmp] = dataCollector(m,n,rns,cRange, muRange, xAxisType , mtxM, @blendSample);
        lgndstr=[lgndstr,'''',plt2LgndStr,''','];%Blendenpik Experiment
        if coloryn==1
            plotstr=[plotstr,'xAxis2,data2,''r+'','];
        else
            plotstr=[plotstr,'xAxis2,data2,''k+'','];
        end
        titlestr=[titlestr, title2,', '];
        yPlotmax=max(yPlotmax,max(data2));
        
        titlestr2=' Failure %, ';
        if coloryn==1
            plotstr2=[plotstr2, 'xfailureP2, failureP2,''r+-'','];
        else
            plotstr2=[plotstr2, 'xfailureP2, failureP2,''k+-'','];
        end
        cohErr=max(cohErr,cohErrTmp);
    elseif isequal(plotTypes(i),'3')
        [xAxis3, data3]=theoDataCollector( m, n, cRange, muRange, delta, xAxisType , @theoBound1);
        if isequal(size(xAxis3,2),0)
            disp(' ');disp('ERROR: None of the test values satisfied the theorem requierements for 3, it will not be plotted.');
            disp('Please note that for plot type 3')
            disp('We must have 0<epsilon<1 and')
            disp('min(n, cc*log(cc/sqrt(delta))) <= c')
            disp('must hold, where')
            disp('cc=96*m*mu')
            disp(' ');
        else
            bndpltcnt=bndpltcnt+1;
            lgndstr=[lgndstr,'''',plt3LgndStr,''','];%Exactly(c) Theo Bnd(3)(bound1)
            if coloryn==1
                plotstr=[plotstr,'xAxis3,data3,''b',plotStylesStore{bndpltcnt},''','];
            else
                plotstr=[plotstr,'xAxis3,data3,''k',plotStylesStore{bndpltcnt},''','];
            end
            titlestr=[titlestr, title3,', '];
            yPlotmax=max(yPlotmax,max(data3));
        end
    elseif isequal(plotTypes(i),'4')
        [xAxis4, data4]=theoDataCollector( m, n, cRange, muRange, delta, xAxisType , @theoBound4);
        if isequal(size(xAxis4,2),0)
            disp(' ');disp('ERROR: None of the test values satisfied the theorem requierements for 4, it will not be plotted.');
            disp('We must have 0<epsilon<1');
            disp(' ');
        else
            bndpltcnt=bndpltcnt+1;
            lgndstr=[lgndstr,'''',plt4LgndStr,''','];%Exactly(c) Theo Bnd(4)(bound4)
            if coloryn==1
                plotstr=[plotstr,'xAxis4,data4,''b',plotStylesStore{bndpltcnt},''','];
            else
                plotstr=[plotstr,'xAxis4,data4,''k',plotStylesStore{bndpltcnt},''','];
            end
            titlestr=[titlestr, title4,', '];
            yPlotmax=max(yPlotmax,max(data4));
        end
    elseif isequal(plotTypes(i),'5')
        [xAxis5, data5]=theoDataCollector( m, n, cRange, muRange, delta, xAxisType , @theoBound6);
        if isequal(size(xAxis5,2),0)
            disp(' ');disp('ERROR: None of the test values satisfied the theorem requierements for 5, it will not be plotted.');
            disp('We must have 0<epsilon<1');
            disp(' ');
        else
            bndpltcnt=bndpltcnt+1;
            lgndstr=[lgndstr,'''',plt5LgndStr,''','];%Blendenpik Theo Bnd(5)(bound6)
            if coloryn==1
                plotstr=[plotstr,'xAxis5,data5,''r',plotStylesStore{bndpltcnt},''','];
            else
                plotstr=[plotstr,'xAxis5,data5,''k',plotStylesStore{bndpltcnt},''','];
            end
            titlestr=[titlestr, title5,', '];
            yPlotmax=max(yPlotmax,max(data5));
        end
    elseif isequal(plotTypes(i),'6')
        [xAxis6, data6]=theoDataCollector( m, n, cRange, muRange, delta, xAxisType , @theoBound3);
        if isequal(size(xAxis6,2),0)
            disp(' ');disp('ERROR: None of the test values satisfied the theorem requierements for 6, it will not be plotted.');
            disp('We must have 0<epsilon<1')
            disp(' ');
        else
            bndpltcnt=bndpltcnt+1;
            lgndstr=[lgndstr,'''',plt6LgndStr,''','];%Exactly(C) Theo Bnd(6)(bound3)
            if coloryn==1
                plotstr=[plotstr,'xAxis6,data6,''b',plotStylesStore{bndpltcnt},''','];
            else
                plotstr=[plotstr,'xAxis6,data6,''k',plotStylesStore{bndpltcnt},''','];
            end
            titlestr=[titlestr, title6,', '];
            yPlotmax=max(yPlotmax,max(data6));
        end
    elseif isequal(plotTypes(i),'7')
        [xAxis7, data7]=theoDataCollector( m, n, cRange, muRange, delta, xAxisType , @theoBound8);
        if isequal(size(xAxis7,2),0)
            disp(' ');disp('ERROR: None of the test values satisfied the theorem requierements for 7, it will not be plotted.');
            disp('We must have 0<epsilon<1')
            disp('we must be able to solve for epsilon')
            disp(' ');
        else
            bndpltcnt=bndpltcnt+1;
            lgndstr=[lgndstr,'''',plt7LgndStr,''','];%Blendenpik Theo Bnd(7)(bound8)
            deltaStr=['+/-',num2str(delta*10^-12),'%'];
            if coloryn==1
                plotstr=[plotstr,'xAxis7,data7,''r',plotStylesStore{bndpltcnt},''','];
            else
                plotstr=[plotstr,'xAxis7,data7,''k',plotStylesStore{bndpltcnt},''','];
            end
            titlestr=[titlestr, title7,', '];
            yPlotmax=max(yPlotmax,max(data7));
        end
%     elseif isequal(plotTypes(i),'8')
%         [xAxis8a,data8a, xfailureP1,failureP1, cohErrTmp] = dataCollector(m,n,rns,cRange, muRange, xAxisType , mtxM, @exactlyCSample);
%         lgndstr=[lgndstr,'''1'','];%Exactly(c) Experiment
%         
%         if coloryn==1
%             plotstr=[plotstr,'xAxis1,data1,''bx'','];
%         else
%             plotstr=[plotstr,'xAxis1,data1,''kx'','];
%         end
%         titlestr=[titlestr, title1,', '];
%         yPlotmax=max(yPlotmax,max(data1));
%         
%         titlestr2=[titlestr2, title1,' Failure %, '];
%         if coloryn==1
%             plotstr2=[plotstr2, 'xfailureP1, failureP1,''bx-'','];
%         else
%             plotstr2=[plotstr2, 'xfailureP1, failureP1,''kx-'','];
%         end
%         cohErr=max(cohErr,cohErrTmp);
%         end
%     else
    end
    delete(wtBr)
end