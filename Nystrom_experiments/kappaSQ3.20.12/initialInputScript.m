try
    load('.kARiPlotterTemp.mat');
    disp('Would you like to reuse your m,n,mu,c from last time?')
    disp('1 - Yes')
    disp('2 - No')
    reuse = condInput('Input choice (default=1): ','reuse','reuse==1 || reuse==2','','1');
    if isequal(reuse,1)
        disp(['m= ',num2str(m)])
        disp(['n= ',num2str(n)])
        if isequal(xAxisType,1)
            disp(['c= ',cVals])
            disp('x-axis is mu')
            disp(['mu= ',muVals])
        elseif isequal(xAxisType,2)
            disp(['mu= ',muVals])
            disp('x-axis is c')
            disp(['c= ',cVals])
        elseif isequal(xAxisType,3)
            disp(['c= ',cVals])
            disp(['mu= ',muVals])
            disp('x-axis is c')
        end
        
    else
        error('collect data')
    end
catch
    disp('Please input integers m and n such that 0<n<m.')
    [m,mStr]=condInput('Input m (default=10^4): ','m','m>0 && isequal(m,round(m))','','10000');
    [n,nStr]=condInput('Input n (default=2): ','n','n>0 && isequal(n,round(n)) && n<m',['m=',mStr,';'],'2');
    disp('1 - mu')
    disp('2 - c')
    disp('3 - c, AND mu is a function of c')
    [xAxisType,~]=condInput('Input x axis (default=2): ','xAxisType','xAxisType==1 || xAxisType==2 || xAxisType==3','','2');
    disp(' ');
    disp('For the rest of the inputs you may enter matlab syntax');
    disp('and use the variables m and n.');disp(' ');
    if isequal(xAxisType,1)
        xlbstr='mu';
        disp('remember, n/m <=mu<=1')
        [muRange,muVals] = condInput('Input mu (default=linspace(n/m,1,5)): ','mu','min(mu)>=n/m && max(mu)<=1 && length(mu)>1',['m=',mStr,';n=',nStr,';'],'linspace(n/m,1,5)');
        xAxis=muRange;
        [cRange,cVals] = condInput('Input c (default=n+round((m-n)/2)): ','c','min(c)>=n && max(c)<=m',['m=',mStr,';n=',nStr,';'],'n+round((m-n)/2)');
    elseif isequal(xAxisType,2)
        xlbstr='c';
        disp('c values must be integers greater or equal to n')
        [cRange,cVals] = condInput('Input c (default=n:n:m): ','c','min(c)>=n && max(c)<=m && length(c)>1',['m=',mStr,';n=',nStr,';'],'n:n:m');
        xAxis=cRange;
        disp('remember, n/m <=mu<=1')
        [muRange,muVals] = condInput('Input mu (default=n/m): ','mu','mu>=n/m && mu<=1',['m=',mStr,';n=',nStr,';'],'n/m');
    elseif isequal(xAxisType,3)
        disp('c values must be integers greater or equal to n')
        [cRange,cVals] = condInput('Input c (default=n:n:m): ','c','min(c)>=n && max(c)<=m && length(c)>1',['m=',mStr,';n=',nStr,';'],'n:n:m');
        disp('remember, n/m <=mu<=1')
        disp('remember, both c and mu are stored as vectors and must be treated as such')
        [muRange,muVals] = condInput('Input mu (default=c/m): ','mu','min(mu)>=n/m && max(mu)<=1 && isequal(size(mu),size(c))',['m=',mStr,';n=',nStr,';c=',cVals,';'],'c/m');
        xlbstr=['c, \mu = ',muVals];
    end
    save('.kARiPlotterTemp.mat')
end