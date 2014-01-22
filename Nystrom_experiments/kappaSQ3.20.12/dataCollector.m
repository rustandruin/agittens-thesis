%License
%     Copyright (C) 2012  Thomas Wentworth
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%------------- BEGIN CODE --------------
function [ xax, out, xfailureP, failureP, cohErr ] = dataCollector( m, n, runs, kRange, muRange, xAxisType , mtxM, sampleM )
%DATACOLLECTOR compiles all of the experiments to be passed back to the
%main function
pdone=0;
cohErr=0;
global wtBr;
skp=0;%this is used to delete the data points where a condn=inf
if xAxisType==1
    out=zeros(size(muRange,2)*runs,1);%run data will be stored here
    xax=out;%xaxis data will be stored here
    failureP=zeros(size(muRange,2),1);
    xfailureP=muRange;
    for i=1:length(muRange)
        [A,cohA]=mtxM(m,n,muRange(i));
    for j=1:runs
        out(j+(i-1)*runs-skp)=runSampling(m,n,kRange, A, sampleM);
        xax(j+(i-1)*runs-skp)=cohA;%muRange(i); INPUT TRUE COHERENC
        if isequal(out(j+(i-1)*runs-skp),inf)
            skp=skp+1;
            failureP(i)=failureP(i)+1;
        end
        if pdone<((i-1)*runs+j)/(length(muRange)*runs)
            pdone=pdone+.01;
            waitbar(((i-1)*runs+j)/(length(muRange)*runs),wtBr)
        end
    end
    end
elseif xAxisType==2
    out=zeros(round(size(kRange,2)*runs),1);%run data will be stored here
    xax=out;%xaxis data will be stored here
    failureP=zeros(size(kRange,2),1);
    xfailureP=kRange;
    cohErr=0;
        %When k is the x-axis we want to have all of the matrices have the same
    %coherenc  This is not always possible, so we will print the error
    for i=1:length(kRange)
        [A,cohA]=mtxM(m,n,muRange);
        cohErr=max(cohErr,abs(cohA-muRange)/muRange*100);
    for j=1:runs
        out(j+(i-1)*runs-skp)=runSampling(m,n,kRange(i), A, sampleM);
        xax(j+(i-1)*runs-skp)=kRange(i);
        if isequal(out(j+(i-1)*runs-skp),inf)
            skp=skp+1;
            failureP(i)=failureP(i)+1;
        end
        if pdone<((i-1)*runs+j)/(length(kRange)*runs)
            pdone=pdone+.01;
            waitbar(((i-1)*runs+j)/(length(kRange)*runs),wtBr)
        end
    end
    end
elseif xAxisType==3
    out=zeros(round(size(kRange,2)*runs),1);%run data will be stored here
    xax=out;%xaxis data will be stored here
    failureP=zeros(size(kRange,2),1);
    xfailureP=kRange;
    cohErr=0;
        %When k is the x-axis we want to have all of the matrices have the same
    %coherenc  This is not always possible, so we will print the error
    for i=1:length(kRange)
        [A,cohA]=mtxM(m,n,muRange(i));
        cohErr=max(cohErr,abs(cohA-muRange(i))/muRange(i)*100);
    for j=1:runs
        out(j+(i-1)*runs-skp)=runSampling(m,n,kRange(i), A, sampleM);
        xax(j+(i-1)*runs-skp)=kRange(i);
        if isequal(out(j+(i-1)*runs-skp),inf)
            skp=skp+1;
            failureP(i)=failureP(i)+1;
        end
        if pdone<((i-1)*runs+j)/(length(kRange)*runs)
            pdone=pdone+.01;
            waitbar(((i-1)*runs+j)/(length(kRange)*runs),wtBr)
        end
    end
    end
else
    error('invalid xAxisType')
end
failureP=failureP/runs;
xax=xax(1:end-skp);
out=out(1:end-skp);


end

