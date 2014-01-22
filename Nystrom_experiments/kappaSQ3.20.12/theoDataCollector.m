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

function [ xax, out ] = theoDataCollector( m, n, kRange, muRange, delta, xAxisType , theoFn )
%THEODATACOLLECTOR compiles all of the theoritical bounds to be passed back to the
%main function
skp=0;%this is used to delete the data points where a theoretical bound could not be established.
pdone=0;
global wtBr;
if xAxisType==1
    out=zeros(size(muRange));%run data will be stored here
    xax=muRange;%xaxis data will be stored here
    for i=1:length(muRange)
        
        out(i-skp)=theoFn(m,n,kRange, muRange(i), delta);
        xax(i-skp)=muRange(i);
        if imag(out(i-skp))>0
            skp=skp+1;
        end
        if pdone<i/length(muRange)
            pdone=pdone+.01;
            waitbar(i/length(muRange),wtBr)
        end
    end
elseif xAxisType==2
    out=zeros(size(kRange));%run data will be stored here
    xax=kRange;%xaxis data will be stored here
    for i=1:length(kRange)
        out(i-skp)=theoFn(m,n,kRange(i), muRange, delta);
        xax(i-skp)=kRange(i);
        if imag(out(i-skp))>0
            skp=skp+1;
        end
        if pdone<i/length(kRange)
            pdone=pdone+.01;
            waitbar(i/length(kRange),wtBr)
        end
    end
elseif xAxisType==3
    out=zeros(size(kRange));%run data will be stored here
    xax=kRange;%xaxis data will be stored here
    for i=1:length(kRange)
        out(i-skp)=theoFn(m,n,kRange(i), muRange(i), delta);
        xax(i-skp)=kRange(i);
        if imag(out(i-skp))>0
            skp=skp+1;
        end
        if pdone<i/length(kRange)
            pdone=pdone+.01;
            waitbar(i/length(kRange),wtBr)
        end
    end
else
    error('invalid xAxisType')
end
xax=xax(1:end-skp);
out=out(1:end-skp);

end

