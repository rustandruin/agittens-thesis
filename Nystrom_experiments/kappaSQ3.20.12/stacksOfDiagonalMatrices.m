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
function [ T, muout ] = stacksOfDiagonalMatrices( m,n,coh)
%This code will create matrices in accordance with sec 5.3, stacks of
%diagonal matrices.

if abs((m/n)-round(m/n))> 10^-14
    error('m/n must be an integer')
end
if (m/n)<2
    error('m/n must be at least 2')
end
if n>m
    error('n<=N')
end
if n==1
    error('N>1')
end

T=makeBaseMtx(m,n,coh);
muout=coherence(T);
end

function T = makeBaseMtx(m,n,coh)
mon=m/n;
% e0=eye(n)/sqrt(mon);   % This is another way to do it
% e1=e0;e1(1,1)=sqrt(coh);
% e2=e0;e2(1,1)=sqrt((1-coh)/(mon-1));


e1=eye(n)*sqrt(coh);
e2=eye(n)*sqrt((1-coh)/(mon-1));
T=e1;
for i=1:mon-1
    T=[T;e2];
end
end