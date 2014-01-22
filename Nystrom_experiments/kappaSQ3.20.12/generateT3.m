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
function [ A, muout ] = generateT3( m,n,coh )

errchk(m,n,coh);
[mt,mt2,nt,nt2]=fun(coh,n);
if mt+mt2>m
    error('I do not know how to compute a matrix for these m and n')
end
if mt==0
    T2=stacksOfDiagonalMatrices(mt2,nt2,coh);
    A=T2;
elseif mt2==0
    T=matricesWithHadamardStructure(mt,nt,coh);
    A=T;
else
    T2=stacksOfDiagonalMatrices(mt2,nt2,coh);
    T=matricesWithHadamardStructure(mt,nt,coh);
    A=blkdiag(T,T2);
end
A=[A;zeros(m-size(A,1),n)];
muout=max(sum(A.^2,2));

end

function [mtout,mt2out,ntout,nt2out] = fun(coh,n)
mmin=inf;
ntout=-1;
nt2out=-1;
for nt2=fliplr([0,2:n])
    nt=n-nt2;
    
    if ceil(1/coh)*nt2<2*nt2;
%     disp('++')
%         m=2*nt2 + 2^(ceil(log2((n-nt2)/coh)))
    mt = 2^(ceil(log2((n-nt2)/coh)));
    mt2=2*nt2;
    
    else
%         disp('--')
%     m=ceil(1/coh)*nt2 + 2^(ceil(log2((n-nt2)/coh)))
    mt = 2^(ceil(log2((n-nt2)/coh)));
    mt2=ceil(1/coh)*nt2;
    end
    m=mt+mt2;
    
    if m<mmin
%         ceil(1/coh)*nt2,2*n
        mmin=m;
        ntout=nt;
        nt2out=nt2;
        mtout=mt;
        mt2out=mt2;
    end
end
end
function [] = errchk(m,n,coh)

if n>m
    error('n<=m')
end
if n==1
    error('m>1')
end
if n/m>coh || 1-coh<-10^-12
    error('n/m<coh1<1')
end

end