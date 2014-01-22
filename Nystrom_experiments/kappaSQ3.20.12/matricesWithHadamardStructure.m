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
function [ T,muout ] = matricesWithHadamardStructure( N, n, coh )
%This code will create a matrix as described in section 5.3, "matrices with
%hadamard structure."
%We perform these operations efficiently in that we only store the
%necessary entries.

%if zeros is 0 then we will try to have as few zero rows as possible, if
%zeros is 1 then we will try to have lots of zero rows

if abs(log2(N)-round(log2(N)))> 10^-14
    disp('hi')
    error('N must be a power of 2, it is possible that there are other sizes (like the hadamard mtx), however I dont know how to make them yet')
end
if n>N-1
    error('n<=N')
end
if N==1
    error('N>1')
end
p=round(log2(N));

alpha=sqrt((coh-(n-1)/(N-1))/(1-(n-1)/(N-1)));
beta=sqrt((1-alpha^2)/(N-1));
T_a2=alpha;
T_b2=beta;
T_a=T_a2;
T_b=T_b2;
for j=2:p
    if size(T_a,2)<n
    T_a=[T_a,-T_b;T_b,T_a];
    T_b=[-T_b,T_b;T_b,T_b];
    else
        T_a=[T_a;T_b];%only store the entries we need
        T_b=[-T_b;T_b];
    end
    if size(T_a,2)>n 
        T_a=T_a(:,1:n);
        T_b=T_b(:,1:n);
    end
end
if size(T_a,2)<n
T=[T_a,-T_b;T_b,T_a];
else
    T=[T_a;T_b];
end
if size(T_a,2)>n 
T=T(:,1:n);
else
end
muout=coherence(T);
end

