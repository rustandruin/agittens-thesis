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
function [ U,mu] = mtxGenMethod2( m,n,mu )
%MTXGENMETHOD1 Generates a mXn matrix of coherence coh
%aus and zus are unsorted, a and z are sorted
%here we recreate the method

U=zeros(m,n);U(1:n,1:n)=eye(n);
zus=ones(m,1)*(n-mu)/(m-1);zus(end)=mu;
aus=zeros(m,1);aus(1:n)=1;
    [z,~]=sort(zus,'ascend');
for ii=1:m-1

[a,aIndex]=sort(aus,'ascend');

amz=a-z;
amz2=(a(1:end-1)-z(2:end));
% amz2=amz;
% [a z sign(a-z)]
% [a(1:end-1) z(2:end)]
k1=find(abs(amz)<10^-15,1,'first');
k2=find(abs(amz)<10^-15,1,'last');
if isempty(k1)
    k1=0;k2=m;
end

j= find((amz(max(k1+1,1):end))>10^-15,1,'first')+k1;
i= find((amz2(1:(k2-1)))<-10^-15,1,'last');
if isempty(j) || isempty(i)
    break;%sometimes we finish early
end
% i= find((a-z)<0,1,'last')
[c,s] = DiffcsSol(U,aIndex(i),aIndex(j),z(i));
% U([aIndex(i) aIndex(j)],:)
U([aIndex(i) aIndex(j)],:) = [c s;-s c]*U([aIndex(i) aIndex(j)],:);
% U([aIndex(i) aIndex(j)],:)
aus(aIndex(i)) = U(aIndex(i),:)*U(aIndex(i),:)';
aus(aIndex(j)) = U(aIndex(j),:)*U(aIndex(j),:)';

end

end
function [c,s]=DiffcsSol(A,i,j,nom)%I previously verified that this is solvable
%here I use the method from dhillon and tropp to solve
a=A(i,:)';
aa=a'*a;
b=A(j,:)';
bb=b'*b;ab=a'*b;
tmp1=-ab/(bb-nom);
tmp2=(sqrt(ab^2-(aa-nom)*(bb-nom)))/(bb-nom);%we define these variables to avoid sign cancellations.
if tmp1<=0
    tp=tmp1 - tmp2;
else
    tp=tmp1+tmp2;
end
% tp=-(ab + sqrt(ab^2-(aa-nom)*(bb-nom)))/(bb-nom);
c=1/sqrt(1+tp^2);
s=c*tp;
end