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
function [ outM, cohA, b ] = specificCohIntRot( A,mu, errTol, b)
%%%%%%
%This function will use a given input matrix A to find a matrix with
%coherence=mu.
%  The algorithm is a modification of minCohIntRot and
%works by first computing the QR decomposition A=UR.  It then rotates 
%two rows of U (i and j) with the largest and the smallest
%row norms.  It will rotate row i and j such that the row norm of row i is
%mu,if it is possible, or as large as possible if mu is not possible.
%It will stop when the correct coherence is obtained (or within errTol) OR
%if m-1 of the rows have been solved for exactly.  If m-1 of the rows have
%a row norm of mu then the last row must also have a row norm of mu (if mu
%=n/m) or less than mu (if mu>n/m).  If mu=n/m then the code will run the
%same was as minCohIntRot.m
%
%Code written by Thomas Wentworth, Feb 15, 2012.
%thomas_wentworth@ncsu.edu
%%%%%%

slvRun=0;
[m,n]=size(A);
if mu>1 || mu<n/m
    error('invalid coherence')
end
if nargin <3
    errTol=0;
end
[A,R]=qr(A,0);  %We begin by orthogonalizing and normalizing the rows of A
rowNorms=sum(A.^2,2); %We compute the row norms
rowNormsUsed=zeros(size(rowNorms));
while slvRun<m%for ii=1:(m-1)% run exact solver at most m-1 times
    [~,i]=max(rowNorms);i=i(1); %find the row with the largest row norm
    if mu>=max(rowNorms)
        [~,j]=min(rowNorms+rowNormsUsed);j=j(1); %find the row with the smallest row norm, but don'y use the same one as last time
        if i==j %we tried them all once already, lets try again
            rowNormsUsed=zeros(size(rowNorms));
            [~,j]=min(rowNorms);j=j(1);
        end
    else
        [~,j]=min(rowNorms);j=j(1); %find the row with the smallest row norm
    end
    [ct,st]=maxacsSol(A,i,j);%find the max possible row norm.
    At=[ct st;-st ct]*A([i,j],:);
    rowNormsAt=sum(At.^2,2);
    if mu<=max(rowNormsAt) %if it is solvable now, do it
        [c,s]=DiffcsSol(A,i,j,mu);
        A([i,j],:)=[c s;-s c]*A([i,j],:);%Rotate rows i and j by thetab radians
        rowNormsUsed=zeros(size(rowNorms));%now we can use all rows again
        slvRun=slvRun+1;
    else
        c=ct;
        s=st;%otherwise make the row a s big as possible
        rowNormsUsed(j)=inf;%don't use the small row again, you already sucked the life out of it.
        if rowNormsAt(1)<rowNormsAt(2)
            ttmmpp=i;%we want the top row to have the largest rownorm
            i=j;
            j=ttmmpp;
        end
        A([i,j],:)=At;%Rotate rows i and j by thetab radians
        %     sum(A([i,j],:).^2,2)
    end
    if nargin >3
        b([i,j],:)=[c s;-s c]*b([i,j],:);%Rotate rows i and j by thetab radians
    end
    rowNorms(j)=sum(A(j,:).^2,2); %We compute the new row norm
    rowNorms(i)=sum(A(i,:).^2,2); %We compute the new row norm
    if slvRun>=m-1 || norm(max(rowNorms)-mu)<errTol
        %if we finish early, lets stop
        break;
    end
end
outM=A*R;
cohA=max(rowNorms);
end

function [c,s]=DiffcsSol(A,i,j,mu)%I previously verified that this is solvable
%here I use the method from dhillon and tropp to solve
a=A(i,:)';aa=a'*a;
b=A(j,:)';bb=b'*b;ab=a'*b;
nom=mu;
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
function [c,s]=maxacsSol(A,i,j)%I used the derivative to find the max possible row norm
a=A(i,:)';aa=a'*a;
b=A(j,:)';bb=b'*b;ab=a'*b;
if norm(ab)>0
    tmp1=-(aa-bb)/(2*ab);
    tmp2=sqrt((aa-bb)^2 + 4*ab^2)/(2*ab);
    if sign(tmp1)==sign(tmp2)
        tp=tmp1+tmp2;
    else
        tp=tmp1-tmp2;
    end
else
    tp=0;
    %note, if aa=bb and ab=0 then any t will work
end
c=1/sqrt(1+tp^2);
s=c*tp;
end