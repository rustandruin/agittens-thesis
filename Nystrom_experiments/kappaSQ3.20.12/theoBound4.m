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
function [ bound ] = theoBound4( m,n,k,mu, delta )
%THEOKAPPAFROMFASTMONTECARLOALGORITHMSFORMATRICES

%This code has been checked on Mar 4, 2012 and is consistent with bound4

if delta<=0 || delta>=1
    error('invalid delta')
end

epsilon = sqrt((m*n*mu)/k) + m*mu*sqrt(8*log(1/delta)/k);

if epsilon<=0 || epsilon>=1
    bound = sqrt(-1);
else

bound=sqrt((1+epsilon)/(1-epsilon));
end

end

