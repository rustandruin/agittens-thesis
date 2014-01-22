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
function [ bound ] = theoBound1(m,~,k,mu, delta  )
%THEOKAPPAFROMFASTERLEASTSQUARESAPPROXIMATION

%This code has been checked on Mar 4, 2012 and is consistent with bound1

if delta<=0 || delta>=1
    error('invalid delta')
end
cc=96*m*mu;
if cc*log(cc/sqrt(delta)) > k
    bound = sqrt(-1);
else
    epsilon = ((m*mu*96)/delta^(1/2))^(1/2)/exp(lambertw(0, k/delta^(1/2))/2);
    if epsilon<=0 || epsilon>=1
        bound = sqrt(-1);
    else
        bound=sqrt((1+epsilon)/(1-epsilon));
    end
end

end