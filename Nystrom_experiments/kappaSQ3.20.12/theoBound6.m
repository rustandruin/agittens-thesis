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
function [ bound ] = theoBound6( m,n,k,mu, delta )
%THEOKAPPAFROMHIGHDIMENSIONALMATCHEDSUBSPACEDETECTIONWHENDATA
% % % % OLD VERSION
% % % % gamma=k/n;%NEEDS CHANGING TO NEW FORM IF REVERTING TO OLD VERSION
% % % alpha = gamma*n/m;
% % % if alpha >=.5
% % %     omega=(3*((4*log(delta/n)^2)/9 - (8*alpha*log(delta/n)^2)/9 + (4*alpha^2*log(delta/n)^2)/9 - (8*alpha*log(2)^2)/9 - (8*log(2)*log(delta/n))/9 + (4*log(2)^2)/9 + (4*alpha^2*log(2)^2)/9 + (16*alpha*log(2)*log(delta/n))/9 - 8*alpha*m*log(delta/n) - (8*alpha^2*log(2)*log(delta/n))/9 + 8*alpha^2*m*log(delta/n) + 8*alpha*m*log(2) - 8*alpha^2*m*log(2))^(1/2) - 2*log(delta/(2*n)) + 2*alpha*log(delta/(2*n)))/(6*alpha);
% % % % -(2*log(delta/(2*n)) + 3*((4*log(delta/n)^2)/9 - (8*alpha*log(delta/n)^2)/9 + (4*alpha^2*log(delta/n)^2)/9 - (8*alpha*log(2)^2)/9 - (8*log(2)*log(delta/n))/9 + (4*log(2)^2)/9 + (4*alpha^2*log(2)^2)/9 + (16*alpha*log(2)*log(delta/n))/9 - 8*alpha*m*log(delta/n) - (8*alpha^2*log(2)*log(delta/n))/9 + 8*alpha^2*m*log(delta/n) + 8*alpha*m*log(2) - 8*alpha^2*m*log(2))^(1/2) - 2*alpha*log(delta/(2*n)))/(6*alpha)
% % % 
% % % else
% % %    % - log(delta/(2*n))/3 - ((4*alpha*log(delta/n)^2 + 72*m*log(2) + 4*alpha*log(2)^2 - 72*m*log(delta/n) - 8*alpha*log(2)*log(delta/n) + 72*alpha*m*log(delta/n) - 72*alpha*m*log(2))/(9*alpha))^(1/2)/2
% % %    omega=((4*alpha*log(delta/n)^2 + 72*m*log(2) + 4*alpha*log(2)^2 - 72*m*log(delta/n) - 8*alpha*log(2)*log(delta/n) + 72*alpha*m*log(delta/n) - 72*alpha*m*log(2))/(9*alpha))^(1/2)/2 - log(delta/(2*n))/3;
% % %  
% % % end
% % % 
% % % epsilon = omega*mu;
% % % 
% % NEW VERSION



%This code has been checked on Mar 4, 2012 and is consistent with bound6
gamma=k/m;
phi = max(1,(1-gamma)/gamma);
l=2/3*log(2*n/delta);
epsilon = mu/2*(phi*l + sqrt(((1-gamma)/gamma)*12*m*l + (phi*l)^2));
% % 

if 0<epsilon && epsilon<1
    bound=sqrt((1+epsilon)/(1-epsilon));
else
    bound=sqrt(-1);
end

if isequal(k,m)
    bound=1; % if k=m, then all rows are sampled and kappa=1
end

end

