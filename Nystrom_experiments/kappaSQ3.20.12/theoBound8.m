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
function [ bound ] = theoBound8( m,n,k,mu, delta )
%NOT YET:This code has been checked on Mar 4, 2012 and is consistent with bound8

if delta<=0 || delta>=1
    error('invalid delta')
end

gamma=k/m;


OptimOptions = optimset('Display','off', 'TolFun',delta*10^-30);
epsilon = fminsearch(@(x) (geTodelta(m,n,gamma,mu,x)-delta)^2,0, OptimOptions); %we find the min value of y
%    epsilon = fminbnd(@(x) (geTodelta(m,n,gamma,mu,x)-delta)^2,0,maxEpsilon); %we find the min value of y


% %These are just some plotting commands I used while writing the function 
% x=linspace(0, max(epsilon+1,1),500);y=x;
% for i=1:length(x);
%     y(i)=(geTodelta(m,n,gamma,mu,x(i))-delta)^2;
% end
% plot(x,y,[epsilon epsilon],[0 1])
% abs(delta-geTodelta(m,n,gamma,mu,epsilon))
% epsilon
% xlabel('epsilon')
% ylabel('delta')
% title('(geTodelta(m,n,gamma,mu,epsilon)-delta)^2')
% pause
% % % abs(delta-geTodelta(m,n,gamma,mu,epsilon))
% % % error('end')
if epsilon<=0 || epsilon>=1%make sure epsilon is of an acceptable value
    bound = sqrt(-1);
else
%check our answer, make sure we obtained delta within 1% of the desired delta
% abs(delta-geTodelta(m,n,gamma,mu,epsilon))/delta
if abs(delta-geTodelta(m,n,gamma,mu,epsilon))/delta<10^-12 %double check answer,  I have an error tolerance for delta which I hardcoded to appear on the plot.
    bound=sqrt((1+epsilon)/(1-epsilon));
else
    bound=sqrt(-1);
end
end

end
function delta = geTodelta(~,n,gamma,mu,epsilon)%THIS is the delta function from bound8
delta = n*((exp(-epsilon)*(1-epsilon)^(epsilon-1))^(gamma/mu) + (exp(epsilon)*(1+epsilon)^(-epsilon-1))^(gamma/mu));
end

