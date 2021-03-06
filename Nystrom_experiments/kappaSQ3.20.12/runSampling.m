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
function [ condn ] = runSampling( m,n, kORgamma, A, sampleM )
%Here we test a matrix generated by mtxM times.

As=sampleM(A,m,n,kORgamma);  %<-- this is SQ

if rank(As)<n
    condn=inf;
else
    condn=cond(As);
end

% [~,R]=qr(As,0);  %If you want to compute cond(AR_s^-1), activate this
% if rank(R)<n     %section of code instead.
%     condn=inf;
% else
% condn=cond(A/R);
% end

end

