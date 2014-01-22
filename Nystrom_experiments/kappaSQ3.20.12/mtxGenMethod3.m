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
function [ A,mufinal ] = mtxGenMethod3( m,n,mu )
mSmall=min(m,ceil(n/mu));%m-mSmall is maximal number of zero rows
[As,mufinal]=mtxGenMethod2(mSmall,n,mu);
A=[As;zeros(m-mSmall,n)];
end

