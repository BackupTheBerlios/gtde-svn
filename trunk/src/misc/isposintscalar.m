function f = isposintscalar(n)

%True if it is a positivive integer scalar
%
% USAGE: f = isposintscalar(n)
%
% PARAMETERS:
%   n ~ the number to test
%
% RETURN VALUE:
%   f ~ boolean
%

% Copyright 2012, Xavier Alameda-Pineda
% INRIA Grenoble Rhône-Alpes
% E-mail: xavi.alameda@gmail.com
% 
% This is part of the gtde program.
% 
% gtde is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

    if numel(n) ~= 1
        error('The input parameter should be one element.');
    end
    f = true;
    if ( (~isreal(n)) || (n<=0) || (round(n) ~= n) )
        f = 0;
    end

end
