function [index sign] = TDEis(mic1,mic2,NMics)

%Index and sign of a TDE given the two microphones' indices
%
% USAGE: [index sign] = TDEis(mic1,mic2,NMics)
%
% PARAMETERS:
%   mic1, mic2 ~ indices of the first and second microphones
%   NMics ~ number of microphones
%
% RETURN VALUE:
%   index, sign ~ self-explanatory
% 
% DESCRIPTION:
%   Computes the index of the tde that corresponds to mic1 and mic2 on a
%   vectorized representation.
%
% EXAMPLE:
%   In case we have five microphones. The vectorized representation of the
%   tde's is:
%       t = (t12 t13 t14 t15 t23 t24 t25 t34 t35 t45)
%   In that case: 
%       [index,sign] = TDEis(2,3,5) gives index = 5 and sign = 1.
%       [index,sign] = TDEis(4,3,5) gives index = 8 and sign = -1.
%       [index,sign] = TDEis(6,1,5) gives an error.

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

    %%% Check input
    if nargin < 3
        error('Usage: [index sign] = TDEis(mic1,mic2,NMics)');
    end
    if numel(mic1)~=1 || numel(mic2)~=1 || numel(NMics)~=1
        error('mic1, mic2 and Nmics should be numbers.');
    end
    if mic1 < 1 || mic1 > NMics
        error('mic1 should be between 1 and NMics');
    end
    if mic2 < 1 || mic2 > NMics
        error('mic1 should be between 1 and NMics');
    end

    %%% Compute the index
    if mic1 < mic2
        index = (mic1-1)*NMics - sum(1:(mic1-1)) + (mic2-mic1);
        sign = 1;
    else
        index = (mic2-1)*NMics - sum(1:(mic2-1)) + (mic1-mic2);
        sign = -1;
    end
end