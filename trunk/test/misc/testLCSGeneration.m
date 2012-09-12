function testLCSGeneration

%Test the generation of LCS
%
%   see also GenerateLCSSignal, GenerateLCSParameters

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

    % Basic parameters
    SamplingPeriod = 1/48000;
    Times = 0:SamplingPeriod:1;
    N = 1000;

    % LCS parameters and LCS signal
    [c,p,f] = GenerateLCSParameters(SamplingPeriod,N);
    Signal = GenerateLCSSignal(Times,c,p,f);

    % Plot
    plot(Times,Signal);

end