function y = GenerateLCSSignal(x,sinusCoefficients,sinusPhases,sinusFrequencies)

%Generate a Linear Combination of Sinusoids (LCS) discrete signal.
%
% USAGE: Signal = GenerateLCSSignal(Times,sinusCoefficients,sinusPhases,sinusFrequencies)
%
% PARAMETERS:
%   Times ~ time instants at which we want to sample the continuous signal.
%   sinusCoefficients ~ coefficients of the linear combination.
%   sinusPhases ~ phase of each of the terms.
%   sinusFrequencies ~ oscillation frequencies of each term.
%
% RETURN VALUE:
%   Signal ~ discrete signal resulting from the sampling.
%
% EXAMPLE:
%   Signal = GenerateLCSSignal(t,[0.2,1],[0,pi/2],[10,20])
%   This will generate a discrete signal which is the sum of two sinusoids
%   sampled at t. The first sinusoid is shifted by 0 radians, at frequency 
%   10 Hz and weitghed by 0.2. The second sinusoid is shifted by pi/2, at
%   frequency 20 Hz and weitghed by 1.
%
%   see also GenerateLCSParameters, SyntheticSignal

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

    % Zero function
    y = zeros(size(x));
    % Add components
    for nC = 1:length(sinusCoefficients),
        y = y + sinusCoefficients(nC)*sin(sinusFrequencies(nC)*x + sinusPhases(nC));
    end
    
end
