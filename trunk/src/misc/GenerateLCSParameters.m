function [sinusCoefficients sinusPhases sinusFrequencies] = ...
    GenerateLCSParameters(samplingPeriod,numberOfFrequencies)

%Generate the parameters of a Linear Combination of Sinusoids (LCS).
%
% USAGE: [sinusCoefficients sinusPhases sinusFrequencies] = GenerateLCSParameters(samplingPeriod,Length)
%
% PARAMETERS:
%   samplingPeriod ~ the sampling period is used not to generate
%     frequencies larger than the Nyquist frequency.
%   Length ~ Desired length of the signal.
%
% RETURN VALUE:
%   sinusCoefficients, sinusPhases, sinusFrequencies ~ Parameters generated
%     for a LCS.
%
%   see also GenerateLCSSignal

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

    % Nyquist frequency
    NyquistFrequency = 1/(2*samplingPeriod);
    MinFrequency = NyquistFrequency/(numberOfFrequencies-1);
    StepFrequency = (NyquistFrequency-MinFrequency)/(numberOfFrequencies-1);
    
    % Generate frequencies, coefficients and samples per each component
    sinusCoefficients = rand(numberOfFrequencies,1);
    sinusPhases = 2*pi*rand(numberOfFrequencies,1);
    sinusFrequencies = MinFrequency:StepFrequency:NyquistFrequency;
    
end
