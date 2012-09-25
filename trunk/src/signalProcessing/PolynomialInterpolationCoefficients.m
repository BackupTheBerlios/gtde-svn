function polynomialCoefficients = PolynomialInterpolationCoefficients(signal,times)

%Compute the polynomial interpolation coefficients
%
% USAGE: polynomialCoefficients = PolynomialInterpolationCoefficients(signal,times)
%
% PARAMETERS:
%  signal ~ the value of the sampled signal
%  times ~ the time values when the signal has been sampled (or the
%     sampling period)
%
% RETURN VALUE:
%  polynomialCoefficients ~ the polynomial coefficients of the interpolated
%     signal
% 
% DESCRIPTION:
%     This function computes the coefficients of the linear interpolation
%     of the signal
%
% REFERENCES:
%     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
%     estimation-based sound source localisation (gTDESSL). Research Report 
%     RR-7988, INRIA, June 2012.
%
%   see also InterpolateCrossCorrelation

% Copyright 2012, Xavier Alameda-Pineda
% INRIA Grenoble Rhone-Alpes
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

    %%% Input check
    if length(times) ~= 1 && length(times) ~= length(signal)
        error('times length need to be either the signal length or 1');
    end
    
    %%% Declare output variables
    polynomialCoefficients = zeros(length(signal)-1,2);
    
    %%% Compute the coefficients
    % Zero order
    polynomialCoefficients(:,1) = signal(1:end-1);
    % First order
    if length(times) == 1
        polynomialCoefficients(:,2) = (signal(2:end)-signal(1:end-1))/times;
    else
        polynomialCoefficients(:,2) = (signal(2:end)-signal(1:end-1))./(times(2:end)-times(1:end-1));
    end
    

end