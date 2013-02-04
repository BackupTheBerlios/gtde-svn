function polynomialCoefficients = PolynomialInterpolationCoefficients(signal,step)

%Compute the polynomial interpolation coefficients
%
% USAGE: polynomialCoefficients = PolynomialInterpolationCoefficients(signal,step)
%
% PARAMETERS:
%  signal ~ the value of the sampled signal
%  step ~ the sampling step
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
    
    %%% State and solve the linear system for natural splines
    % Vector
    b = signal(1:end-2) - 2*signal(2:end-1) + signal(3:end);
    b = b*6/step^2;
    % Matrix
    A = zeros(length(b));
    A(1:end-1,2:end) = A(1:end-1,2:end) + eye(length(b)-1);
    A(2:end,1:end-1) = A(2:end,1:end-1) + eye(length(b)-1);
    A = A + 4*eye(length(b));
    % Solve system
    M = A\b';
    % Add the two zeros to make natural splines
    M = [0 M' 0];
    
    %%% Allocate and compute polynomial coefficients
    polynomialCoefficients = zeros(length(signal)-1,4);
    % Zero order
    polynomialCoefficients(:,1) = signal(1:end-1);
    % First order
    polynomialCoefficients(:,2) = (signal(2:end)-signal(1:end-1))/step - (M(2:end)+2*M(1:end-1))*step/6;
    % Second order
    polynomialCoefficients(:,3) = M(1:end-1)/2; 
    % Third order
    polynomialCoefficients(:,4) = (M(2:end)-M(1:end-1))/(6*step);
    
%     %%% Linear interpolation
%     polynomialCoefficients = zeros(length(signal)-1,2);
%     polynomialCoefficients(:,1) = signal(1:end-1);
%     polynomialCoefficients(:,2) = (signal(2:end)-signal(1:end-1))/step;
return

