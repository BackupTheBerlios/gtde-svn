function [maxCCValue, maxdCCValue] = MaxCCValues(PCCC,T)

%InterpolatieCrossCorrelation interpolates the cross-correlation function 
%of two discrete signals
%
% USAGE: [maxCCValue, maxdCCValue] = MaxCCValues(PCCC,T)
%
% PARAMETERS:
%   FirstSignal ~ values of the first signal
%   SecondSignal ~ values of the second signal
%   Delays ~ delays in which we want to evaluate the cross-correlation function
%   SamplingPeriod ~ sampling period of the signals
%
% RETURN VALUE:
%   CrossCorrelation ~ values of the cross-correlation function at delays
%   Derivatives ~ values of its derivative at delays
%   Curvature ~ values of its second derivative at delays
% 
% DESCRIPTION:
%   This function interpolates the cross-correlation function at times Delays
%   of the signals FirstSignal and SecondSignal (sampled at SamplingPeriod). 
%   It does that assuming some polynomial interpolation at each time interval.
% 
%   This function will then compute:
%   1) The sequence of coefficients of each polynomial interpolation.
%   2) The cross-correlation of these sequences of coefficients.
%   3) The interpolation of the cross-correlation function at each of the Delays.
%
% REFERENCES:
%     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
%     estimation-based sound source localisation (gTDESSL). Research Report 
%     RR-7988, INRIA, June 2012.
%
%   see also PolynomialInterpolationCoefficients, PolynomialCoefficientsCrossCorrelation and performCCInterpolation

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
    if nargin < 2
        error('Usage: [maxCCValue, maxdCCValue] = MaxCCValues(PCCC,T)');
    end

    % Coefficients
    a = PCCC(4,:);
    b = PCCC(3,:);
    c = PCCC(2,:);
    d = PCCC(1,:);
    
    % Max CC
    % Potential maxima
    sp = (-b + sqrt(b.^2-4*a.*c)) ./ (3*a);
    sn = (-b - sqrt(b.^2-4*a.*c)) ./ (3*a);
    % The right ones
    indp = 0 <= sp & sp < T & isreal(sp);
    indn = 0 <= sn & sn < T & isreal(sn);
    % Maximum
    pot = cat(2,d,...
        a(indp).*sp(indp).^3 + b(indp).*sp(indp).^2 + c(indp).*sp(indp) + d(indp),...
        a(indn).*sn(indn).^3 + b(indn).*sn(indn).^2 + c(indn).*sn(indn) + d(indn));
    maxCCValue = max(abs(pot));
    
    % Max CCd
    % Potential maxima
    s = -b ./ (3*a);
    % The right ones
    ind = 0 <= s & s < T & isreal(s);
    % Potential
    pot = cat(2,c,...
        3*a(ind).*s(ind).^2 + 2*b(ind).*s(ind) + c(ind));
    maxdCCValue = max(abs(pot));
    
%     %%% Max value of the cross-correlation function
%     %%% TODO Change from 1 to better.
%     maxCCValue = 1;
%     
%     % Coefficient cross-correlations
%     a = PCCC(:,:,2:end);
%     b = PCCC(:,:,1:end-1);
% 
%     % Coefficients of the second order derivative of the
%     % cross-correlation
%     sndODC = zeros(6,size(PCCC,3)-1);
%     sndODC(6,:) = b(2,1,:) - b(1,2,:) + (2 * b(3,1,:) * T) - (2 * a(3,1,:) * T) + (a(3,2,:) * T ^ 2) + (a(4,2,:) * T ^ 3) - a(2,1,:) + 0.2e1 / 0.3e1 * b(3,3,:) * (T ^ 3) + (a(2,2,:) * T) + a(1,2,:) + 0.6e1 / 0.5e1 * b(4,4,:) * (T ^ 5) - (3 * a(4,1,:) * T ^ 2) + (2 * b(4,2,:) * T ^ 3) + (b(3,2,:) * T ^ 2) + (b(3,4,:) * T ^ 4) / 0.2e1 + (3 * b(4,1,:) * T ^ 2) + 0.3e1 / 0.2e1 * b(4,3,:) * (T ^ 4);
%     sndODC(5,:) = (2 * a(3,1,:)) - (2 * a(3,2,:) * T) - (3 * a(4,2,:) * T ^ 2) - (2 * b(1,3,:)) + (2 * a(4,3,:) * T ^ 3) - (3 * b(4,2,:) * T ^ 2) - a(2,2,:) + (2 * a(2,3,:) * T) + (6 * a(4,1,:) * T) + (2 * a(3,3,:) * T ^ 2) - (2 * b(4,3,:) * T ^ 3) - (6 * b(4,1,:) * T) - 0.3e1 / 0.2e1 * b(4,4,:) * (T ^ 4) - (2 * b(3,1,:)) + b(2,2,:) + (2 * a(1,3,:));
%     sndODC(4,:) = a(3,2,:) + 3 * a(2,4,:) * T + 3 * a(3,4,:) * T ^ 2 - 2 * a(3,3,:) * T + 3 * a(4,2,:) * T - b(3,2,:) + b(2,3,:) - 3 * a(4,3,:) * T ^ 2 + 3 * a(4,4,:) * T ^ 3 - 3 * b(1,4,:) - 3 * a(4,1,:) + 3 * b(4,1,:) + 3 * a(1,4,:) - a(2,3,:);
%     sndODC(3,:) = -(2 * a(3,4,:) * T) + b(2,4,:) + 0.2e1 / 0.3e1 * a(3,3,:) - a(2,4,:) - a(4,2,:) + (2 * a(4,3,:) * T) + b(4,2,:) - (3 * a(4,4,:) * T ^ 2) - 0.2e1 / 0.3e1 * b(3,3,:);
%     sndODC(2,:) = a(3,4,:) / 0.2e1 - a(4,3,:) / 0.2e1 + 0.3e1 / 0.2e1 * a(4,4,:) * T - b(3,4,:) / 0.2e1 + b(4,3,:) / 0.2e1;
%     sndODC(1,:) = 0.3e1 / 0.10e2 * b(4,4,:) - 0.3e1 / 0.10e2 * a(4,4,:);
%     %%% Max value of the derivative of the cross-correlation function
%     % Loop on the cross-correlation function and keep potential maxs
%     potentialMax = T*(0:size(PCCC,3)-2);
%     for m = 1:size(PCCC,3)-1,
%         % Solving for polynomial roots
%         polRoots = roots(sndODC(:,m));
%         
%         % Check for roots in the interval
%         polRoots = polRoots(isreal(polRoots));
%         polRoots = polRoots( 0 <= polRoots & polRoots < T );
%         
%         % Concatenate
%         potentialMax = cat(1,potentialMax,polRoots + (m-1)*T );
%     end
% 
%     %%% Compute the cross-correlation derivatives
%     [~, dCC] = C2CrossCorrelationInterpolation(PCCC,potentialMax,T,1);
%     
%     %%% Maximum
%     maxdCCValue = max(abs(dCC));
    
end
    