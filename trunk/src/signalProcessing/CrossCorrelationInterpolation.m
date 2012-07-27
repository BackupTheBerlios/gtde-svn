function [cc ccd ccdd] = CCInterpolation(PCCC,Delay,SamplingPeriod,ZeroIndex)

%ccInterpolation Cross-correlation function interpolation
%
% USAGE: [cc ccd ccdd] = CCInterpolation(PCCC,Delay,T[,ZI])
%
% PARAMETERS:
%  PCCC ~ cross-correlation signals of the polynomial coefficients.
%  Delay ~ delay is the point(s) in which we want to estimate the signals'
%     cross-correlation function.
%  T ~ the signals' sampling period.
%  ZI ~ the index on the PCCC representing no delay (the default value is ceil(L/2), where L is the
%     size of the coefficients cross-correlation function.)
% 
% RETURN VALUE:
%  cc,ccd,ccdd ~ The value of the cross-correlation function, its first and
%     its second derivative.
% 
% DESCRIPTION:
%     This function estimates the value of the cross-correlation function
%     at Delay from from the polynomial coefficients' cross-correlation 
%     signals PCCC.
%    
% REFERENCES:
%     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
%     estimation-based sound source localisation (gTDESSL). Research Report 
%     RR-7988, INRIA, June 2012.
%
%     see also K1, K2, InterpolateCrossCorrelation

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

    % Input check
    if nargin < 3
        error('Usage: [cc ccd ccdd] = CCInterpolation(PCCC,Delay,SamplingPeriod[,ZeroIndex])');
    end
    if nargin < 4
        ZeroIndex = ceil(size(PCCC,3)/2);
    end

    % Output value
    cc = zeros(size(Delay));
    ccd = zeros(size(Delay));
    ccdd = zeros(size(Delay));
    
    % Index in the cross-correlation
    ccIndex = floor(Delay/SamplingPeriod);
    % Remaining of the division (0 <= tau < SamplingPeriod)
    tau = Delay - ccIndex*SamplingPeriod;
    
    % Interpolation polynomial degrees
    D1 = size(PCCC,1)-1;
    D2 = size(PCCC,2)-1;
    
    for d1 = 0:D1,
        for d2 = 0:D2,
            try
                % Value of the first interpolating function
                [k1 k1d k1dd] = K1(d1,d2,tau,SamplingPeriod);
                % Value of the second interpolation function
                [k2 k2d k2dd] = K2(d1,d2,tau,SamplingPeriod);
                % First cross-correlation value
                R1 = (squeeze(PCCC(d1+1,d2+1,ccIndex+1+ZeroIndex))');
                % Second cross-correlation value
                R2 = (squeeze(PCCC(d1+1,d2+1,ccIndex+ZeroIndex))');
                % Cross-correlation's value
                cc = cc + R1.*k1 + R2.*k2;
                % Cross-correlation's derivative
                ccd = ccd + R1.*k1d + R2.*k2d;
                % Cross-correlation's second derivative
                ccdd = ccdd + R1.*k1dd + R2.*k2dd;
            catch err
%                 disp('Exception');
            end
        end
    end

end