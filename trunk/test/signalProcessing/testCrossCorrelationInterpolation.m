function testCrossCorrelationInterpolation

%test the computation of the polynomial coefficients

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

    clc
    % Plotting period
    tauinc = 0.001;
    % Samplingperiod
    T = 0.01;
    % Length
    Length = 1;
    % Sinus frequency
    sinFreq = 10;
    % Plotting times
    tau = -5*T:tauinc:5*T;
    % sampling times
    samplingTimes = 0:T:Length;

    signal = sin(sinFreq*samplingTimes).*exp(-3*samplingTimes);

    polynomialCoefficients = PolynomialInterpolationCoefficients(signal,T);

    maxLag = 20;
    PCCC = PolynomialCoefficientsCrossCorrelation(polynomialCoefficients,polynomialCoefficients,maxLag);
    
    [cc ccd ccdd] = CrossCorrelationInterpolation(PCCC,tau,T);
% 
    plot(tau,cc,'ob');
    hold on
    plot(tau,ccd,'gs');
    plot(tau,ccdd,'rd');
end