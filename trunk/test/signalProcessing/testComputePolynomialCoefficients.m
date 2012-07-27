function testComputePolynomialCoefficients

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
    plottingPeriod = 1e-3;
    % Samplingperiod
    samplingPeriod = 1e-2;
    % Length
    Length = 0.5;
    % Sinus frequency
    sinFreq = 100;
    % Plotting times
    plottingTimes = 0:plottingPeriod:Length;
    % sampling times
    samplingTimes = 0:samplingPeriod:Length;

    signal = sin(sinFreq*samplingTimes);

    polynomialCoefficients = PolynomialInterpolationCoefficients(signal,samplingTimes);

    interpolatedSignal = InterpolateSignal(polynomialCoefficients,samplingTimes,plottingTimes);

    plot(plottingTimes,interpolatedSignal);
    hold on
    plot(samplingTimes,signal,'or');
    plot(samplingTimes(1:end-1),polynomialCoefficients(:,1),'gx');
    hold off
    
end