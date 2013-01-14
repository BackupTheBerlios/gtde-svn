function testInterpolateCrossCorrelation

% test the interpolation of the cross-correlation function

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

    close all; clear; clc

    % Scale parameters
    a = 0.1;
    c = a / sqrt(2);
    % Microphones in a regular tetrahedron
    MICS = [ a   0   -c ;...
            -a   0   -c ;...
             0    a   c ;...
             0   -a   c];

    % Generate potision
    positionOptions = struct('coordinateSystem','spherical');
    positionOptions.bounds = cell(1,3);
        positionOptions.bounds{1} = [1 5]; 
    positionOptions.numberOfIntervals = zeros(1,3);
        positionOptions.numberOfIntervals(1) = 4;
    positionOptions.useExtrema = true(1,3);
        positionOptions.useExtrema(3) = false;
    Positions = GeneratePositions(3,positionOptions);

    % Choose the signal
    F = 8000;
    L = 200;

    % Sinusoidal with exponential decay
%    myFun = @(x) SyntheticSignal(x,2,[F,L]);

    % Sampling frequency
    samplingPeriod = 1/48000;

    % Length
    Length = 0.01;

    positionIndex = 10;
    % Generate the signal
    signals = GenerateDiscreteSignals(Positions(positionIndex,:),...
                                      MICS,...
                                      samplingPeriod,...
                                      Length,...
                                      myFun);

    % Plotting period
    plottingPeriod = samplingPeriod/50;
    h = 0.9*plottingPeriod;
    % Plotting times
    maxTDE = TDEmax(MICS);
    maxTDE = maxTDE(1);
    maxLAG = ceil(maxTDE/samplingPeriod);
    maxTDE = maxLAG*samplingPeriod;
    plottingTimes = -maxTDE:plottingPeriod:maxTDE;
    % sampling times
    samplingTimes = 0:samplingPeriod:Length;
    [myCorr myDerivative my2Derivative] = InterpolateCrossCorrelation(signals(1,:),signals(2,:),plottingTimes,samplingPeriod,maxTDE);
%     myCorr = myCorr / max(myCorr);
%     myDerivative = myDerivative / max(myDerivative);
%     my2Derivative = my2Derivative / max(my2Derivative);

    % Discrete correlation
%     [discreteCorr lags] = xcorr(signals(1,:),signals(2,:),maxLAG,'unbiased');
%     discreteCorr = discreteCorr / max(discreteCorr);
    
    t1 = norm(Positions(positionIndex,:)-MICS(1,:),2)/343.2;
    t2 = norm(Positions(positionIndex,:)-MICS(2,:),2)/343.2;
    fprintf('TDE: %1.8e\n',-(t2-t1));
    
    aCC = analyticCrossCorrelation(F,1/L,Length,t1,t2,plottingTimes);
%     aCC = aCC / max(aCC);

    % First plot: signals
    subplot(2,2,1);
    plot(samplingTimes,signals(1,:));
    hold on
    plot(samplingTimes,signals(2,:),'r');
    hold off
    legend('signal1','signal2');
    
    % Second plot: cross-correlation
    subplot(2,2,2);
%     plot( -maxTDE:samplingPeriod:maxTDE, discreteCorr, 'ob');
    hold on
    plot( plottingTimes, myCorr, 'gx');
    plot(plottingTimes,aCC,'k');
    hold off
    legend('Mycorr','Analytic autocorrelation');
    
    % Thir plot: xcorr's derivative
    subplot(2,2,3);
%     dDerivative = discreteCorr(2:end)-discreteCorr(1:end-1);
%     dDerivative = dDerivative / max(dDerivative);
%     plot( - maxTDE:samplingPeriod:(maxTDE-samplingPeriod), dDerivative,'ob');
    hold on
    plot(plottingTimes, myDerivative,'gx');
    hold off
    legend('My derivative');
    
    % Fourth: 2nd derivative
    subplot(2,2,4);
%     ddDerivative = dDerivative(2:end)-dDerivative(1:end-1);
%     ddDerivative = ddDerivative / max(ddDerivative);
%     plot( - maxTDE:samplingPeriod:(maxTDE-2*samplingPeriod), ddDerivative,'ob');
    hold on
    plot(plottingTimes, my2Derivative,'gx');
    hold off
    legend('My 2nd derivative');
return

function acc = analyticCrossCorrelation(F,lambda,length,t1,t2,tau)
    
    M = min(length+t1,length+t2+tau);
    m = max(t1,t2+tau);

    acc = (lambda/4)*exp((t1+t2+tau)/lambda).*( ( exp(-2*M/lambda).*( cos(F*(2*M-t1-tau-t2)) - F*lambda*sin(F*(2*M-t1-tau-t2)) ) - ...
        exp(-2*m/lambda).*( cos(F*(2*m-t1-t2-tau)) - F*lambda*sin(F*(2*m-t1-t2-tau)) ) )/(1+(lambda*F).^2) - cos(F*(t2+tau-t1)).*( exp(-2*M/lambda) - exp(-2*m/lambda) ) );

return
