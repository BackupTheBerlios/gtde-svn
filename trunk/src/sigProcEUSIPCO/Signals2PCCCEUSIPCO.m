function [PCCC] = Signals2PCCCEUSIPCO(signals,MICS,samplingFrequency)

%Geometric-constrained time difference estimation
%
% USAGE:  PCCC = Signals2PCCC(signals,MICS,samplingFrequency)
%
% PARAMETERS:
%     sigPCCC ~ polynomial cross correlation of the signals' coefficients
%     microphones ~ positions of the microphones
%     samplingPeriod ~ sampling period of the signal
%     x0 ~ initialization points of the local minimizer
%
% RESULT:
%     The time delay estimates.
%
% DESCRIPTION:
%     Estimates the time difference of signals (or their Polynomial 
%     Coefficients Cross-Correlation) acquired by microphones 
%     sampled at samplingPeriod. The method is based on a geometric-
%     constrained maximization of the continuous estimation of the cross-
%     correlation functions. 'signals' is a M-by-N matrix with one row per 
%     acquired signal. 'microphoned' is a M-by-3 matrix with one row per 
%     microphone's position. TDE is a M(M-1)/2 vector with the TDEs. 
% 
%   see also PolynomialInterpolationCoefficients, ipsolver_parallel, gTDELogCriterion


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

    % Check input
    if nargin < 3
        error('Usage:  [PCCC EPCCC] = Signals2PCCC(signals,MICS,samplingFrequency).');
    end
    
    % Number of microphones
    NMics = size(signals,1);
    
    % Compute maxlag
    maxLag = TDEmax(MICS);
    maxLag = ceil(samplingFrequency*max(maxLag));
    % Allocate
    PC = cell(NMics,1);
    % Compute
    for ss = 1:NMics,
        PC{ss} = PolynomialInterpolationCoefficientsEUSIPCO(signals(ss,:),1/samplingFrequency);
%         figure;
%         subplot(5,1,1), plot(signals(ss,:));
%         subplot(5,1,2), plot(PC{ss}(:,1));
%         subplot(5,1,3), plot(PC{ss}(:,2));
%         subplot(5,1,4), plot(PC{ss}(:,3));
%         subplot(5,1,5), plot(PC{ss}(:,4));
    end
    % Compute the Cross-correlation of the interpolation coefficients
    % Allocate
    PCCC = cell(NMics);
    % Compute
    for mic1 = 1:NMics,
        for mic2 = mic1:NMics,
            PCCC{mic1,mic2} = PolynomialCoefficientsCrossCorrelationEUSIPCO(PC{mic1},PC{mic2},maxLag);
        end
    end
    
end