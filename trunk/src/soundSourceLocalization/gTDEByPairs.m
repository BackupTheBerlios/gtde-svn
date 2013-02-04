function [solverOutput] = gTDEByPairs(PCCC, microphones, samplingPeriod, x0, methodOptions)

%By-pair time difference estimation
%
% USAGE: [TDE fVal cVal] = gTDEByPairs(sigPCCC, microphones, samplingPeriod, x0)
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
%     sampled at samplingPeriod. The method maximizes each
%     cross-correlation function independently of signals' received at
%     other microphones.
% 
%   see also PolynomialInterpolationCoefficients


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
    
    %%% Input check
    if nargin < 4
        error('Usage TDE = gTDE(signals, microphones, samplingPeriod, x0[, verbose])');
    end
    if size(sigPCCC,1) ~= size(microphones,1)
        error('There should be as many signals as microphones.');
    end
    if size(microphones,2) ~= 3 && size(microphones,2) ~= 2
        error('The system works in 2D or in 3D.');
    end
    if samplingPeriod <= 0
        error('The sampling period should be positive.');
    end
    
    %%% Declare solver output
    if methodOptions.output > 0
        solverOutput.x = cell(NMics-1,1);
        solverOutput.f = cell(NMics-1,1);
    end
    TDE = zeros(1,NMics-1);
    
    %%% Look for the minima in each TDE
    for mic= 2:NMics,
        % Declare the objective function
        objFunction = @(x) CrossCorrelationInterpolation(PCCC{1,mic},x,samplingPeriod);
        % Evaluate at X0
        f = objFunction(x0{mic-1}');
        % Compute the minimum now
        [~,index] = max(f);
        TDE(mic-1) = x0{mic-1}(index);
        if methodOptions.output > 0
            % Store other results
            outputSolver.x{mic-1} = x0{mic-1};
            outputSolver.f{mic-1} = f;
        end
    end
    
    % Save optimal
    outputSolver.xstar = TDE;
end

% function f = objective(PCCC,x,samplingPeriod)
%     % Compute
%     f = CrossCorrelationInterpolation(PCCC,x,samplingPeriod);
% end