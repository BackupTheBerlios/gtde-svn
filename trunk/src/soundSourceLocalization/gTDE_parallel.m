function [solverOutput] = gTDE_DIP(sigPCCC, microphones, samplingPeriod, x0)

%Geometric-constrained time difference estimation
%
% USAGE: [TDE fVal cVal] = gTDE_DIP(sigPCCC, microphones, samplingPeriod, x0)
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
    
    %%% Input check
    if nargin < 4
        error('Usage TDE = gTDE_DIP(signals, microphones, samplingPeriod, x0[, verbose])');
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
    
    %%% General variables
    % Number of microphones
    NMics = size(microphones,1);
    % Maximum TDE values from the microphones
    maxTDE = TDEmax(microphones)/samplingPeriod;
    
    % Check if the input are signals or PCCC
    if ~iscell(sigPCCC)
        % Number of samples
        NSamples = size(sigPCCC,2);
        % Sampling times
        samplingTimes = 0:samplingPeriod:(NSamples-1)*samplingPeriod;
        % Maximum cross-correlation lag
        maxTDESamples = ceil(1.1*maxTDE/samplingPeriod);
        maxLAG = 2*max(maxTDESamples);
        
        % Compute the interpolation coefficients
        PC = cell(NMics,1);
        for ss = 1:NMics,
            PC{ss} = PolynomialInterpolationCoefficients(sigPCCC(ss,:),samplingTimes);
        end

        % Compute the Cross-correlation of the interpolation coefficients
        PCCC = cell(NMics);
        for mic1 = 1:NMics,
            for mic2 = mic1:NMics,
                PCCC{mic1,mic2} = PolynomialCoefficientsCrossCorrelation(PC{mic1},PC{mic2});
            end
        end
    else
        PCCC = sigPCCC;
    end
    
    % Declare the objective function
%     objFunction = @(x) gTDECriterion(x,PCCC,microphones,samplingPeriod);
    objFunction = @(x) gTDECriterion(x,PCCC,microphones,samplingPeriod);
    consFunction = @(x,z) constraints(x,z,microphones,samplingPeriod);

    % Run the parallel local optimization
    [solverOutput, experimentOptions.methodOptions] = ipsolver_parallel(x0,objFunction,consFunction,experimentOptions.methodOptions);
    
    
end

function [c, J, W] = constraints(TDEs,z,microphones,samplingPeriod)
    % Compute the restriction (is a positive restriction)
    [c, J, W] = TDEDiscriminant(TDEs,microphones,samplingPeriod);
    % Change the sign of the restriction and its gradient, and set the
    % equatlity constratints to null
    c = -c;
    % The computed J is the gradient, and the solver expects the
    % Jacobian...
    J = -J;
    if ~isempty(z)
        W = -z(1)*W;
    end
end
