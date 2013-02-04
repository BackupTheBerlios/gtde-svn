function [solverOutput, methodOptions] = ngTDE_DIP(PCCC, microphones, samplingPeriod, x0, methodOptions)

%Unconstrained time difference estimation
%
% USAGE: [TDE fVal cVal] = ngTDE_parallel(sigPCCC, microphones, samplingPeriod, x0)
%
% PARAMETERS:
%     PCCC ~ polynomial cross correlation of the signals' coefficients
%     microphones ~ positions of the microphones
%     samplingPeriod ~ sampling period of the signal
%     x0 ~ initialization points of the local minimizer
%     methodOptions
%
% RESULT:
%     The time delay estimates.
%
% DESCRIPTION:
%     Estimates the time difference of signals (or their Polynomial 
%     Coefficients Cross-Correlation) acquired by microphones 
%     sampled at samplingPeriod. The method is based on an 
%     unconstrained maximization of the continuous estimation of the cross-
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
        error('Usage [solverOutput, methodOptions] = ngTDE_DIP(PCCC, microphones, samplingPeriod, x0, methodOptions)');
    end
    if size(PCCC,1) ~= size(microphones,1)
        error('There should be as many signals as microphones.');
    end
    if size(microphones,2) ~= 3 && size(microphones,2) ~= 2
        error('The system works in 2D or in 3D.');
    end
    if samplingPeriod <= 0
        error('The sampling period should be positive.');
    end
    
    % Declare the objective function
    objFunction = @(x) gTDECriterion(x,PCCC,microphones,samplingPeriod);
    consFunction = @(x,z) constraints(x,z);

    % Run the parallel local optimization
    [solverOutput, methodOptions] = ipsolver_parallel(x0,objFunction,consFunction,methodOptions);
end

function [c, J, W] = constraints(TDEs,z)
    % This is a dummy restriction, -inf everywhere with 0 gradient and
    % hessian
    c = ones(1,size(TDEs,2));
    c = -c;
    % The computed J is the gradient, and the solver expects the
    % Jacobian...
    J = zeros(size(TDEs));
    if ~isempty(z)
        W = zeros(size(TDEs,1),size(TDEs,1),size(TDEs,2));
    end
end
