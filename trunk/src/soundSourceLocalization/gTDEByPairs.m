function [TDE fVal cVal] = gTDEByPairs(sigPCCC, microphones, samplingPeriod, x0)

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
    
    %%% Look for the minima in each TDE
%     maxTDE = TDEmax(microphones);
    TDE = zeros(size(1,NMics-1));
    fVal = zeros(size(1,NMics-1));
    cVal = zeros(size(1,NMics-1));
    for mic= 2:NMics,
        objFunction = @(x) objective(PCCC{1,mic},x,samplingPeriod);
%         consFunction = @(x,z) constraints(x,z);
%         [t f] = ipsolver_parallel(x0{mic-1}'/samplingPeriod,objFunction,consFunction,'newton',1e-8,15,false,true);
%         % Trash out those who are out of range
%         indices = t > maxTDE(mic-1) | t < -maxTDE(mic-1);
%         f(indices) = max(f);
        f = objFunction(-x0{mic-1}'/samplingPeriod);
        % Compute the minimum now
        [~,index] = max(f);
        TDE(mic-1) = x0{mic-1}(index);
        fVal(mic-1) = f(index);
        cVal(mic-1) = 1;
    end

end

function [f g H] = objective(PCCC,x,samplingPeriod)
    x = x*samplingPeriod;
    % Compute
    [f gp Hp] = CrossCorrelationInterpolation(PCCC,x,samplingPeriod);
    % Allocate as the need of the algorithm
    g = -gp*(samplingPeriod);
    H = zeros(size(Hp,1),size(Hp,1),size(Hp,2));
    H(1,1,:) = Hp*(samplingPeriod^2);
end

function [c J W] = constraints(TDEs,z)
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