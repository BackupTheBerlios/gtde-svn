function [TDE fVal cVal] = gTDEByPairs(sigPCCC, microphones, samplingPeriod, x0)

% gXCE geometric-constrained time difference estimation
% 
%     TDE = gTDE(sigPCCC, microphones, samplingPeriod, x0, verbose)
%     Estimates the time difference of signals (or their Polynomial 
%     Coefficients Cross-Correlation) acquired by microphones 
%     sampled at samplingPeriod. The method is based on a geometric-
%     constrained maximization of the continuous estimation of the cross-
%     correlation functions. 'signals' is a M-by-N matrix with one row per 
%     acquired signal. 'microphoned' is a M-by-3 matrix with one row per 
%     microphone's position. TDE is a M(M-1)/2 vector with the TDEs. 
%
%     The reference microphone will be the first one, although this 
%     should not matter.
%
%     see also gTDESSL, PolynomialInterpolationCoefficients,
%     PolynomialCoefficientsCrossCorrelation
    
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