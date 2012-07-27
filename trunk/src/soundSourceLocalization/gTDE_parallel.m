function [TDE fVal cVal] = gTDE_parallel(sigPCCC, microphones, samplingPeriod, x0)

% gTDE geometric-constrained time difference estimation
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
    
    % Declare the objective function
%     objFunction = @(x) gTDECriterion(x,PCCC,microphones,samplingPeriod);
    objFunction = @(x) gTDELogCriterion(x,PCCC,microphones,samplingPeriod);
    consFunction = @(x,z) constraints(x,z,microphones,samplingPeriod);
%     hessFunction = @(x,lambda) gTDEHessianFunction(x,lambda,PCCC,microphones,samplingPeriod,true);
    
    % if verbose, declare verbose options
%     if verbose
%         options = optimset('GradObj','on',...
%                            'GradConstr','on',...
%                            'Hessian','user-supplied',...
%                            'Algorithm','interior-point',...
%                            'DerivativeCheck','on',...
%                            'HessFcn',hessFunction,...
%                            'TolX',1e-14,...
%                            'TolFun',1e-8,...
%                            'PlotFcns',{...
%                                @optimplotx,...
%                                @optimplotfunccount,...
%                                @optimplotfval,...
%                                @optimplotconstrviolation,...
%                                @optimplotstepsize,...
%                                @optimplotfirstorderopt...
%                                @myPlotFunction...                                   
%                                },...
%                            'Display','on'...
%                            );        
%     else
%         options = optimset('GradObj','on',...
%                            'GradConstr','on',...
%                            'Hessian','user-supplied',...
%                            'Algorithm','interior-point',...
%                            'HessFcn',hessFunction,...
%                            'TolX',1e-14,...
%                            'TolFun',1e-8,...
%                            'MaxIter',10,...
%                            'Display','off'...
%                            );
%     end

    % Run the optimization method
%     [TDE fVal] = My_fmincon(objFunction, x0,...
%                                   [], [], [], [],...
%                                   -maxTDE(1:NMics-1), maxTDE(1:NMics-1),... 
%                                   consFunction,...
%                                   options);
    [TDE fVal cVal] = ipsolver_parallel(x0,objFunction,consFunction,'newton',1e-8,15,false);
end

function [c J W] = constraints(TDEs,z,microphones,samplingPeriod)
    % Compute the restriction (is a positive restriction)
    [c J W] = TDEDiscriminant(TDEs,microphones,samplingPeriod);
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
