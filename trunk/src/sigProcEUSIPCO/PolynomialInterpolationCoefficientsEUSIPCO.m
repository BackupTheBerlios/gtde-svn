function polynomialCoefficients = PolynomialInterpolationCoefficientsEUSIPCO(signal,times)

% PolynomialInterpolationCoefficients Computes the coefficients of the polynomials
% corresponding to the interpolation of the discrete signal. Right now it is
% linear interpolation.
% 
%   polynomialCoefficients = PolynomialInterpolationCoefficients(signal,times) 
%   computes the coefficients of the linear interpolation of signal sampled
%   at times. The variable times may be the sampling times or the sampling
%   period.
%   
%   see also InterpolateCrossCorrelation

    %%% Input check
    if length(times) ~= 1 && length(times) ~= length(signal)
        error('times length need to be either the signal length or 1');
    end
    
    %%% Declare output variables
    polynomialCoefficients = zeros(length(signal)-1,2);
    
    %%% Compute the coefficients
    % Zero order
    polynomialCoefficients(:,1) = signal(1:end-1);
    % First order
    if length(times) == 1
        polynomialCoefficients(:,2) = (signal(2:end)-signal(1:end-1))/times;
    else
        polynomialCoefficients(:,2) = (signal(2:end)-signal(1:end-1))./(times(2:end)-times(1:end-1));
    end
    

end