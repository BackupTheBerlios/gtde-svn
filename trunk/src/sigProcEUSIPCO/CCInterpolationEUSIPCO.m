function [ccValue, ccd, ccdd] = CCInterpolationEUSIPCO(PCCC,Delay,SamplingPeriod,ZeroIndex)

% ccValue = CCInterpolation(PCCC,Delay,T,ZI)
% 
%     This function estimates the value of the cross-correlation function from
%     a polynomial interpolation of the two signals. PCCC is the cross-correlation
%     function of the polynomial coefficients from the two signals (the output
%     of the computeCoefficientsCrossCorrelation). Delay is he point at which
%     we want to estimate the signals' cross-correlation function. T is the
%     signals' sampling period. ZI is the index representing the position in 
%     the PCCC with no delay (the default value is ceil(L/2), where L is the
%     size of the coefficients cross-correlation function.
%     
%     see also K1, K2, InterpolateCrossCorrelation

    % Input check
    if nargin < 3
        error('Usage: ccValue = CCInterpolation(PCCC,Delay,SamplingPeriod[,ZeroIndex])');
    end
    if nargin < 4
        ZeroIndex = ceil(size(PCCC,3)/2);
    end

    % Output value
    ccValue = zeros(size(Delay));
    ccd = zeros(size(Delay));
    ccdd = zeros(size(Delay));
    
    % Index in the cross-correlation
    ccIndex = floor(Delay/SamplingPeriod);
    % Remaining of the division (0 <= tau < SamplingPeriod)
    tau = Delay - ccIndex*SamplingPeriod;
    
    % Interpolation polynomial degrees
    D1 = size(PCCC,1)-1;
    D2 = size(PCCC,2)-1;
    
    for d1 = 0:D1,
        for d2 = 0:D2,
            try
                % Value of the first interpolating function
                [k1, k1d, k1dd] = K1EUSIPCO(d1,d2,tau,SamplingPeriod);
                % Value of the second interpolation function
                [k2, k2d, k2dd] = K2EUSIPCO(d1,d2,tau,SamplingPeriod);
                % First cross-correlation value
                R1 = (squeeze(PCCC(d1+1,d2+1,ccIndex+1+ZeroIndex))');
                % Second cross-correlation value
                R2 = (squeeze(PCCC(d1+1,d2+1,ccIndex+ZeroIndex))');
                % Cross-correlation's value
                ccValue = ccValue + R1.*k1 + R2.*k2;
                % Cross-correlation's derivative
                ccd = ccd + R1.*k1d + R2.*k2d;
                % Cross-correlation's second derivative
                ccdd = ccdd + R1.*k1dd + R2.*k2dd;
            catch err
%                 disp('Exception');
            end
        end
    end

end