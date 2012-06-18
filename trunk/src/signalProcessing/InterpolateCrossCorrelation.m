function [CrossCorrelation Derivative Curvature] = InterpolateCrossCorrelation(FirstSignal,SecondSignal,Delays,SamplingPeriod,MaxDelay)

% InterpolatieCrossCorrelation interpolates the cross-correlation of two
% discrete signals
% 
%   CrossCorrelation = InterpolateCrossCorrelation(FirstSignal,SecondSignal,Delays,SamplingPeriod)
% 
%   This function interpolates the cross-correlation function at times Delays
%   of the signals FirstSignal and SecondSignal (sampled at SamplingPeriod). 
%   It does that assuming some polynomial interpolation at each time interval.
% 
%   This function will then compute:
%   1) The sequence of coefficients of each polynomial interpolation.
%   2) The cross-correlation of these sequences of coefficients.
%   3) The interpolation of the cross-correlation function at each of the Delays.
% 
%   see also PolynomialInterpolationCoefficients, PolynomialCoefficientsCrossCorrelation and performCCInterpolation

    %%% Input check
    if nargin < 4
        error('Usage: [CrossCorrelation Derivative] = Interpolation(FirstSignal,SecondSignal,Delays,SamplingPeriod[,MaxDelay]');
    end

    %%% Compute the interpolation polynomial coefficients of the signals
    [FirstPolynomialCoefficients] = PolynomialInterpolationCoefficients(FirstSignal,...
                                                                        SamplingPeriod);
    [SecondPolynomialCoefficients] = PolynomialInterpolationCoefficients(SecondSignal,...
                                                                         SamplingPeriod);
    
    %%% Compute the cross-correlation of the coefficients
    if nargin < 5
        [PCCC] = PolynomialCoefficientsCrossCorrelation(FirstPolynomialCoefficients,...
                                                        SecondPolynomialCoefficients);
    else
        MaxLag = ceil(MaxDelay/SamplingPeriod);
        [PCCC] = PolynomialCoefficientsCrossCorrelation(FirstPolynomialCoefficients,...
                                                        SecondPolynomialCoefficients,...
                                                        MaxLag);
    end
                                                        
    
    %%% Interpolate the cross-correlation at "Delays"
    CrossCorrelation = zeros(size(Delays));
    Derivative = zeros(size(Delays));
    Curvature = zeros(size(Delays));
    
    for dd = 1:numel(Delays),
        if dd == 2651,
            disp('stop');
        end
        [CrossCorrelation(dd) Derivative(dd) Curvature(dd)] = CCInterpolation(PCCC,...
                                                                Delays(dd),...
                                                                SamplingPeriod);
    end
    
end
    