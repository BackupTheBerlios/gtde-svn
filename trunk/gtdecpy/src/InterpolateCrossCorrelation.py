
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def InterpolateCrossCorrelation(FirstSignal, SecondSignal, Delays, SamplingPeriod, MaxDelay):

    # Local Variables: SecondPolynomialCoefficients, FirstPolynomialCoefficients, MaxDelay, SecondSignal, MaxLag, PCCC, Delays, Curvature, SamplingPeriod, CrossCorrelation, FirstSignal, Derivative, dd
    # Function calls: disp, PolynomialCoefficientsCrossCorrelation, numel, PolynomialInterpolationCoefficients, ceil, nargin, zeros, InterpolateCrossCorrelation, error, size
    #%InterpolatieCrossCorrelation interpolates the cross-correlation function 
    #%of two discrete signals
    #%
    #% USAGE: [CrossCorrelation Derivative Curvature] = 
    #%   InterpolateCrossCorrelation(FirstSignal,SecondSignal,Delays,SamplingPeriod)
    #%
    #% PARAMETERS:
    #%   FirstSignal ~ values of the first signal
    #%   SecondSignal ~ values of the second signal
    #%   Delays ~ delays in which we want to evaluate the cross-correlation function
    #%   SamplingPeriod ~ sampling period of the signals
    #%
    #% RETURN VALUE:
    #%   CrossCorrelation ~ values of the cross-correlation function at delays
    #%   Derivatives ~ values of its derivative at delays
    #%   Curvature ~ values of its second derivative at delays
    #% 
    #% DESCRIPTION:
    #%   This function interpolates the cross-correlation function at times Delays
    #%   of the signals FirstSignal and SecondSignal (sampled at SamplingPeriod). 
    #%   It does that assuming some polynomial interpolation at each time interval.
    #% 
    #%   This function will then compute:
    #%   1) The sequence of coefficients of each polynomial interpolation.
    #%   2) The cross-correlation of these sequences of coefficients.
    #%   3) The interpolation of the cross-correlation function at each of the Delays.
    #%
    #% REFERENCES:
    #%     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
    #%     estimation-based sound source localisation (gTDESSL). Research Report 
    #%     RR-7988, INRIA, June 2012.
    #%
    #%   see also PolynomialInterpolationCoefficients, PolynomialCoefficientsCrossCorrelation and performCCInterpolation
    #% Copyright 2012, Xavier Alameda-Pineda
    #% INRIA Grenoble Rhone-Alpes
    #% E-mail: xavi.alameda@gmail.com
    #% 
    #% This is part of the gtde program.
    #% 
    #% gtde is free software: you can redistribute it and/or modify
    #% it under the terms of the GNU General Public License as published by
    #% the Free Software Foundation, either version 3 of the License, or
    #% (at your option) any later version.
    #% 
    #% This program is distributed in the hope that it will be useful,
    #% but WITHOUT ANY WARRANTY; without even the implied warranty of
    #% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #% GNU General Public License for more details.
    #% 
    #% You should have received a copy of the GNU General Public License
    #% along with this program.  If not, see <http://www.gnu.org/licenses/>.
    #%%% Input check
    if nargin<4.:
        matcompat.error('Usage: [CrossCorrelation Derivative] = Interpolation(FirstSignal,SecondSignal,Delays,SamplingPeriod[,MaxDelay]')
    
    
    #%%% Compute the interpolation polynomial coefficients of the signals
    [FirstPolynomialCoefficients] = PolynomialInterpolationCoefficients(FirstSignal, SamplingPeriod)
    [SecondPolynomialCoefficients] = PolynomialInterpolationCoefficients(SecondSignal, SamplingPeriod)
    #%%% Compute the cross-correlation of the coefficients
    if nargin<5.:
        [PCCC] = PolynomialCoefficientsCrossCorrelation(FirstPolynomialCoefficients, SecondPolynomialCoefficients)
    else:
        MaxLag = np.ceil(matdiv(MaxDelay, SamplingPeriod))
        [PCCC] = PolynomialCoefficientsCrossCorrelation(FirstPolynomialCoefficients, SecondPolynomialCoefficients, MaxLag)
        
    
    #%%% Interpolate the cross-correlation at "Delays"
    CrossCorrelation = np.zeros(matcompat.size(Delays))
    Derivative = np.zeros(matcompat.size(Delays))
    Curvature = np.zeros(matcompat.size(Delays))
    for dd in np.arange(1., (numel(Delays))+1):
        if dd == 2651.:
            np.disp('stop')
        
        
        Delays[int(dd)-1]
        
    return [CrossCorrelation, Derivative, Curvature]