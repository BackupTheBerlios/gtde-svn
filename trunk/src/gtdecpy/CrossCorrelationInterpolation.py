
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def CrossCorrelationInterpolation(PCCC, Delay, SamplingPeriod, ZeroIndex):

    # Local Variables: ccIndex, R1, d2, tau, cc, PCCC, k1dd, ZeroIndex, R2, k1, ccdd, Delay, k2, SamplingPeriod, k1d, ccd, d1, D2, k2dd, k2d, D1
    # Function calls: floor, ceil, nargin, K1, CrossCorrelationInterpolation, zeros, error, squeeze, K2, size
    #%ccInterpolation Cross-correlation function interpolation
    #%
    #% USAGE: [cc ccd ccdd] = CrossCorrelationInterpolation(PCCC,Delay,T[,ZI])
    #%
    #% PARAMETERS:
    #%  PCCC ~ cross-correlation signals of the polynomial coefficients.
    #%  Delay ~ delay is the point(s) in which we want to estimate the signals'
    #%     cross-correlation function.
    #%  T ~ the signals' sampling period.
    #%  ZI ~ the index on the PCCC representing no delay (the default value is ceil(L/2), where L is the
    #%     size of the coefficients cross-correlation function.)
    #% 
    #% RETURN VALUE:
    #%  cc,ccd,ccdd ~ The value of the cross-correlation function, its first and
    #%     its second derivative.
    #% 
    #% DESCRIPTION:
    #%     This function estimates the value of the cross-correlation function
    #%     at Delay from from the polynomial coefficients' cross-correlation 
    #%     signals PCCC.
    #%    
    #% REFERENCES:
    #%     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
    #%     estimation-based sound source localisation (gTDESSL). Research Report 
    #%     RR-7988, INRIA, June 2012.
    #%
    #%     see also K1, K2, InterpolateCrossCorrelation
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
    #% Input check
    if nargin<3.:
        matcompat.error('Usage: [cc ccd ccdd] = CCInterpolation(PCCC,Delay,SamplingPeriod[,ZeroIndex])')
    
    
    if nargin<4.:
        ZeroIndex = np.ceil((matcompat.size(PCCC, 3.)/2.))
    
    
    #% Output value
    cc = np.zeros(matcompat.size(Delay))
    ccd = np.zeros(matcompat.size(Delay))
    ccdd = np.zeros(matcompat.size(Delay))
    #% Index in the cross-correlation
    ccIndex = np.floor(matdiv(Delay, SamplingPeriod))
    #% Remaining of the division (0 <= tau < SamplingPeriod)
    tau = Delay-np.dot(ccIndex, SamplingPeriod)
    #% Interpolation polynomial degrees
    D1 = matcompat.size(PCCC, 1.)-1.
    D2 = matcompat.size(PCCC, 2.)-1.
    for d1 in np.arange(0., (D1)+1):
        for d2 in np.arange(0., (D2)+1):
            #%            try
            
        
    return [cc, ccd, ccdd]