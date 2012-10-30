
import numpy as np
import scipy

from gtdecpy.src import K1
from gtdecpy.src import K2

def CrossCorrelationInterpolation(PCCC, Delay, SamplingPeriod, ZeroIndex=-1):

#CrossCorrelationInterpolation Cross-correlation function interpolation
#
# USAGE: [cc ccd ccdd] = CrossCorrelationInterpolation(PCCC,Delay,T[,ZI])
#
# PARAMETERS:
#  PCCC ~ cross-correlation signals of the polynomial coefficients.
#  Delay ~ delay is the point(s) in which we want to estimate the signals'
#     cross-correlation function.
#  T ~ the signals' sampling period.
#  ZI ~ the index on the PCCC representing no delay (the default value is ceil(L/2), where L is the
#     size of the coefficients cross-correlation function.)
# 
# RETURN VALUE:
#  cc,ccd,ccdd ~ The value of the cross-correlation function, its first and
#     its second derivative.
# 
# DESCRIPTION:
#     This function estimates the value of the cross-correlation function
#     at Delay from from the polynomial coefficients' cross-correlation 
#     signals PCCC.
#    
# REFERENCES:
#     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
#     estimation-based sound source localisation (gTDESSL). Research Report 
#     RR-7988, INRIA, June 2012.
#
#     see also K1, K2, InterpolateCrossCorrelation

# Copyright 2012, Xavier Alameda-Pineda
# INRIA Grenoble Rhone-Alpes
# E-mail: xavi.alameda@gmail.com
# 
# This is part of the gtde program.
# 
# gtde is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    # Input check
    if ZeroIndex == -1:
        ZeroIndex = np.ceil(PCCC.shape[2]/2.)
    ZeroIndex = int(ZeroIndex)
    
    # Output value
    cc = np.zeros(Delay.size)
    ccd = np.zeros(Delay.size)
    ccdd = np.zeros(Delay.size)
    
    # Index in the cross-correlation
    ccIndex = np.floor( Delay/SamplingPeriod ).astype(int)
    
    # Remaining of the division (0 <= tau < SamplingPeriod)
    tau = Delay - ccIndex*SamplingPeriod
    
#    for ii in np.arange(0,ccIndex.size):
#        print Delay[ii], ccIndex[ii], tau[ii]
    
    # Interpolation polynomial degrees
    D1 = PCCC.shape[0]-1.
    D2 = PCCC.shape[1]-1.
    
    # Loop on the degree
    for d1 in np.arange(0., (D1)+1):
#        print d1
        for d2 in np.arange(0., (D2)+1):
#                print d2
                # Value of the first interpolating function
                K1V = K1.K1(int(d1),int(d2),tau,SamplingPeriod)
                # Value of the second interpolation function
                K2V = K2.K2(int(d1),int(d2),tau,SamplingPeriod)
                # First cross-correlation value
                R1 = PCCC[d1,d2,ccIndex+ZeroIndex]
                # Second cross-correlation value
                R2 = PCCC[d1,d2,ccIndex-1+ZeroIndex]
                # Cross-correlation's value
                cc = cc + R1*K1V[0] + R2*K2V[0]
                # Cross-correlation's derivative
                ccd = ccd + R1*K1V[1] + R2*K2V[1]
                # Cross-correlation's second derivative
                ccdd = ccdd + R1*K1V[2] + R2*K2V[2]        

    return [cc, ccd, ccdd]