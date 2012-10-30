
import numpy as np
import scipy
import sys
sys.path.append("/local_scratch/alamedap/Software/gtde/trunk")
from gtdecpy.src import DiscreteCrossCorrelation

def PolynomialCoefficientsCrossCorrelation(firstCoefficients, secondCoefficients, maxLag=-1):

#Compute the cross-correlation of two sets of coefficients.
#
# USAGE: PCCC = PolynomialCoefficientsCrossCorrelation(firstCoefficients,secondCoefficients,MaxLag)
#
# PARAMETERS:
#   firstCoefficients ~ set of coefficients of the first signal
#   secondCoefficients ~ set of coefficients of the second signal
#   MaxLag ~ maximum lag value to take into account when computing the cross
#     correlation
#
# RETURN VALUE:
#   PCCC ~ discrete cross-correlation of the polynomial coefficients.
# 
# DESCRIPTION:
#     Outputs a (D1+1)x(D2+1)x(2*MaxLag+1) matrix, where D1 and D2 are the 
#     degrees of the polynomial interpolations and MaxLag is the size of the 
#     desired output correlation. The default value for MaxLag is NSamples-1, 
#     where NSamples is the length of firstCoefficients (or secondCoefficients, 
#     which must be the same).
#
# REFERENCES:
#     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
#     estimation-based sound source localisation (gTDESSL). Research Report 
#     RR-7988, INRIA, June 2012.
#
#   see also InterpolateCrossCorrelation, DiscreteCrossCorrelation

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
    
    ### Input check
    if firstCoefficients.ndim != 2:
        raise Exception("firstCoefficients should be 2D.")
    if secondCoefficients.ndim != 2:
        raise Exception("firstCoefficients should be 2D.")    
    
    ### General variables
    # Compute the output correlation size
    if maxLag == -1:
        maxLag = min(firstCoefficients.shape[1],secondCoefficients.shape[1])-1.
    
    
    # Polynomial degress
    M1 = firstCoefficients.shape[0]
    M2 = secondCoefficients.shape[0]
    # Declare output variable
    PCCC = np.zeros( (M1, M2, (2.*maxLag+1.)) )
    ### Computation
    for m1 in np.arange(1., (M1)+1):
        for m2 in np.arange(1., (M2)+1):
            PCCC[int(m1)-1,int(m2)-1,:] = DiscreteCrossCorrelation.DiscreteCrossCorrelation(firstCoefficients[int(m1)-1,:], secondCoefficients[int(m2)-1,:], maxLag)            
        
    return PCCC