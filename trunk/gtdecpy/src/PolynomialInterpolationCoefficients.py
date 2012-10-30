
import numpy as np
import scipy

def PolynomialInterpolationCoefficients(signal, step):

#Compute the polynomial interpolation coefficients
#
# USAGE: polynomialCoefficients = PolynomialInterpolationCoefficients(signal,step)
#
# PARAMETERS:
#  signal ~ the value of the sampled signal
#  step ~ the sampling step
#
# RETURN VALUE:
#  polynomialCoefficients ~ the polynomial coefficients of the interpolated
#     signal
# 
# DESCRIPTION:
#     This function computes the coefficients of the linear interpolation
#     of the signal
#
# REFERENCES:
#     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
#     estimation-based sound source localisation (gTDESSL). Research Report 
#     RR-7988, INRIA, June 2012.
#
#   see also InterpolateCrossCorrelation

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
 
    ### Declare output, cubic interpolation means 4 coefficients
    polynomialCoefficients = np.zeros((4,signal.size-1))
    # Signal size shorcut
    L = signal.size
    
    ### State and solve the linear system for the natural splines
    # Independent vector
    bVector = signal[0:L-2] - 2*signal[1:L-1] + signal[2:L]
    bVector *= 6./step**2
    # Matrix
    AMatrix = np.zeros((L-2,L-2))
    for ii in np.arange(0,L-3):
        AMatrix[ii,ii+1] = 1
        AMatrix[ii+1,ii] = 1
        AMatrix[ii,ii] = 4
    AMatrix[L-3,L-3] = 4
    # Solve the system
    mVector = np.dot(np.linalg.inv(AMatrix),bVector)
    # Add the two zeros to produce natural splines
    mVector = np.concatenate(([0], mVector, [0]))
    
    ### Compute the spline coefficients
    # Zero order
    polynomialCoefficients[0,:] = signal[0:signal.size-1]
    # First order
    polynomialCoefficients[1,:] = (signal[1:L]-signal[0:L-1])/step - (mVector[1:L]-mVector[0:L-1])*step/6
    # Second order
    polynomialCoefficients[2,:] = mVector[0:L-1]/2
    # Third order
    polynomialCoefficients[3,:] = (mVector[1:L]-mVector[0:L-1])/(6*step)
        
    ### Return value
    return polynomialCoefficients