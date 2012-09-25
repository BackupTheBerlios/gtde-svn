
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def DiscreteCrossCorrelation(signal1, signal2, maxLag):

    # Local Variables: xcPos, maxLag, chunkSize, signal2, signal1, corr
    # Function calls: DiscreteCrossCorrelation, sum, nargin, length, zeros, error
    #%Discrete cross-correlation
    #%
    #% USAGE: Cross-Correlation = DiscreteCrossCorrelation(FirstSignal,SecondSignal,MaxLag)
    #%
    #% PARAMETERS:
    #%   FirstSignal ~ values of the first signal
    #%   SecondSignal ~ values of the second signal
    #%   MaxLag ~ the maximum value of the lag we want to compute
    #%
    #% RETURN VALUE:
    #%   Cross-Correlation ~ values of the cross-correlation signal
    #% 
    #% DESCRIPTION:
    #%   This function computes the discrete cross-corrlation of two 1D signals.
    #%
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
    #%%% Check input
    if nargin<3.:
        matcompat.error('Usage: Cross-Correlation = DiscreteCrossCorrelation(FirstSignal,SecondSignal,MaxLag)')
    
    
    #%%% Initialize output
    corr = np.zeros(1., (2.*maxLag+1.))
    #% How many samples we use to compute each of the cross-correlation
    #% values
    chunkSize = length(signal1)-maxLag
    #%%% Computation
    #% Positive values of the delay
    for xcPos in np.arange(2., (maxLag+1.)+1):
        corr[int((maxLag+xcPos))-1] = np.sum((signal1[int(xcPos)-1:chunkSize-1.+xcPos]*signal2[0:chunkSize]))
        
    #% Negative values of the delay
    for xcPos in np.arange(2., (maxLag+1.)+1):
        corr[int((maxLag+2.-xcPos))-1] = np.sum((signal2[int(xcPos)-1:chunkSize-1.+xcPos]*signal1[0:chunkSize]))
        
    #% At 0
    corr[int((maxLag+1.))-1] = np.sum((signal1[0:chunkSize]*signal2[0:chunkSize]))
    return [corr]