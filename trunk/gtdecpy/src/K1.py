
import numpy as np
import scipy
import gmpy as gp

def K1(p, q, tau, T):

#K1 function
#
# USAGE: [k1 k1d k1dd] = K1(p,q,tau,T)
#
# PARAMETERS:
#  p,q ~ Degree of the interpolation polynomial of the first and second
#       signals.
#  tau ~ Point (modulo T) in which the cross-correlation function is 
#       evaluated.
#  T ~ Sampling period of the discrete signals.
#
# RETURN VALUE:
#  k1,k1d,k1dd ~ Value of the function, the first and the second derivative
#       at tau respectively.
# 
# DESCRIPTION:
#     The K1 function corresponds to the following formula
#     k1 = sum_{p=0}^q choose(q,k) * (T-tau)^(q-k) * tau^(p+k+1) / (p+k+1)
#     and its derivatives (1st and 2nd) wigh respect to tau. T and tau may 
#     be two vectors of the same size, or vector and number.
#
# REFERENCES:
#     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
#     estimation-based sound source localisation (gTDESSL). Research Report 
#     RR-7988, INRIA, June 2012.
#
#     see also K2, CCInterpolation
# Copyright 2012, Xavier Alameda-Pineda
# INRIA Grenoble Rhone-Alpes
# E-mail: xavi.alameda@gmail.com
 
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

    ### Compute the function's value
    k1 = np.zeros((1,tau.size))
    for r in np.arange(0., q+1):
        k1 = k1 + gp.comb(q, int(r)) * (T-tau)**(q-r) * tau**(p+r+1.) / float(p+r+1.)
        
    ### Compute its first derivative
    k1d = np.zeros((1,tau.size))
    for r in np.arange(0., q):
        k1d = k1d + gp.comb(q, int(r)) * (T-tau)**(q-r-1.) * tau**(r+p) * ( T-(1.+p+q)*tau / float(r+p+1.) )
    # r = q may give some numerical problems, special formula
    k1d = k1d + tau**(q+p)
    
    
    ### Compute its second derivative
    k1dd = np.zeros((1,tau.size))
    # -----
    # r = 0
    # -----
    if p >= 1.:
        if q >= 2.:
            # Same formula
            r = 0.
            k1dd = k1dd + gp.comb(q, int(r)) * (T-tau)**(q-r-2.) * tau**(r+p-1.) * ( (r+p)*T**2. - 2.*(q+p)*T*tau + (q+p)*(q+p+1)*tau**2. / float(r+p+1.) )
        elif q == 1.:
            k1dd = k1dd - tau**(p-1) * ( p*(tau-T) + 2*tau )
            
        else:
            k1dd = k1dd + p*tau**(p-1)
    else:
        # Special cases
        if q >= 2.:
            k1dd = k1dd + q*(T-tau)**(q-2) * ( (q+1)*tau - 2*T )
        elif q == 1.:
            k1dd = k1dd-2.
            
    # -------
    # r = q-1, exists if q >= 1, but for q = 1, the r = q will take
    # care of it
    # -------
    if q >= 2.:
        k1dd = k1dd + q*tau**(q+p-2) * ( (q+p)*(T-tau) - (T+tau) )
    
    # -----
    # r = q, it exists for q >= 1, since for q = 0, r = q will take
    # care of it
    # -----
    if q >= 1.:
        k1dd = k1dd + (q+p)*tau**(q+p-1)
        
    # Intermediate terms
    for r in np.arange(1., (q-2.)+1):
        k1dd = k1dd + gp.comb(q, int(r)) * (T-tau)**(q-r-2.) * tau**(r+p-1.) * ( float((r+p)*T**2.) - 2*(q+p)*T*tau + (q+p)*(q+p+1.)*tau**2/float(r+p+1) )
            
    return [k1, k1d, k1dd]