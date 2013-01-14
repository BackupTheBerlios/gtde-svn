function [k2,k2d,k2dd] = K2(p,q,tau,T)

    // Output variables initialisation (not found in input variables)
    k2=[];
    k2d=[];
    k2dd=[];
    
    // Number of arguments in function call
    [%nargout,%nargin] = argn(0)
    
    // Display mode
    mode(0);
    
    // Display warning for floating point exception
    ieee(1);
    
    
    //K2   K2 function
    // 
    // USAGE: [k2 k2d k2dd] = K2(p,q,tau,T)
    // 
    // PARAMETERS:
    //  p,q ~ Degree of the interpolation polynomial of the first and second
    //       signals.
    //  tau ~ Point (modulo T) in which the cross-correlation function is 
    //       evaluated.
    //  T ~ Sampling period of the discrete signals.
    // 
    // RETURN VALUE:
    //  k2,k2d,k2dd ~ Value of the function, the first and the second derivative
    //       at tau respectively.
    // 
    // DESCRIPTION:
    //     The K2 function corresponds to the following formula
    //     K2 = sum_(k=0).entries^q choose(q,k) * (tau)^(q-k) * ( T^(k+p+1) - tau^(k+p+1)) / (k+p+1)
    //     and its derivatives (1st and 2nd) wigh respect to tau. T and tau may 
    //     be two vectors of the same size, or vector and number.
    // 
    // REFERENCES:
    //     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
    //     estimation-based sound source localisation (gTDESSL). Research Report 
    //     RR-7988, INRIA, June 2012.
    // 
    //     see also K1, CCInterpolation
    
    // Copyright 2012, Xavier Alameda-Pineda
    // INRIA Grenoble Rhone-Alpes
    // E-mail: xavi.alameda@gmail.com
    // 
    // This is part of the gtde program.
    // 
    // gtde is free software: you can redistribute it and/or modify
    // it under the terms of the GNU General Public License as published by
    // the Free Software Foundation, either version 3 of the License, or
    // (at your option) any later version.
    // 
    // This program is distributed in the hope that it will be useful,
    // but WITHOUT ANY WARRANTY; without even the implied warranty of
    // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    // GNU General Public License for more details.
    // 
    // You should have received a copy of the GNU General Public License
    // along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    // Compute the function''s value
    k2 = 0;
    for r = mtlb_imp(0,q)
        k2 = mtlb_a(k2,((combinations(q,r) .*((-tau) .^mtlb_s(q,r))) .*mtlb_s(T .^mtlb_a(mtlb_a(p,r),1),tau .^mtlb_a(mtlb_a(p,r),1))) ./mtlb_a(mtlb_a(p,r),1));
    end;
    
    // If asked, compute its derivative
    if %nargout>1 then
        k2d = 0;
        // The general formula works for r < q
        for r = mtlb_imp(0,mtlb_s(q,1))
            k2d = mtlb_a(k2d,(combinations(q,r)*((-tau) .^mtlb_s(mtlb_s(q,r),1))) .*mtlb_a((mtlb_s(r,q)*mtlb_s(T .^mtlb_a(mtlb_a(r,p),1),tau .^mtlb_a(mtlb_a(r,p),1)))/mtlb_a(mtlb_a(r,p),1),tau .^mtlb_a(mtlb_a(r,p),1)));
      end;
            k2d = mtlb_s(k2d,tau .^mtlb_a(q,p));
        end;
    
        // If asked, compute its second derivative
        if %nargout>2 then
            %v0_1 = size(tau);  k2dd = zeros(%v0_1(1),%v0_1(2));
            // The general formula works for r < q-1
            for r = mtlb_imp(0,mtlb_s(q,2))
                k2dd = mtlb_a(k2dd,(combinations(q,r)*((-1) .^mtlb_s(mtlb_s(q,r),2))) .*mtlb_a(((mtlb_s(1+r,q)*mtlb_s(r,q))*mtlb_s(((T .^p) .*(tau .^mtlb_s(q,1))) .*((T ./tau) .^(r+1)),tau .^mtlb_s(mtlb_a(q,p),1)))/mtlb_a(mtlb_a(r,p),1),mtlb_s(mtlb_s(r,p),2*q)*(tau .^mtlb_s(mtlb_a(q,p),1))));
        end;
        // If both q and p are 0, do not add anything
        if mtlb_logic(mtlb_a(q,p),">",0) then
            // Term for r = q-1
            k2dd = mtlb_a(k2dd,(q*mtlb_a(mtlb_a(q,p),1))*(tau .^mtlb_s(mtlb_a(q,p),1)));
            // Term for r = q
            k2dd = mtlb_s(k2dd,mtlb_a(q,p)*(tau .^mtlb_s(mtlb_a(q,p),1)));
        end;
    end;
return;
endfunction
