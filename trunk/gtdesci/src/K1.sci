function [k1,k1d,k1dd] = K1(p,q,tau,T)

    // Output variables initialisation (not found in input variables)
    k1=[];
    k1d=[];
    k1dd=[];
    
    // Number of arguments in function call
    [%nargout,%nargin] = argn(0)
    
    // Display mode
    mode(0);
    
    // Display warning for floating point exception
    ieee(1);
    
    
    //K1   K1 function
    // 
    // USAGE: [k1 k1d k1dd] = K1(p,q,tau,T)
    // 
    // PARAMETERS:
    //  p,q ~ Degree of the interpolation polynomial of the first and second
    //       signals.
    //  tau ~ Point (modulo T) in which the cross-correlation function is 
    //       evaluated.
    //  T ~ Sampling period of the discrete signals.
    // 
    // RETURN VALUE:
    //  k1,k1d,k1dd ~ Value of the function, the first and the second derivative
    //       at tau respectively.
    // 
    // DESCRIPTION:
    //     The K1 function corresponds to the following formula
    //     k1 = sum_(p=0).entries^q choose(q,k) * (T-tau)^(q-k) * tau^(p+k+1) / (p+k+1)
    //     and its derivatives (1st and 2nd) wigh respect to tau. T and tau may 
    //     be two vectors of the same size, or vector and number.
    // 
    // REFERENCES:
    //     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
    //     estimation-based sound source localisation (gTDESSL). Research Report 
    //     RR-7988, INRIA, June 2012.
    // 
    //     see also K2, CCInterpolation
    
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
    k1 = 0;
    for r = mtlb_imp(0,q)
        k1 = mtlb_a(k1,((combinations(q,r)*(mtlb_s(T,tau) .^mtlb_s(q,r))) .*(tau .^mtlb_a(mtlb_a(p,r),1)))/mtlb_a(mtlb_a(p,r),1));
    end;
    
    // If asked, compute its derivative
    if %nargout>1 then
        k1d = 0;
        for r = mtlb_imp(0,mtlb_s(q,1))
            k1d = mtlb_a(k1d,((combinations(q,r)*(mtlb_s(T,tau) .^mtlb_s(mtlb_s(q,r),1))) .*(tau .^mtlb_a(r,p))) .*mtlb_s(T,(mtlb_a(mtlb_a(1,p),q)*tau)/mtlb_a(mtlb_a(r,p),1)));
        end;
        // r = q may give some numerical problems, special formula
        k1d = mtlb_a(k1d,tau .^mtlb_a(q,p));
    end;
    
    // If asked, compute its second derivative
    if %nargout>2 then
         %v0_1 = size(tau);  k1dd = zeros(%v0_1(1),%v0_1(2));
        // -----
        // r = 0
        // -----
        if mtlb_logic(p,">=",1) then
            if mtlb_logic(q,">=",2) then
                // Same formula
                r = 0;
                k1dd = mtlb_a(k1dd,((combinations(q,r)*(mtlb_s(T,tau) .^mtlb_s(mtlb_s(q,r),2))) .*(tau .^mtlb_s(mtlb_a(r,p),1))) .*mtlb_a(mtlb_s(mtlb_a(r,p)*(T .^2),((2*mtlb_a(q,p))*T) .*tau),((mtlb_a(q,p)*mtlb_a(mtlb_a(q,p),1))*(tau .^2))/mtlb_a(mtlb_a(r,p),1)));
            elseif mtlb_logic(q,"==",1) then
                k1dd = mtlb_s(k1dd,(tau .^mtlb_s(p,1)) .*mtlb_a(p*mtlb_s(tau,T),2*tau));
            else
                k1dd = mtlb_a(k1dd,p*(tau .^mtlb_s(p,1)));
            end;
        else
            // Special cases
            if mtlb_logic(q,">=",2) then
                k1dd = mtlb_a(k1dd,(q*(mtlb_s(T,tau) .^mtlb_s(q,2))) .*mtlb_s(mtlb_a(q,1)*tau,2*T));
            elseif mtlb_logic(q,"==",1) then
                k1dd = mtlb_s(k1dd,2);
            end;
        end;
        // -------
        // r = q-1, exists if q >= 1, but for q = 1, the r = q will take
        // care of it
        // -------
        if mtlb_logic(q,">=",2) then
            k1dd = mtlb_a(k1dd,(q*(tau .^mtlb_s(mtlb_a(q,p),2))) .*mtlb_s(mtlb_a(q,p)*mtlb_s(T,tau),mtlb_a(T,tau)));
        end;
        // -----
        // r = q, it exists for q >= 1, since for q = 0, r = q will take
        // care of it
        // -----
        if mtlb_logic(q,">=",1) then
            k1dd = mtlb_a(k1dd,mtlb_a(q,p)*(tau .^mtlb_s(mtlb_a(q,p),1)));
        end;
        // Intermediate terms
        for r = mtlb_imp(1,mtlb_s(q,2))
            k1dd = mtlb_a(k1dd,((combinations(q,r)*(mtlb_s(T,tau) .^mtlb_s(mtlb_s(q,r),2))) .*(tau .^mtlb_s(mtlb_a(r,p),1))) .*mtlb_a(mtlb_s(mtlb_a(r,p)*(T .^2),((2*mtlb_a(q,p))*T) .*tau),((mtlb_a(q,p)*mtlb_a(mtlb_a(q,p),1))*(tau .^2))/mtlb_a(mtlb_a(r,p),1)));
        end;
    end;
return;
endfunction
