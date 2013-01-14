function [k1 k1d k1dd] = K1(p,q,tau,T)

%K1   K1 function
%
% USAGE: [k1 k1d k1dd] = K1(p,q,tau,T)
%
% PARAMETERS:
%  p,q ~ Degree of the interpolation polynomial of the first and second
%       signals.
%  tau ~ Point (modulo T) in which the cross-correlation function is 
%       evaluated.
%  T ~ Sampling period of the discrete signals.
%
% RETURN VALUE:
%  k1,k1d,k1dd ~ Value of the function, the first and the second derivative
%       at tau respectively.
% 
% DESCRIPTION:
%     The K1 function corresponds to the following formula
%     k1 = sum_{p=0}^q choose(q,k) * (T-tau)^(q-k) * tau^(p+k+1) / (p+k+1)
%     and its derivatives (1st and 2nd) wigh respect to tau. T and tau may 
%     be two vectors of the same size, or vector and number.
%
% REFERENCES:
%     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
%     estimation-based sound source localisation (gTDESSL). Research Report 
%     RR-7988, INRIA, June 2012.
%
%     see also K2, CCInterpolation

% Copyright 2012, Xavier Alameda-Pineda
% INRIA Grenoble Rhone-Alpes
% E-mail: xavi.alameda@gmail.com
% 
% This is part of the gtde program.
% 
% gtde is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

    
    % Compute the function's value
    k1 = 0;
    for r = 0:q,
        k1 = k1 + nchoosek(q,r) * (T-tau).^(q-r) .* tau.^(p+r+1) / (p+r+1);
    end
    
    % If asked, compute its derivative
    if nargout > 1
        k1d = 0;
        for r = 0:q-1,
            k1d = k1d + nchoosek(q,r) * (T-tau).^(q-r-1) .* (tau).^(r+p) .* (T - (1+p+q)*tau/(r+p+1));
        end
        % r = q may give some numerical problems, special formula
        k1d = k1d + tau.^(q+p);
    end
    
    % If asked, compute its second derivative
    if nargout > 2
        k1dd = zeros(size(tau));
        % -----
        % r = 0
        % -----
        if p >= 1
            if q >= 2
                % Same formula
                r = 0;
                k1dd = k1dd + nchoosek(q,r) * (T-tau).^(q-r-2) .* (tau).^(r+p-1) .*( (r+p)*T.^2 - 2*(q+p)*T.*tau + (q+p)*(q+p+1)*tau.^2/(r+p+1));
            elseif q == 1
                k1dd = k1dd - tau.^(p-1).*( p*(tau-T) + 2*tau);
            else
                k1dd = k1dd + p*tau.^(p-1);
            end
        else
            % Special cases
            if q >= 2
                k1dd = k1dd + q*(T-tau).^(q-2).*( (q+1)*tau - 2*T );
            elseif q == 1
                k1dd = k1dd - 2;
            end
        end
        % -------
        % r = q-1, exists if q >= 1, but for q = 1, the r = q will take
        % care of it
        % -------
        if q >= 2
            k1dd = k1dd + q*tau.^(q+p-2).*( (q+p)*(T-tau) - (T+tau) );
        end
        % -----
        % r = q, it exists for q >= 1, since for q = 0, r = q will take
        % care of it
        % -----
        if q >= 1
            k1dd = k1dd + (q+p)*tau.^(q+p-1);
        end
        % Intermediate terms
        for r = 1:q-2,
            k1dd = k1dd + nchoosek(q,r) * (T-tau).^(q-r-2) .* (tau).^(r+p-1) .*( (r+p)*T.^2 - 2*(q+p)*T.*tau + (q+p)*(q+p+1)*tau.^2/(r+p+1));
        end
    end
return
