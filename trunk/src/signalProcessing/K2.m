function [k2 k2d k2dd] = K2(p,q,tau,T)

%K2   K2 function
%
% USAGE: [k2 k2d k2dd] = K2(p,q,tau,T)
%
% PARAMETERS:
%  p,q ~ Degree of the interpolation polynomial of the first and second
%       signals.
%  tau ~ Point (modulo T) in which the cross-correlation function is 
%       evaluated.
%  T ~ Sampling period of the discrete signals.
%
% RETURN VALUE:
%  k2,k2d,k2dd ~ Value of the function, the first and the second derivative
%       at tau respectively.
%
% DESCRIPTION:
%     The K2 function corresponds to the following formula
%     K2 = sum_{k=0}^q choose(q,k) * (tau)^(q-k) * ( T^(k+p+1) - tau^(k+p+1)) / (k+p+1)
%     and its derivatives (1st and 2nd) wigh respect to tau. T and tau may 
%     be two vectors of the same size, or vector and number.
%
% REFERENCES:
%     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
%     estimation-based sound source localisation (gTDESSL). Research Report 
%     RR-7988, INRIA, June 2012.
%
%     see also K1, CCInterpolation

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
    k2 = 0;
    for r = 0:q,
        k2 = k2 + nchoosek(q,r) .* (-tau).^(q-r) .* ( T.^(p+r+1)-tau.^(p+r+1) ) ./ (p+r+1);
    end
    
    % If asked, compute its derivative
    if nargout > 1
        k2d = 0;
        % The general formula works for r < q
        for r = 0:q-1,
            k2d = k2d + nchoosek(q,r) * (-tau).^(q-r-1) .* ( (r-q)*(T.^(r+p+1)-tau.^(r+p+1))/(r+p+1) + tau.^(r+p+1));
        end
        k2d = k2d - tau.^(q+p);
    end
    
    % If asked, compute its second derivative
    if nargout > 2
        k2dd = zeros(size(tau));
        % The general formula works for r < q-1
        for r = 0:q-2,
            k2dd = k2dd + nchoosek(q,r) * (-1).^(q-r-2).*( (1+r-q)*(r-q)*(T.^p.*tau.^(q-1).*(T./tau).^(r+1)-tau.^(q+p-1))/(r+p+1) + (r-p-2*q)*(tau).^(q+p-1));
        end
        % If both q and p are 0, do not add anything
        if q+p > 0
            % Term for r = q-1
            k2dd = k2dd + q*(q+p+1)*tau.^(q+p-1);
            % Term for r = q
            k2dd = k2dd - (q+p)*tau.^(q+p-1);
        end
    end
return

