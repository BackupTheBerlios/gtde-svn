function [k1, k1d, k1dd] = K1EUSIPCO(p,q,tau,T)

% [k1 k1d] = K1(p,q,tau,T)
% 
%     This function will compute the following formula
%     k1 = sum_{p=0}^q choose(q,k) * (T-tau)^(q-k) * tau^(p+k+1) / (p+k+1)
%     and its derivatives (1st and 2nd) wigh respect to tau. T and tau may 
%     be two vectors of the same size, or vector & number.
%
%     see also K2, CCInterpolation
    
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
%             % Particular case q = 1, p = 0
%             if q == 1 && p == 0
%                 k1dd = k1dd - 2*tau;
%             else
%                 k1dd = k1dd + q *tau.^(q+p-2).*( (q-1)*T - (q+1)*tau);
%             end
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
end