function [k2, k2d, k2dd] = K2EUSIPCO(p,q,tau,T)

% [k2 k2d] = K2(p,q,tau,T)
% 
%     This function will compute the following formula
%     sum_{k=0}^q choose(q,k) * (tau)^(q-k) * ( T^(k+p+1) - tau^(k+p+1)) / (k+p+1)
%     and its first and second derivatives wigh respect to tau. T and tau
%     maybe two vectors of the same size, of a vector and a number.
%
%     see also K1, CCInterpolation

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
end