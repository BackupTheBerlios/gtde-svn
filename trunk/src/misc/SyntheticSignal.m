function y = SyntheticSignal(x,choice,params)

    switch choice
        case 1
            y = sin(params(1)*x);
        case 2
            y = zeros(size(x));
            y(x>=0) = exp(-abs(params(2)*x(x>=0))).*sin(params(1)*x(x>=0));
%             y = exp(-abs(params(2)*x)).*sin(params(1)*x);
        case 3
            y = zeros(size(x));
            y(x>=0) = exp(-abs(params(1)*x(x>=0)));
%             y = exp(-abs(params(1)*x));
        otherwise
            error('Not implemented');
    end

end