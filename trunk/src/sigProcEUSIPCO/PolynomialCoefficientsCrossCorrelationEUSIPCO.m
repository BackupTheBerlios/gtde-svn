function PCCC = PolynomialCoefficientsCrossCorrelationEUSIPCO(firstCoefficients,secondCoefficients,MaxLag)

% PolynomialCoefficientsCrossCorrelation computes the cross-correlation of 
% the coefficients of the polynomial interpolation of two signals.
% 
%   PCCC = PolynomialCoefficientsCrossCorrelation(firstCoefficients,secondCoefficients)
%   outputs a (D1+1)x(D2+1)x(2*MaxLag+1) matrix, where D1 and D2 are the degrees of the 
%   polynomial interpolations and MaxLag is the size of the desired output correlation.
%   The default value for MaxLag is NSamples-1, where NSamples is the
%   length of firstCoefficients (or secondCoefficients, which must be the
%   same).
%   
%   see also InterpolateCrossCorrelation, xcorr

    %%% Input check
    if nargin < 2
        error('Usage: PCCC = PolynomialCoefficientsCrossCorrelation(firstCoefficients,secondCoefficients[,MaxLag])');
    end
    if size(firstCoefficients,1) ~= size(secondCoefficients,1)
        error('Polynomial coefficients should have the same size.');
    end
    
    %%% General variables
    % Compute the output correlation size
    signalSize = size(firstCoefficients,1);
    if nargin < 3
        MaxLag = signalSize-1;
    end

    % Polynomial degress
    M1 = size(firstCoefficients,2);
    M2 = size(secondCoefficients,2);

    % Declare output variable
    PCCC = zeros(M1,M2,2*MaxLag+1);

    %%% Computation
    for m1 = 1:M1,
        for m2 = 1:M2,
            PCCC(m1,m2,:) = xcorrEUSIPCO(firstCoefficients(:,m1),...
                                  secondCoefficients(:,m2),...
                                  MaxLag,...
                                  'none');
        end
    end

end