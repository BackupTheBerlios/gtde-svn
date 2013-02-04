function corr = DiscreteCrossCorrelation(signal1,signal2,maxLag)

%Discrete cross-correlation
%
% USAGE: Cross-Correlation = DiscreteCrossCorrelation(FirstSignal,SecondSignal,MaxLag)
%
% PARAMETERS:
%   FirstSignal ~ values of the first signal
%   SecondSignal ~ values of the second signal
%   MaxLag ~ the maximum value of the lag we want to compute
%
% RETURN VALUE:
%   Cross-Correlation ~ values of the cross-correlation signal
% 
% DESCRIPTION:
%   This function computes the discrete cross-corrlation of two 1D signals.
%

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

    %%% Check input
    if nargin < 3
        error('Usage: Cross-Correlation = DiscreteCrossCorrelation(FirstSignal,SecondSignal,MaxLag)');
    end
    
    %%% Initialize output
    corr = zeros(1,2*maxLag+1);
    % How many samples we use to compute each of the cross-correlation
    % values
    chunkSize = length(signal1)-maxLag;
    
    %%% Computation
    % Positive values of the delay
    for xcPos = 2:maxLag+1,
        nXCSamples = 0;
        for shift = 0:maxLag+1-xcPos,
            corr(maxLag+xcPos) = corr(maxLag+xcPos) + ...
                sum(signal1(xcPos+shift:chunkSize-1+xcPos+shift).*signal2(shift+1:chunkSize+shift));
            nXCSamples = nXCSamples + 1;
        end
        corr(maxLag+xcPos) = corr(maxLag+xcPos) / nXCSamples;
    end
    % Negative values of the delay
    for xcPos = 2:maxLag+1,
        nXCSamples = 0;
        for shift = 0:maxLag+1-xcPos,
            corr(maxLag+2-xcPos) = corr(maxLag+2-xcPos) + ...
                sum(signal2(xcPos+shift:chunkSize-1+xcPos+shift).*signal1(shift+1:chunkSize+shift));
            nXCSamples = nXCSamples + 1;
        end
        corr(maxLag+2-xcPos) = corr(maxLag+2-xcPos) / nXCSamples;
    end
    % At 0
    nXCSamples = 0;
    for shift = 0:maxLag,
        corr(maxLag+1) = corr(maxLag+1) + ...
            sum(signal1(shift+1:shift+chunkSize).*signal2(shift+1:shift+chunkSize));
        nXCSamples = nXCSamples + 1;
    end
    corr(maxLag+1) = corr(maxLag+1) / nXCSamples;

%     %%% Computation
%     % Positive values of the delay
%     for xcPos = 2:maxLag+1,
%         corr(maxLag+xcPos) = sum( signal1(1:end-xcPos+1) .* signal2(xcPos:end) );
%         corr(maxLag+2-xcPos) = sum( signal2(1:end-xcPos+1) .* signal1(xcPos:end) );
%     end
%     corr(maxLag+1) = sum(signal1.*signal2);
    
return

