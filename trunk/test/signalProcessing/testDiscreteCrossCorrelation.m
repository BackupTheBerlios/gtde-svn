function testDiscreteCrossCorrelation

%test the DiscreteCrossCorrelationFunction

% Copyright 2012, Xavier Alameda-Pineda
% INRIA Grenoble Rhône-Alpes
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

%     % Scale parameters
%     a = 0.1;
%     c = a / sqrt(2);
%     % Microphones in a regular tetrahedron
%     MICS = [ a   0   -c ;...
%             -a   0   -c ;...
%              0    a   c ;...
%              0   -a   c];
% 
%     % Frequency and decaying factor
%     F = 800;
%     T = 10;
% 
%     myFun = @(x) SyntheticSignal(x,2,[F,T]);
% 
%     % Sampling frequency
%     samplingPeriod = 1/48000;
% 
%     % Index of the used positions
%     truePosition = 1;
% 
%     % Length
%     Length = 0.1;
% 
%     % Generate the signals
%     [signals samplingTimes] = GenerateDiscreteSignals(Positions(truePosition,:),...
%                                   MICS,...
%                                   samplingPeriod,...
%                                   Length,...
%                                   myFun);
    maxLag = 20;
    times = 0:0.01:1;
    signal1 = sin(10*times);
    signal2 = sin(-30*times);
    
    corr = DiscreteCrossCorrelation(signal1,signal2,maxLag);
    
    figure
    plot(times,signal1);
    figure
    plot(times,signal2);
    figure
    plot(-maxLag:maxLag,corr);
    
end