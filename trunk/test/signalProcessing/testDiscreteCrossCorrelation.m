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

    close all;

    %%% Generate DAta
    
    % Scale parameters
    a = 0.1;
    c = a / sqrt(2);
    % Microphones in a regular tetrahedron
    MICS = [ a   0   -c ;...
            -a   0   -c ;...
             0    a   c ;...
             0   -a   c];
    % Random rotation
%     Vx = rand(3,1)-0.5;
%     Vy = rand(3,1)-0.5;
%     alpha = 2*pi*rand(1);
%     r = rotation(Vx,Vy,alpha);
%     MICS = MICS * r;

    % Generate potision2
    positionOptions = struct('coordinateSystem','cartesian');
    % Declare the bounds in the three variables
    positionOptions.bounds = cell(1,3);
    % Set the bounds for the radius, the azimuth and elevation will de the
    % defaults
    positionOptions.bounds{1} = [-2 2];
    positionOptions.bounds{2} = [-2 2];
    positionOptions.bounds{3} = [-2 2];
    % Number of intervals
    positionOptions.numberOfIntervals = [8, 8, 8];
    % Do not use the evelation's extrema
    positionOptions.useExtrema = true(1,3);
    positionOptions.useExtrema(3) = false;
    % Generate positions
    global Positions;
    Positions = GeneratePositions(3, positionOptions);
    Positions(284,:) = [];
    

    % Choose the position
    truePosition = 56;

    % Sinusoidal
    % 120Hz
    % myFun = @(x) SyntheticSignal(x,1,120);
    % Sinusoidal with exponential decay
    % 120Hz and decay constant = 3
    % myFun = @(x) SyntheticSignal(x,2,[120,3]);
    % Exponential decay
    myFun = @(x) SyntheticSignal(x,2,[8000,1]);

    % Sampling frequency
    samplingPeriod = 1/48000;

    % Length
    Length = 0.1;

    % Generate the signals
    [signals] = GenerateDiscreteSignals(Positions(truePosition,:),...
                                      MICS,...
                                      samplingPeriod,...
                                      Length,...
                                      myFun);
%     signals = signals + 0.2*(rand(size(signals))-0.5);

    maxLag = 40;
%     times = 0:0.01:1;
%     signal1 = sin(10*times);
%     signal2 = sin(-30*times);

    signal1 = signals(1,:);
    signal2 = signals(2,:);
    
    corr = DiscreteCrossCorrelation(signal1,signal2,maxLag);
    
    figure
    plot(signal1);
    figure
    plot(signal2);
    figure
    plot(-maxLag:maxLag,corr);
    
return
