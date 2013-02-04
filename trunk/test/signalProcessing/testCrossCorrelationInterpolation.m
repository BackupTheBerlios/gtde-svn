function testCrossCorrelationInterpolation

%test the computation of the polynomial coefficients

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

    clc; close all;

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
    myFun = @(x) SyntheticSignal(x,2,[20000,200]);

    % Sampling frequency
    samplingPeriod = 1/48000;

    % Length
    Length = 0.1;

    % Generate the signals
    [rawSignals] = GenerateDiscreteSignals(Positions(truePosition,:),...
                                      MICS,...
                                      samplingPeriod,...
                                      Length,...
                                      myFun);
%     signals = signals + 0.2*(rand(size(signals))-0.5);
%     windowLenght = round(size(rawSignals,2)/5);
%     myW = blackman(windowLenght);
%     for ii=1:size(rawSignals,1),
%         signals(ii,:) = conv(myW,rawSignals(ii,:));
%     end
    signals = rawSignals;


    %%% General variables
    % Number of microphones
    NMics = size(MICS,1);
    % Number of samples
    NSamples = size(signals,2);
    % Maximum TDE values from the microphones
    maxTDE = TDEmax(MICS);
    maxTDESamples = ceil(1.1*maxTDE/samplingPeriod);
    maxLAG = max(maxTDESamples);
    
    %%% Compute the interpolation coefficients
    PCCC = Signals2PCCC(signals,MICS,1/samplingPeriod);
    
    %%% tau and step
    tau = -5*samplingPeriod:0.01*samplingPeriod:5*samplingPeriod;
%     tau = tau + 0.0001*samplingPeriod*(rand(size(tau))-0.5);
%     tau = -10.^(-50:0)*samplingPeriod;
    step = 1e-6;
    
    for mic1=1:NMics,
        for mic2=mic1+1:NMics
            % Derivatives
%             showName = ['Derivative MICS=(' num2str(mic1) ',' num2str(mic2) ')'];
%             myFun = @(x) CrossCorrelationInterpolation(PCCC{mic1,mic2},x,samplingPeriod);
%             checkDerivatives(myFun,tau,step,showName);
%             % Seconds
%             showName = ['Second derivative MICS=(' num2str(mic1) ',' num2str(mic2) ')'];
%             myFun = @(x) HessCrossCorrelation(PCCC{mic1,mic2},x,samplingPeriod);
%             checkDerivatives(myFun,tau,step,showName);
            % Plot
            [cc, ccd, ccdd] = CrossCorrelationInterpolation(PCCC{mic1,mic2},tau,samplingPeriod);
%             close all;
%             figure;
%             for ii=1:4,
%                 for jj=1:4,
%                     subplot(4,4,4*(ii-1)+jj);
%                     plot(squeeze(PCCC{mic1,mic2}(ii,jj,:)),'bo-');
%                 end
%             end
%             figure;
%             subplot(2,1,1), plot(signals(mic1,:));
%             subplot(2,1,2), plot(signals(mic2,:));
            figure(1);
            plot(tau,cc,'ob');
            figure(2);
            plot(tau,ccd,'gs');
            figure(3);
            plot(tau,ccdd,'rd');
        end
    end
end

function [g,h] = HessCrossCorrelation(PCCC,x,samplingPeriod)
    [~,g,h] = CrossCorrelationInterpolation(PCCC,x,samplingPeriod);
end

function checkDerivatives(fun, X0, step, name)
    % Check out the dimension
    Dimension = size(X0,1);
    % Compute the values provided by the function
    [~, Gradient] = fun(X0);
    % Check the gradient
    NGradient = zeros(size(Gradient));
    for ii=1:Dimension
        % Direction
        Direction = zeros(size(Gradient));
        Direction(ii,:) = 1;
        % Forward point
        XF = X0 + step*Direction;
        % Backward point
        XB = X0 - step*Direction;
        % Compute the gradient
        NGradient(ii,:) = (fun(XF)-fun(XB));
    end
    % Normalize
    NGradient = NGradient / (2*step);
    % Printout
    fprintf('========\n');
    fprintf('%s\n',name);
    fprintf('========\n');
    for ii=1:Dimension,
        AbsDiff = abs(Gradient(ii,:)-NGradient(ii,:));
        RelDiff = zeros(size(AbsDiff));
        RelDiff(Gradient(ii,:)~=0) = AbsDiff(Gradient(ii,:)~=0)./abs(Gradient(ii,Gradient(ii,:)~=0));
        fprintf('Abs ~ (%1.5e,%1.5e) \t Rel: (%d;%1.5e)\n',min(AbsDiff),max(AbsDiff),sum(Gradient(ii,:)~=0),max(RelDiff));
        figure;
        plot(X0,AbsDiff,'kx');
    end
end