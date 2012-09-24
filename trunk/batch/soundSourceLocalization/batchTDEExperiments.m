function batchTDEExperiments(experimentOptions)

%Batch the experiments on geometrically constrained time delay estimation
%
% USAGE: batchTDEExperiments(experimentOptions)
%
% PARAMETERS:
%   experimentOptions ~ structure containing several structures describing 
%     the options of the experiment. The sub-structures are described here.
%
%     rootFolder ~ string: the root folder for data and results.
%     dimension ~ dimension of the experiment (2 or 3).
%     ism ~ boolean indicating the use of the ISM model.
%     sourcePositionOptions ~ options for the generation of the sound source 
%         positions. Please see SpecifySourcePositionOptions.
%     microphonePositionOptions ~ options to specify the microphone positions.
%         Please see SpecifyMicrophonesOptions.
%     ismOptions ~ options to specify the ISM filters. Please see SpecifyISMFilterOptions.
%     dataOptions ~ options to specify the used data. Please see SpecifyDataOptions.
%     methodOptions ~ options to specify the used method. Please see SpecifyMethodOptions.
% 
%   see also ProcessExperimentOptions, SpecifySourcePositionOptions, SpecifyMicrophonesOptions
%   SpecifyISMFilterOptions, SpecifyDataOptions, SpecifyMethodOptions

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

    %%% General variables
    % Check options
    experimentOptions = ProcessExperimentOptions(experimentOptions);
    % Get Source positions
    experimentOptions = GetSourcePositions(experimentOptions);
    % Get microphones positions
    experimentOptions = GetMicrophonePositions(experimentOptions);
    % Get ISM Data
    experimentOptions = GetISMData(experimentOptions);
    % Generate Initialization positions
    experimentOptions = GetInitializationPositions(experimentOptions);
    % Get signals
    experimentOptions = GetSignals(experimentOptions);
    
    %%% Useful variables
    % Number of sound source positions
    nSourcePositions = size(experimentOptions.sourcePositions,1);
    % Number of SNR values
    nSNRValues = length(experimentOptions.dataOptions.snrValues);
    % Number of T60 values
    nT60 = length(experimentOptions.ismOptions.t60);
    % Number of signals
    nSignals = numel(experimentOptions.signals);
    % Store the results
    foundTDEs = cell(nSignals,...
                     nSourcePositions,...
                     nSNRValues,...
                     nT60);
    % Control the time
    donePart = 0;
    tStart = clock;
    
    % Compute the ISM data if needed.
    if experimentOptions.ismOptions.generate
        fprintf('Computation of the ISM models...\n');
        experimentOptions = GenerateISMData(experimentOptions);
        fprintf('ISM Computation done.\n');
    end

    
    %%% Looping over all the parameters
    % For each position
    for sPosition = 1:1%nSourcePositions,
        fprintf('Position %d/%d:\n',sPosition,nSourcePositions);
        % For each signal
        for sSignal = 1:nSignals
            fprintf('    Signal %d/%d   <',sSignal,nSignals);
            % For each SNR value
            for sSNR = 1:nSNRValues,
                % For each T60 value,
                for sT60 = 1:nT60,
                    
                    % Get the signals
                    [signals fs] = GetSignalRealization(experimentOptions,sSignal,sPosition,sSNR,sT60);
                    % Reassign the sampling frequency
                    experimentOptions.samplingFrequency = fs;
                    % Allocate results
                    foundTDEs{sSignal,sPosition,sSNR,sT60} = zeros(numel(signals),experimentOptions.dimension);
                    
                    % Loop for each partial signal
                    for subs = 1:numel(signals),
                        
                        % Get partial signal
                        signalPiece = signals{subs};
                        % Apply the method
                        foundTDE = EstimateTimeDelays(signalPiece,experimentOptions);
                        % Save the result
                        foundTDEs{sSignal,sPosition,sSNR,sT60}(subs,:) = foundTDE; 
                    end
                    
                    % Notification counter
                    donePart = donePart + 100/(nSignals*nSourcePositions*nSNRValues*nT60);
                    fprintf('.');
                    
                end % sT60
                fprintf('-');
            end % sSNR
            fprintf('|');
            
            % Notify
            tEnd = clock;
            eTime = etime(tEnd,tStart);
            rTime = eTime*(100-donePart)/donePart;
            fprintf('>   [Done by: %s]\n',datestr(now+(rTime/86400),0));
            
        end % sPosition
    end % sSignal
    
    % Experiment number
    directoryStructure = dir(experimentOptions.resultOptions.folder);
    expNum = numel(directoryStructure) - 2;
    % End time
    experimentOptions.endingTime = datestr(clock,'yyyy-mm-dd-HH-MM-SS');
    experimentOptions.expNum = expNum;
    % Set the file name and save
    fileName = strcat(experimentOptions.resultOptions.folder,...
                      num2str(expNum),'.mat');
    save(fileName,'experimentOptions','foundTDEs');
    
end % function