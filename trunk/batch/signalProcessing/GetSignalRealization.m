function [signals, fs] = GetSignalRealization(experimentOptions,sSignal,sPosition,sSNR,t60)

%Get the realization of a signal
%
% [signals, fs] = GetSignalRealization(experimentOptions,sSignal,sPosition,sSNR,sT60)
%
% PARAMETERS:
%   experimentOptions ~ options of the experiment, see batchTDEExperiments
%   sSignal ~ index of the signal
%   sPosition ~ index of the sound source positions
%   sSNR ~ index of the SNR value
%   sT60 ~ index of the t60 value
%
% RETURN VALUE:
%   signals ~ signals received at the microphone(s)
%   fs ~ sampling frequency of the signals
%
% DESCRIPTION: Adds gaussian white noise to a signal
% 
%   see also batchTDEExperiments

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
%

% Depending on the data
    % In case of synthetic data, they are functions
    switch experimentOptions.dataOptions.type
        case 'synthetic'
            fs = experimentOptions.dataOptions.samplingFrequency;
            % ISM case
            if experimentOptions.ism
                % Generate signal
                signal = GenerateDiscreteSignals(1./experimentOptions.dataOptions.samplingFrequency,...
                                                 experimentOptions.dataOptions.cutLength,...
                                                 experimentOptions.signals{sSignal});
                % Sampling frequency index
                sT60 = find(experimentOptions.ismOptions.t60 == t60);
                if isempty(sT60)
                    error(['ISM filters were not computed for a T60 of ' num2str(experimentOptions.dataOptions.t60(sT60)) '.']);
                end
                sFS = find(experimentOptions.ismOptions.samplingFrequencies == experimentOptions.dataOptions.samplingFrequency);
                if isempty(sFS)
                    error(['ISM filters were not computed for sampling frequency of ' num2str(experimentOptions.dataOptions.samplingFrequency) '.']);
                end
                % ISM File Name
                ISMFileName = GetISMFileName(experimentOptions,sT60,sPosition,sFS);
                if ~(exist(ISMFileName,'file') == 2)
                    error(strcat('ISM filter does not exist (',ISMFileName,').'));
                end
                % Apply RIR filter
                partialSignals = ISM_AudioData(ISMFileName,signal,'SilentFlag',1);
                partialSignals = partialSignals'; 
            % Non ISM case
            else
                % Generate the signals
                partialSignals = GenerateDiscreteSignals(...
                    experimentOptions.sourcePositions(sPosition,:),...
                    experimentOptions.microphonesPositions,...
                    1./experimentOptions.dataOptions.samplingFrequency,...
                    experimentOptions.dataOptions.cutLength,...
                    experimentOptions.signals{sSignal}...
                    );
            end
        case 'simulated'
            %%%%%%%%%%%%% WARNING %%%%%%%%%%%%%%%
            % THE NON REVERBERANT CASE IS T60 = 0
            %%%%%%%%%%%%% WARNING %%%%%%%%%%%%%%%
            % Checking the installation of ISM
            if isempty(which('ISM_AudioData'))
                fprintf('The software to compute ISM is not installed.\n');
                fprintf('Please go to http://www.eric-lehmann.com/ism code.html to download and install.\n');
                error('Software missing.');
            end
            % Load the signal
            if sPosition > numel(experimentOptions.signals)
                error('sSignal should be a valid signal index.');
            end
            fprintf('%s',experimentOptions.signals{sPosition});
            [signal, fs] = wavread(experimentOptions.signals{sPosition});
            % Sampling frequency index
            sFS = find(experimentOptions.ismOptions.samplingFrequencies == fs);
            if isempty(sFS)
                error(['ISM filters were not computed for sampling frequency of ' num2str(fs) '.']);
            end
            % Sampling frequency index
            sT60 = find(experimentOptions.ismOptions.t60 == t60);
            if isempty(sT60)
                error(['ISM filters were not computed for a T60 of ' num2str(experimentOptions.dataOptions.t60(sT60)) '.']);
            end
            % File name
            ISMFileName = GetISMFileName(experimentOptions,sT60,sPosition,sFS);
            % Check the for the file
            if ~(exist(ISMFileName,'file') == 2)
                error(strcat('ISM filter does not exist (',ISMFileName,').'));
            end
            % Apply RIR filter
            partialSignals = ISM_AudioData(ISMFileName,signal,'SilentFlag',1);
            partialSignals = partialSignals';
            if exist(experimentOptions.microphonePositionOptions,'subSet'),
                partialSignals = partialSignals(experimentOptions.microphonePositionOptions.subSet,:);
            end
        case 'real'
            % Compute the recording index, from sPos
    %         numTilts = numel(experimentOptions.tiltIndices);
    %         numPans = numel(experimentOptions.panIndices);
            sPos = sPosition-1;
            tiltIndex = 3*floor(sPos/6)+1;
            panIndex = 6*(sPos-6*fix(sPos/6))+1;
            wavIndex = 37*(tiltIndex-1)+panIndex;
    %         fprintf('%s\n',strcat(experimentOptions.signals,num2str(wavIndex),'.wav'));
            [partialSignals, fs] = wavread(strcat(experimentOptions.signals,num2str(wavIndex),'.wav'));
            partialSignals = partialSignals';
    end
    
    %%% Store the signals for iorgos
    % Shift'em
    if exist(experimentOptions.microphonePositionOptions,'subSet'),
        tdeSamples = TDEGeometricDirect(experimentOptions.sourcePositions(sPosition,:),experimentOptions.microphonePositions(experimentOptions.microphonePositionOptions.subSet,:));
    else
        tdeSamples = TDEGeometricDirect(experimentOptions.sourcePositions(sPosition,:),experimentOptions.microphonePositions);
    end
    tdeSamples = round(tdeSamples*fs);
    shiftedSignals = zeros(size(partialSignals));
    shiftedSignals(1,:) = partialSignals(1,:);
    for mm=2:size(experimentOptions.microphonePositions,1),
%         fprintf('%d ',tdeSamples(mm-1));
        if tdeSamples(mm-1) < 0
            shiftedSignals(mm,-tdeSamples(mm-1):end) = partialSignals(mm,1:end+tdeSamples(mm-1)+1);
        elseif tdeSamples(mm-1) > 0
            shiftedSignals(mm,1:end-tdeSamples(mm-1)+1) = partialSignals(mm,tdeSamples(mm-1):end);
        else
            shiftedSignals(mm,:) = partialSignals(mm,:);
        end
    end
    partialSignals = shiftedSignals;
       
    %%% Cut the signals
    pieceLength = experimentOptions.dataOptions.cutLength*fs; 
    NPieces = floor(size(partialSignals,2)/pieceLength);
    signals = cell(NPieces,1);

    for np = 1:NPieces,
        signals{np} = partialSignals( :, ((np-1)*pieceLength+1):(np*pieceLength) );
        if ~strcmp(experimentOptions.dataOptions.type,'real'),
            %%% Add noise
            for channel = 1:size(partialSignals,1),
                signals{np}(channel,:) = AddWhiteNoise(signals{np}(channel,:),experimentOptions.dataOptions.snrValues(sSNR));
            end
        end
    end
    
%     %%% Save'em
%     save(strcat('/local_scratch/alamedap/dataGeorge/UpToSubsample_',num2str(sSignal),'_',num2str(sPosition),'_',num2str(sSNR),'_',num2str(sT60),'.mat'),'signals');
end