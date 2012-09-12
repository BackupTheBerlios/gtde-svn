function [signals fs] = GetSignalRealization(experimentOptions,sSignal,sPosition,sSNR,sT60)

%Get the realization of a signal
%
% [signals fs] = GetSignalRealization(experimentOptions,sSignal,sPosition,sSNR,sT60)
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
            fs = experimentOptions.samplingFrequency;
            % ISM case
            if experimentOptions.ism
                % Generate signal
                signal = GenerateDiscreteSignals(1./experimentOptions.samplingFrequency,...
                                                 experimentOptions.length,...
                                                 experimentOptions.signals{sSignal});
                % ISM File Name
                ISMFileName = getISMFileName(experimentOptions,...
                    experimentOptions.T60(sT60),...
                    experimentOptions.sourcePositions(sPosition,:),...
                    fs);
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
                    1./experimentOptions.samplingFrequency,...
                    experimentOptions.length,...
                    experimentOptions.signals{sSignal}...
                    );
            end
        case 'simulated'
            %%%%%%%%%%%%% WARNING %%%%%%%%%%%%%%%
            % THE NON REVERBERANT CASE IS T60 = 0
            %%%%%%%%%%%%% WARNING %%%%%%%%%%%%%%%
            % Load the signal
            switch sSignal,
                case 1,
                    [signal fs] = wavread(experimentOptions.signals{2});
                case 2,
                    [signal fs] = wavread(experimentOptions.signals{9});
                case 3,
                    index = randi(numel(experimentOptions.signals),1);
                    [signal fs] = wavread(experimentOptions.signals{index});
                otherwise
                    error('Signal should be either 1 (caugh), 2 (female speech) or 3 (random sound).');
            end
            % ISM File Name
            sFS = find(experimentOptions.ismOptions.samplingFrequencies == fs);
            if isempty(sFS)
                error(['ISM filters were not computed for sampling frequency of ' num2str(fs) '.']);
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
        case 'real'
            % Compute the recording index, from sPos
    %         numTilts = numel(experimentOptions.tiltIndices);
    %         numPans = numel(experimentOptions.panIndices);
            sPos = sPosition-1;
            tiltIndex = 3*floor(sPos/6)+1;
            panIndex = 6*(sPos-6*fix(sPos/6))+1;
            wavIndex = 37*(tiltIndex-1)+panIndex;
    %         fprintf('%s\n',strcat(experimentOptions.signals,num2str(wavIndex),'.wav'));
            [partialSignals fs] = wavread(strcat(experimentOptions.signals,num2str(wavIndex),'.wav'));
            partialSignals = partialSignals';
    end
        
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
end