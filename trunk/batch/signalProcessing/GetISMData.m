function experimentOptions = GetISMData(experimentOptions)

%Generate and save the data of the ISM Models
%
% USAGE: experimentOptions = GenerateISMData(exprimentOptions)
%
% PARAMETERS:
%   experimentOptions ~ structure containing the options of the experiment

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

    % If we use real data, just ignore it
    if ~strcmp(experimentOptions.dataOptions.type,'real') && experimentOptions.ismOptions.generate,
        % Checking the installation of ISM
        if isempty(which('ISM_RIR_bank'))
            fprintf('The software to compute ISM is not installed.\n');
            fprintf('Please go to http://www.eric-lehmann.com/ism code.html to download and install.\n');
            error('Software missing.');
        end
        % Check microphones are in the room
        checkMicrophonesRoom(experimentOptions);
        % If you want to generate, create the folder and generate data
        [status] = system(['mkdir -p ' experimentOptions.ismOptions.folder]);
        if status
            error(['An error ocurred while creating the ismFolder: ', experimentOptions.ismOptions.ismFolder]);
        end
        % Sampling Frequencies
        samplingFrequencies = experimentOptions.ismOptions.samplingFrequencies;
        donePart = 0;
        % Done part done
        % For each value of T60
        for nt = 1:numel(experimentOptions.ismOptions.t60)
            t60 = experimentOptions.ismOptions.t60(nt);
            % For each sampling frequency
            for nSF = 1:numel(samplingFrequencies),
            SF = samplingFrequencies(nSF);
                % For each sound source positions
                for np = 1:size(experimentOptions.sourcePositions,1)
                    % Extract position
                    position = experimentOptions.sourcePositions(np,:);
                    % ISM file name
                    ISMFileName = GetISMFileName(experimentOptions,nt,np,nSF);
                    % If the file does not exist, create it
                    if ~(exist(ISMFileName,'file') == 2)
                        ISMStruct = gTDE_ISM_setup();
                        % Modify it
                        ISMStruct.Fs = SF;
                        ISMStruct.room = experimentOptions.ismOptions.room;
                        ISMStruct.mic_pos = experimentOptions.microphonePositions;
                        ISMStruct.src_traj = position;
                        ISMStruct.T60 = t60;
                        ISMStruct.abs_weights = experimentOptions.ismOptions.absorptionWeights;
                        % Create the RIR filter
                        ISM_RIR_bank(ISMStruct,ISMFileName);
                    end
                    donePart = donePart + 100/(numel(experimentOptions.ismOptions.t60)*size(experimentOptions.sourcePositions,1)*numel(samplingFrequencies));
                    fprintf('T60: %1.1f %d [x= %1.3f, y= %1.3f, z= %1.3f] (%1.1f %%)\n',t60,SF,position(1),position(2),position(3),donePart);
                end %for np
            end % for nSF
        end % for nt
        % Save experiment options
        save(strcat(experimentOptions.ismOptions.folder,'experimentOptions.mat'),'experimentOptions');
    end
end

function checkMicrophonesRoom(experimentOptions)
    % For each microphone
    for mic = 1:size(experimentOptions.microphonePositionOptions,1)
        % For each dimension
        for d = 1:experimentOptions.dimension,
            % Check it to be positive and less than the room sizes
            if experimentOptions.microphonePositions(mic,d) <= 0 ||...
               experimentOptions.microphonePositions(mic,d) > experimentOptions.ismOptions.room(d)
               error('Microphones should be INSIDE the room.');
            end
        end
    end
end