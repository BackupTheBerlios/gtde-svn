function experimentOptions = GetSignals(experimentOptions)

%Check the options of the experiment
%
% USAGE: experimentOptions = GetSignals(experimentOptions)
%
% PARAMETERS:
%   experimentOptions ~ structure containing the options of the experiment
%
% DESCRIPTION: get the signals to experiment with
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
            NSignals = 1;
            signals = cell(NSignals,1);
            % Synthetic signals
            for ns = 1:NSignals,
                [sinusCoefficients sinusPhases sinusFrequencies] = GenerateLCSParameters(1/experimentOptions.dataOptions.samplingFrequency,200);
                signals{ns} = @(x) GenerateLCSSignal(x,sinusCoefficients,sinusPhases,sinusFrequencies);
            end  
        case 'simulated'
            % Retrieve folder
            [~, signalNames] = system(['find ' experimentOptions.dataOptions.wavFolder ' | grep .wav']);
            % Cut in names
            indices = strfind(signalNames,'.wav');
            % Number of signals, and cell declaration
            NSignals = numel(indices);
            signals = cell(NSignals,1);
            % Complete to ease the loop
            indices = cat(2,-4,indices);
            for kk = 1:numel(indices)-1,
                signals{kk} = signalNames(indices(kk)+5:indices(kk+1)+3);
            end
            rI = randi(NSignals,1);
            signals = signals([2,4,rI]);
        case 'real'
            signals = strcat('/scratch/pictor/deleforg/the_CASA_REDMINE/audio_recordings/',experimentOptions.field,'/sound12/Recorded/');
        otherwise
            error('Not ready yet!');
    end
    
    % Save it
    experimentOptions.signals = signals;
end