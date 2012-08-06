function experimentOptions = SpecifyMicrophonePositionOptions(experimentOptions)

%Check the source position options of the experiment
%
% USAGE: experimentOptions = SpecifyMicrophonePositionOptions(experimentOptions)
%
% PARAMETERS:
%   experimentOptions ~ structure containing the options of the experiment.
%     The structure should containt a sub-structure called
%     "dataOptions". That sub-structure must contain the following points:
%     type ~ either synthetic, simulated or real. Describes the data used.
%
%     if type == synthetic, the field samplingFrequency is mandatory.
%     if type == simulated, the field wavFolder is mandatory. This is the 
%       folder relative to rootFolder containing the emitted signals.
%     if type == real, the structure experimentOptions.sourcePositionOptions 
%       must contain the field, the panIndices and the tiltIndices to 
%       retrieve the sound source position from a file.
%
% DESCRIPTION: the options are checked, if something fails an error is
%   reported.

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

    % Input check
    if ~isfield(experimentOptions,'dataOptions')
        error('The field "dataOptions" should be specified in the experimentOtions structure.');
    end
    
    % Check the data specified
    mandatoryFields = {'type'};
    for f = 1:numel(mandatoryFields),
        if ~isfield(experimentOptions,mandatoryFields{f})
            error(['Field ' mandatoryFields{f} ' is mandatory in the structure "dataOptions".']);
        end
    end
    
    % Loop searching for the data type choice
    dataTypeOptions = {'tetrahedron'};
    dataFlag = 0;
    for d = 1:numel(dataTypeOptions)
        dataFlag = dataFlag + strcmp(experimentOptions.dataOptions.type,dataTypeOptions{d});
    end
    if dataFlag == 0
        error('Microphone array type specified not known.');
    end
    
    % Particular fields needed
    if strcmp(experimentOptions.dataOptions.type,'tetrahedron')
        mandatoryFields = {'scale','offset'};
        for f = 1:numel(mandatoryFields),
            if ~isfield(experimentOptions.dataOptions,mandatoryFields{f})
                error(['Field ' mandatoryFields{f} ' is mandatory when using synthetic data.']);
            end
        end
    end
end