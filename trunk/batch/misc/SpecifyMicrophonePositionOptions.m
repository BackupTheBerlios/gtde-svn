function experimentOptions = SpecifyMicrophonePositionOptions(experimentOptions)

%Check the source position options of the experiment
%
% USAGE: experimentOptions = SpecifyMicrophonePositionOptions(experimentOptions)
%
% PARAMETERS:
%   experimentOptions ~ structure containing the options of the experiment.
%     The structure should containt a sub-structure called
%     'microphonePositionOptions' with the field 'type' set to
%     'tetrahedron', and the fields scale and offset properly set.
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
    if ~isfield(experimentOptions,'microphonePositionOptions')
        error('The field "dataOptions" should be specified in the experimentOtions structure.');
    end
    
    % Check the data specified
    mandatoryFields = {'type'};
    for f = 1:numel(mandatoryFields),
        if ~isfield(experimentOptions.microphonePositionOptions,mandatoryFields{f})
            error(['Field ' mandatoryFields{f} ' is mandatory in the structure "microphonePositionOptions".']);
        end
    end
    
    % Loop searching for the data type choice
    typeOptions= {'tetrahedron'};
    dataFlag = 0;
    for d = 1:numel(typeOptions)
        dataFlag = dataFlag + strcmp(experimentOptions.microphonePositionOptions.type,typeOptions{d});
    end
    if dataFlag == 0
        error('Microphone array type specified not known.');
    end
    
    % Particular fields needed
    if strcmp(experimentOptions.microphonePositionOptions.type,'tetrahedron')
        mandatoryFields = {'scale','offset'};
        for f = 1:numel(mandatoryFields),
            if ~isfield(experimentOptions.microphonePositionOptions,mandatoryFields{f})
                error(['Field ' mandatoryFields{f} ' is mandatory when using a tetrahedron array.']);
            end
        end
    end
end