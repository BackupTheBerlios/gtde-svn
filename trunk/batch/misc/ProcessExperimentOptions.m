function experimentOptions = ProcessExperimentOptions(experimentOptions)

%Check the options of the experiment
%
% USAGE: experimentOptions = ProcessExperimentOptions(experimentOptions)
%
% PARAMETERS:
%   experimentOptions ~ structure containing the options of the experiment
%
% DESCRIPTION: the options are checked, if something fails an error is
%   reported.
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

    % Check that there is the kind of data used and the type of sensor
    mandatoryFields = {'dimension','ism','rootFolder',...
        'sourcePositionOptions','microphonePositionOptions',...
        'ismOptions','dataOptions','methodOptions'};
    for f = 1:numel(mandatoryFields),
        if ~isfield(experimentOptions,mandatoryFields{f})
            error(['Field ' mandatoryFields{f} ' is mandatory.']);
        end
    end
    
    % Check the dimension
    if experimentOptions.dimension ~= 2 && experimentOptions.dimension ~= 3
        error('Dimension should be 2 or 3.');
    end
    
    % Source position options
    experimentOptions = SpecifySourcePositionOptions(experimentOptions);

    % Data options
    experimentOptions = SpecifyDataOptions(experimentOptions);
    
    % Microphone position options
    experimentOptions = SpecifyMicrophonePositionOptions(experimentOptions);
    
    % ISM Filtering options
    experimentOptions = SpecifyISMFilterOptions(experimentOptions);
    
    % Method options
    experimentOptions = SpecifyMethodOptions(experimentOptions);
    

end