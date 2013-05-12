function experimentOptions = SpecifyISMFilterOptions(experimentOptions)

%Check the source position options of the experiment
%
% USAGE: experimentOptions = SpecifyISMFilterOptions(experimentOptions)
%
% PARAMETERS:
%   experimentOptions ~ structure containing the options of the experiment.
%     The structure should containt a sub-structure called
%     "ismOptions". That sub-structure must contain the following points:
%     
%       ismFolder [M] ~ folder in which the ism data is recorded, not to
%         compute it every time, if it exists, we load the ISM data from
%         this folder, otherwise we compute ISM data (absolute path)
%       absorptionWeights ~ the relative absorption weights of the walls 
%       room ~ size of the room
%       t60 ~ values for t60
%
% DESCRIPTION: the options are checked, if something fails an error is
%   reported.
%
% REFERENCES:
%   [1] E. A. Lehmann. Matlab code for image-source model in room acoustics.
%       http://www.eric-lehmann.com/ism code.html.
%   [2] E. A. Lehmann and A. M. Johansson. Prediction of energy decay in 
%       room impulse responses simulated with an image-source model. The 
%       Journal of the Acoustical Society of, 124(1):269–277, 2008.

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

    % If ISM is not used just return
    if ~experimentOptions.ism
        return;
    end
    
    % If real data do not bother with ISM
    if strcmp(experimentOptions.dataOptions.type,'real'),
        fprintf('Warning: ISM filters are not used in real data.\n');
        experimentOptions.ism = false;
        return;
    end
    
    % Checking for the structure
    if ~isfield(experimentOptions,'ismOptions')
        error('The use of the ISM model requires to specify the ismOptions substructure.');
    end
    
    if experimentOptions.dimension == 2 && experimentOptions.ism
        error('Not using ISM in the 2D case.');
    end
    
    % If the ismFolder exists, we load the ism data from it
    if ~isfield(experimentOptions.ismOptions,'folder'),
        error('The "folder" field should be specified in the "ismOptions" structure.');
    end
    
    % Do we have to generate ISM
    generateISM = false;
    if ~exist(experimentOptions.ismOptions.folder,'dir'),
        generateISM = true;
    end
    
    % Check what to do in both cases
    if ~generateISM
        % The ISM folder exist (the user wants to use some precomputed ISM filters)
        tmp = load(strcat(experimentOptions.ismOptions.folder,'experimentOptions.mat'));
        % Copy the ISM, source positions and microphone positions
        experimentOptions.microphonePositionOptions = tmp.experimentOptions.microphonePositionOptions;
        experimentOptions.sourcePositionOptions = tmp.experimentOptions.sourcePositionOptions;
        experimentOptions.microphonePositions = tmp.experimentOptions.microphonePositions;
        experimentOptions.sourcePositions = tmp.experimentOptions.sourcePositions;
        % Keep the wanted folder (for different machine processing
        % purposes).
        tmp.experimentOptions.ismOptions.folder = experimentOptions.ismOptions.folder;
        experimentOptions.ismOptions = tmp.experimentOptions.ismOptions;
    else
        % The ISM folder does not exist (we will need to compute them)
        % Check some mandatory fields
        mandatoryFields = {'absorptionWeights','t60','room','samplingFrequencies'};
        for f = 1:numel(mandatoryFields),
            if ~isfield(experimentOptions.ismOptions,mandatoryFields{f})
                error(['Field ' mandatoryFields{f} ' is mandatory in the structure "ismOptions".']);
            end
        end
    end
    experimentOptions.ismOptions.generate = generateISM;

end