function experimentOptions = GetSourcePositions(experimentOptions)

%Get the source positions depending on the options of the experiment.
%
% USAGE: experimentOptions = GetSourcePositions(experimentOptions)
%
% PARAMETERS:
%   experimentOptions ~ the options of the experiment. If
%     experimentOptions.sourcePositions does not exist, you should specify
%     experimentOptions.sourcePositionOptions.
% 
%   see also GeneratePositions, GeneratePositionsCasaRedmine

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

    if isfield(experimentOptions,'sourcePositions')
        fprintf('Warning: The field "sourcePositions" is set, clearing the field "sourcePositionOptions".\n');
        experimentOptions.sourcePositionOptions = [];
    else
        % Choose depending on the data used for the experiment
        switch experimentOptions.dataOptions.type
            case 'real'
                experimentOptions.sourcePositions = GeneratePositionsCasaRedmine(experimentOptions);
            case 'simulated'
                if isfield(experimentOptions.sourcePositionOptions,'file'),
                    tmp = importdata(experimentOptions.sourcePositionOptions.file,'\t',2);
                    experimentOptions.sourcePositions = tmp.data(:,2:4) + repmat(experimentOptions.microphonePositionOptions.offset,size(tmp.data,1),1);
                else
                    experimentOptions.sourcePositions = GeneratePositions(experimentOptions.dimension,experimentOptions.sourcePositionOptions);
                end
            case 'synthetic'
                experimentOptions.sourcePositions = GeneratePositions(experimentOptions.dimension,experimentOptions.sourcePositionOptions);
        end
    end

end