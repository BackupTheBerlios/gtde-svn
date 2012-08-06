function experimentOptions = SpecifySourcePositionOptions(experimentOptions)

%Check the source position options of the experiment
%
% USAGE: experimentOptions = SpecifySourcePositionOptions(experimentOptions)
%
% PARAMETERS:
%   experimentOptions ~ structure containing the options of the experiment.
%     The structure should containt a sub-structure called
%     "sourcePositionOptions". This structure follows the syntax specified
%     in GeneratePositions. The only check done here is the existence of
%     such structure and the dimension coherence.
%
% DESCRIPTION: the options are checked, if something fails an error is
%   reported.
%
% see also GeneratePositions

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
    if ~isfield(experimentOptions,'sourcePositionOptions')
        error('The field "sourcePositionOptions" should be specified in the experimentOtions structure.');
    end
    
    % Check the dimension
    if ~isfield(experimentOptions.sourcePositionOptions,'dimension')
        error('The field "dimension" is not specified in "experimentOptions.sourcePositionOptions", see GeneratePositions for more information.');
    end
    if experimentOptions.dimension ~= experimentOptions.sourcePositionOptions.dimension
        fprintf('Warning: The dimenion of the experiment and of the sound source space are not the same. Taking the dimension of the experiment.');
        experimentOptions.sourcePositionOptions.dimension = experimentOptions.dimension;
    end

end