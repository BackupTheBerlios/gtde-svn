function experimentOptions = SpecifyMethodOptions(experimentOptions)

%Check the source position options of the experiment
%
% USAGE: experimentOptions = SpecifyMethodOptions(experimentOptions)
%
% PARAMETERS:
%   experimentOptions ~ structure containing the options of the experiment.
%     The structure should containt a sub-structure called
%     "methodOptions". That sub-structure must contain the following points:
%     
%       type ~ describing the method used: 'gtde' for
%       geometrically-constrained tde, 'tde' for unconstrained tde, 'init'
%       withough local minimization, 'bypairs' independent tde estimation.
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

    % If ISM is not used just return
    if ~isfield(experimentOptions,'methodOptions')
        error('The "methodOptions" substructure is needed.');
    end
    
    % Mandatory fields
    if ~isfield(experimentOptions.methodOptions,'type')
        error('The method type should be specified in the method options substructure.');
    end
    
    method =  {'gtde','tde','init','bypairs','truth'};
    % Check the method specified
    methodFlag = 0;
    for d = 1:numel(method)
        methodFlag = methodFlag + strcmp(experimentOptions.methodOptions.type,method{d});
    end
    if methodFlag == 0
        error('Method specified not known.');
    end
    
    % Grid
    if ~isfield(experimentOptions.methodOptions,'gridSize'),
        experimentOptions.methodOptions.gridSize = 15;
    end
end
