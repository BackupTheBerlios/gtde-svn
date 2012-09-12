function experimentOptions = SpecifyResultOptions(experimentOptions)

%Check the source position options of the experiment
%
% USAGE: experimentOptions = SpecifyResultOptions(experimentOptions)
%
% PARAMETERS:
%   experimentOptions ~ structure containing the options of the experiment.
%     The structure should containt a sub-structure called
%     "methodOptions". That sub-structure must contain the following points:
%     
%       folder ~ absolute folder [Default = strcat(rootFolder,'../results/')]
%       suffix ~ string that will be used as a suffix [Default = '']
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
    if ~isfield(experimentOptions,'resultOptions')
        % Default options
        resultOptions.folder = strcat(experimentOptions.rootFolder,'results/');
        resultOptions.suffix = '';
        experimentOptions.resultOptions = resultOptions;
    else
        % If some of the field of the strucutre are not set, choose the
        % default values
        if ~isfield(experimentOptions.resultOptions,'folder')
            experimentOptions.resultOptions.foder = strcat(experimentOptions.rootFolder,'results/');
        end
        if ~isfield(experimentOptions.resultOptions,'suffix')
            experimentOptions.resultOptions.suffix = '';
        end
    end
    
    % Check if directory exists
    resultsFolder = strcat(experimentOptions.resultOptions.folder);
    if ~exist(resultsFolder,'dir'),
        fprintf('The results folder does not exist, please create it.\n%s',resultsFolder);
        error('Missing resutls folder.');
    end
    
    % Check for permissions
    currentDir = cd;
    cd(resultsFolder);
    [success message] = fileattrib;
    if success
        if message.UserWrite ~= true
            error(['You do not have writing permissions in the results directory. ' ...
                   'The computation will be lost. ' ...
                   'Please check with your system administrator or change the results folder.']);
        end
    else
        fprintf('Warning: an error ocurred retrieving the properties of the results folder. The computation may be lost.');
    end
    cd(currentDir);
    
end