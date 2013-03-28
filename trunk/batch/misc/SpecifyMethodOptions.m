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
    
    method =  {'dip','init','bypairs','truth','fminunc','eusipco','bab','bp2dip','bab2dip'};
    % Check the method specified
    methodFlag = 0;
    for d = 1:numel(method)
        methodFlag = methodFlag + strcmp(experimentOptions.methodOptions.type,method{d});
    end
    if methodFlag == 0
        error('Method specified not known.');
    end
    
    % If sqplab, check software
    if strcmp(experimentOptions.methodOptions.type,'sqplab'),
        if isempty(which('sqplab'))
            fprintf('SQPLAB not found. Please go to the following link and install it or add the right path.\n');
            fprintf('https://who.rocq.inria.fr/Jean-Charles.Gilbert/modulopt/optimization-routines/sqplab/sqplab.html');
            error('SQPLAB Software not found.');
        end
    end
    
    % If m1qn3, check software
    if strcmp(experimentOptions.methodOptions.type,'m1qn3'),
        if isempty(which('m1qn3'))
            fprintf('M1QN3 not found. Please go to gtde/trunk/m1qn3_wrapper and run "sh compile.sh".\n');
            fprintf('Set the LD_LIBRARY_PATH properly and restart MATLAB.');
            error('M1QN3 Software not found.');
        end
    end
    
    % Check for the regularization constant on fminunc
    if strcmp(experimentOptions.methodOptions.type,'fminunc'),
        if ~isfield(experimentOptions.methodOptions,'regConstant')
            experimentOptions.methodOptions.regConstant = 10;
        end
    end
    
    % Grid
    if ~isfield(experimentOptions.methodOptions,'gridSize'),
        experimentOptions.methodOptions.gridSize = 15;
    end
    % Verbose
    if ~isfield(experimentOptions.methodOptions,'verbose'),
        experimentOptions.methodOptions.verbose = 0;
    end
    % Output
    if ~isfield(experimentOptions.methodOptions,'output'),
        experimentOptions.methodOptions.output = 0;
    end
end
