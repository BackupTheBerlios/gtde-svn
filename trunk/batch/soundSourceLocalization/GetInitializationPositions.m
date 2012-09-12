function experimentOptions = GetInitializationPositions(experimentOptions)

%Get the initialization positions for the optimization procedure
%
% USAGE: experimentOptions = GetInitializationPositions(experimentOptions)
%
% PARAMETERS:
%   experimentOptions ~ structure containing the options of the experiment
%
% DESCRIPTION: the user can specify the grid, in methodOptions.
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

    % Dimension
    initializationPositionOptions.dimension = size(experimentOptions.microphonePositions,1)-1;
    % Coordinate system
    initializationPositionOptions.coordinateSystem = 'cartesian';
    % Bounds
    maxTDEs = TDEmax(experimentOptions.microphonePositions);
    initializationPositionOptions.bounds = cell(initializationPositionOptions.dimension,1);
    for d = 1:initializationPositionOptions.dimension
        initializationPositionOptions.bounds{d} = [-maxTDEs(d),maxTDEs(d)];
    end
    % Number of intervals
    initializationPositionOptions.numberOfIntervals = experimentOptions.methodOptions.gridSize;
    % Generate positions
    initializationPositions = GeneratePositions(initializationPositionOptions.dimension,initializationPositionOptions);

    % Modify the initial positions depending on the method used
    switch experimentOptions.methodOptions.type
        case 'gtde'
            Constraint = TDEDiscriminant(initializationPositions',experimentOptions.microphonesPositions);
            initializationPositions = initializationPositions(Constraint>0,:);
        case 'tde'
        case 'init'
        case 'bypairs'
            auxInit = initializationPositions;
            % Reallocate them
            initializationPositions = cell(1,experimentOptions.dimension);
            % Chose the unique values for each column
            for d = 1:experimentOptions.dimension,
                initializationPositions{d} = unique(auxInit(:,d));
            end
    end
    
    % Save it!
    experimentOptions.initializationPositions = initializationPositions;
end