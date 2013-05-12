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
    initializationPositionOptions.useExtrema=false;
    % Number of intervals
    initializationPositionOptions.numberOfIntervals = experimentOptions.methodOptions.gridSize;
    % Generate positions
    initializationPositions = GeneratePositions(initializationPositionOptions.dimension,initializationPositionOptions)';
    % Remove those who are outside the bounds
    initializationPositions = DiscardOutBounds(initializationPositions,experimentOptions.microphonePositions);

    % Modify the initial positions depending on the method used
    switch experimentOptions.methodOptions.type
        % Nothing to do for these methods
        case {'init','tde','m1qn3'}
        % Remove those that do not satisfy the constraint
        case {'dip','sqplab'}
            Constraint = TDEDiscriminant(initializationPositions,experimentOptions.microphonePositions);
            initializationPositions = initializationPositions(:,Constraint>0);
        % The essence is the same as in the first case, the structure
        % different.
        case {'bypairs','bp2dip'}
            auxInit = initializationPositions;
            % Reallocate them
            initializationPositions = cell(1,experimentOptions.dimension);
            % Chose the unique values for each column
            for d = 1:experimentOptions.dimension,
                initializationPositions{d} = unique(auxInit(d,:));
            end
        % Much more dense data
        case 'truth'
            % Interval resolution, the error should be lower than half of this quantity
            resolution = 2e-5;
            % Number of intervals (per TDE pair)
            initializationPositionOptions.numberOfIntervals = zeros(initializationPositionOptions.dimension,1);
            for d = 1:initializationPositionOptions.dimension
                initializationPositionOptions.numberOfIntervals(d) = ceil(2*maxTDEs(d)/resolution);
            end
            initializationPositions = GeneratePositions(initializationPositionOptions.dimension,initializationPositionOptions)';
            % Remove those who are outside the bounds
            initializationPositions = DiscardOutBounds(initializationPositions,experimentOptions.microphonePositions);
    end
    
    % Save it!
    experimentOptions.initializationPositions = initializationPositions;
end
