function [Positions options] = GeneratePositions(D,options)

%GeneratePositions generate positions in the D-dimensional space
%
% USAGE: [Positions, options] = GeneratePositions(D[, options])
%
% PARAMETERS:
%   D ~ dimension of the space
%   options ~ structure to specify some possible options:
%     coordinateSystem is a string with the values 'cartesian', 'polar',
%       'spherical' or 'cylindrical'. 'cartesian' is the default value, 
%       valid for any value of D. 'polar' is valid for D=2 (radius, angle)
%       and 'spherical' (radius, azimuth, elevation) and 'cylindrical' 
%       (planar radius, azimuth, high) are valid for D=3.
%
%     bounds is either a 1-element cell or a D-element cell, determining
%       the bounds of the points (in the chosen coordinate system, of
%       course). The default if [0,1] for the cartesian case, [0.1,1] for 
%       the radius [0,2*pi] for the polar angle and [-pi/2,pi/2] for the 
%       elevation angle. In the case of D-element cell, the user may
%       determine some of the bounds and let the other be the default 
%       ones. If the users wants the default value for some variables,
%       just let the d-cell empty.
%
%     numberOfIntervals is either a 1-element cell or a D-element cell,
%       determining the number of intervals PER COORDINATE PARAMETER. The
%       default number is 10 and it follows the same logic as 'bounds'.
%       If the users wants the default value for some variables, set
%       the d-th position of the array to 0.
%
%     uniformMaps is either a 1-element cell a D-element cell, determining
%       the space in which the intervals will be equidistant (by means of 
%       a 1D-to-1D function). The default is the identity and it follows 
%       the same logic as 'bounds'. NOTE: the function does not validate 
%       if the mapping can be applied to the specified bounds. If the
%       user wants the default value for some variables, let the d-cell
%       empty.
%
%     useExtrema is either a 1-element boolean or a D-element boolean,
%       determining if the extrema should be included or not. If the
%       extrema are included we will have NI+1 points, otherwise NI-1,
%       where NI is the number of intervals for this parameter. The
%       default value is true except for the elevation, which is false.
%
% RETURN VALUE:
%   Positions ~ NxD matrix in which each row is a position on the 
%     D-dimensional space
%   options ~ the options struct completed. That will include the following
%     fields:
%       uniformBounds is a D-element cell describing the bounds in the
%       uniforming space, i.e.:
%           uniformBounds{d} = uniformMaps{d}(bounds{d})
%
%   see also GenerateCartesianPositions, GenerateSphericalPositions,
%   GeneratePolarPositions, GenerateCylindricalPositions

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

    %%% Check input
    if nargin < 1
        error('Usage: Positions = GeneratePositions(D [, options])');
    end
    
    if nargin < 2
        options = struct([]);
    end

    if ~isposintscalar(D)
        error('The dimension parameter should be a positive integer.');
    end
    
    %%% Evaluate compatibility and set default values
    
    % Coordinate systems
    options = checkCoordinateSystem(D,options);
    
    % Bounds
    options = checkBounds(D,options);
    
    % Number of points
    options = checkNumberOfIntervals(D,options);
    
    % Uniform Map
    options = checkUniformMap(D,options);
    
    % useExtrema
    options = checkUseExtrema(D,options);
    
    %%% Compute the points
    
    % compute the bounds in the uniform space
    options = computeUniformBounds(D,options);
    
    % compute the parameter points
    options = computeParameterPoints(D,options);
    
    % map to the coordinate system
    switch options.coordinateSystem
        case 'polar'
            Positions = GeneratePolarPositions(options.parameterPoints{1},options.parameterPoints{2});
        case 'spherical'
            Positions = GenerateSphericalPositions(options.parameterPoints{1},options.parameterPoints{2},options.parameterPoints{3});
        case 'cylindrical'
            Positions = GenerateCylindricalPositions(options.parameterPoints{1},options.parameterPoints{2},options.parameterPoints{3});
        case 'cartesian'
            Positions = GenerateCartesianPositions(options.parameterPoints);
        otherwise
            error('Coordinate system not found.');
    end
    
    %%% Return
end


function options = checkCoordinateSystem(D,options)
    % Check if the coordinate system corresponds to the dimension
    if isfield(options,'coordinateSystem')
        switch options.coordinateSystem
            case 'polar'
                if D ~= 2
                    error('Polar coordinate systems are accepted at dimension D=2.');
                end
            case 'spherical'
                if D ~= 3
                    error('Spherical coordinate systems are accepted at dimension D=3.');
                end
            case 'cylindrical'
                if D ~= 3
                    error('Cylindrical coordinate systems are accepted at dimension D=3.');
                end
            case 'cartesian'
            otherwise
                error('coordinateSystem should be either "cartesian", "polar", "spherical" or "cylindrical".');
        end
    else
        % Default value
        options = struct('coordinateSystem','cartesian');
    end
end

function options = checkBounds(D,options)
    % Default options
    defaultRadia = [0.1 1];
    defaultCartesian = [0 1];
    defaultAzimuth = [0 2*pi];
    defaultElevation = [-pi/2 pi/2];
    % Check the size of bounds and set the defaults, if any
    if isfield(options,'bounds')
        if ~iscell(options.bounds)
            error('bounds should be of type cell');
        else
            if numel(options.bounds) ~= 1 && numel(options.bounds) ~= D
                error('bounds should be of size 1 or D.');
            elseif numel(options.bounds) == 1
                auxBounds = options.bounds{1};
                if ~checkB([auxBounds(1),auxBounds(2)])
                    error('The chosen bounds are not valid.');
                end
                options.bounds = cell(1,D);
                for d = 1:D
                    options.bounds{d} = auxBounds;
                end
            else
                switch options.coordinateSystem
                    % Check if defaults are needed
                    case 'polar'
                        if isempty(options.bounds{1})
                            options.bounds{1} = defaultRadia;
                        elseif numel(options.bounds{1}) ~= 2
                            error('When specifying the bounds of the radial component: it should be [lowerBound upperBound]');
                        elseif ~checkB(options.bounds{1})
                            error('The chosen bounds are not valid.');
                        end
                        if isempty(options.bounds{2})
                            options.bounds{2} = defaultAzimuth;
                        elseif numel(options.bounds{2}) ~= 2
                            error('When specifying the bounds of the angular component: it should be [lowerBound upperBound]');
                        elseif ~checkB(options.bounds{2})
                            error('The chosen bounds are not valid.');
                        end
                    case 'spherical'
                        if isempty(options.bounds{1})
                            options.bounds{1} = defaultRadia;
                        elseif numel(options.bounds{1}) ~= 2
                            error('When specifying the bounds of the radial component: it should be [lowerBound upperBound]');
                        elseif ~checkB(options.bounds{1})
                            error('The chosen bounds are not valid.');
                        end
                        if isempty(options.bounds{2})
                            options.bounds{2} = defaultAzimuth;
                        elseif numel(options.bounds{2}) ~= 2
                            error('When specifying the bounds of the azimuth component: it should be [lowerBound upperBound]');
                        elseif ~checkB(options.bounds{2})
                            error('The chosen bounds are not valid.');
                        end
                        if isempty(options.bounds{3})
                            options.bounds{3} = defaultElevation;
                        elseif numel(options.bounds{3}) ~= 2
                            error('When specifying the bounds of the elevation component: it should be [lowerBound upperBound]');
                        elseif ~checkB(options.bounds{3})
                            error('The chosen bounds are not valid.');
                        end
                    case 'cylindrical'
                        if isempty(options.bounds{1})
                            options.bounds{1} = defaultRadia;
                        elseif numel(options.bounds{1}) ~= 2
                            error('When specifying the bounds of the radial component: it should be [lowerBound upperBound]');
                        elseif ~checkB(options.bounds{1})
                            error('The chosen bounds are not valid.');
                        end
                        if isempty(options.bounds{2})
                            options.bounds{2} = defaultAzimuth;
                        elseif numel(options.bounds{2}) ~= 2
                            error('When specifying the bounds of the angular component: it should be [lowerBound upperBound]');
                        elseif ~checkB(options.bounds{2})
                            error('The chosen bounds are not valid.');
                        end
                        if isempty(options.bounds{3})
                            options.bounds{3} = defaultCartesian;
                        elseif numel(options.bounds{3}) ~= 2
                            error('When specifying the bounds of the high component: it should be [lowerBound upperBound]');
                        elseif ~checkB(options.bounds{3})
                            error('The chosen bounds are not valid.');
                        end
                    case 'cartesian'
                        for d = 1:D
                            if isempty(options.bounds{d})
                                options.bounds{d} = defaultCartesian;
                            elseif numel(options.bounds{d}) ~= 2
                                error(snprintf('When specifying the bounds of the %d-th component: it should be [lowerBound upperBound]',d));
                            elseif ~checkB(options.bounds{d})
                                error('The chosen bounds are not valid.');
                            end
                        end
                end
            end
        end
    else
        % Default values depending on the coordinate system
        switch options.coordinateSystem
            case 'polar'
                options.bounds = cell(1,2);
                options.bounds{1} = defaultRadia;
                options.bounds{2} = defaultAzimuth;
            case 'spherical'
                options.bounds = cell(1,3);
                options.bounds{1} = defaultRadia;
                options.bounds{2} = defaultAzimuth;
                options.bounds{3} = defaultElevation;
            case 'cylindrical'
                options.bounds = cell(1,3);
                options.bounds{1} = defaultRadia;
                options.bounds{2} = defaultAzimuth;
                options.bounds{3} = defaultCartesian;
            case 'cartesian'
                options.bounds = cell(1,D);
                for d = 1:D,
                    options.bounds{d} = defaultCartesian;
                end
        end
    end
end

function valid = checkB(bounds)
    valid = bounds(1) < bounds(2);
end

function options = checkNumberOfIntervals(D,options)
    % Default number of points
    defaultNumber = 10;
    % Check number of points and set defaults
    if isfield(options,'numberOfIntervals')
        if numel(options.numberOfIntervals) ~= 1 && numel(options.numberOfIntervals) ~= D
            error('NumberOfIntervals should be of size 1 or D.');
        elseif numel(options.numberOfIntervals) == 1
            auxNumber = options.numberOfIntervals;
            if ~isposintscalar(auxNumber)
                error('The specified number of points should be a positive integer (array).');
            end
            options.numberOfIntervals = auxNumber*ones(1,D);
        else
            for d = 1:D,
                if options.numberOfIntervals(d) == 0
                    options.numberOfIntervals(d) = defaultNumber;
                elseif ~isposintscalar(options.numberOfIntervals(d))
                    error('NumberOfIntervals should be a positive integer (array).');
                end
            end
        end
    else
        options.numberOfIntervals = defaultNumber*ones(1,D);
    end
end

function options = checkUniformMap(D,options)
    % Default map
    defaultMap = {@(x) x; @(x) x};
    % Check uniform maps
    if isfield(options,'uniformMaps')
        if ~iscell(options.uniformMaps)
            error('uniformMaps should be a cell');
        else
            if numel(options.uniformMaps) ~= 2 && numel(options.uniformMaps) ~= 2*D
                error('uniformMaps should be of size 2 or 2*D.');
            elseif numel(options.uniformMaps) == 2
                auxMap = options.uniformMaps;
                options.uniformMaps = cell(2,D);
                for d = 1:D,
                    options.uniformMaps{1,d} = auxMap{1};
                    options.uniformMaps{2,d} = auxMap{2};
                end
            else
                for d = 1:D
                    if isempty(options.uniformMaps{d})
                        options.uniformMaps{1,d} = defaultMap{1};
                        options.uniformMaps{2,d} = defaultMap{2};
                    end
                end
            end
        end
    else
        options.uniformMaps = cell(2,D);
        for d = 1:D
            options.uniformMaps{1,d} = defaultMap{1};
            options.uniformMaps{2,d} = defaultMap{2};
        end
    end
end

function options = checkUseExtrema(D,options)
    % Check use extrema
    if isfield(options,'useExtrema'),
        if ~isa(options.useExtrema,'logical')
            error('useExtrema must be a logical (array).');
        end
        if numel(options.useExtrema) ~= 1 && numel(options.useExtrema) ~= D
            error('useExtrema should be of size 1 or D.');
        elseif numel(options.useExtrema) == 1
            auxUseExtrema = options.useExtrema(1);
            if auxUseExtrema == false && sum(options.numberOfIntervals == 1)
                error('You need to use the extrema for those dimensions with one interval.');
            end
            options.useExtrema = false(1,D);
            for d = 1:D
                options.useExtrema(d) = auxUseExtrema;
            end            
        end
    else
        % Default use extrema
        options.useExtrema = true(1,D);
    end
end

function options = computeUniformBounds(D,options)
    % Set the structure field
    options.uniformBounds = cell(1,D);
    % Compute them
    for d = 1:D,
        uniformMap = options.uniformMaps{1,d};
        options.uniformBounds{d} = [uniformMap(options.bounds{d}(1)) uniformMap(options.bounds{d}(2))];
    end
end

function options = computeParameterPoints(D,options)
    % Set the structure field
    options.parameterPoints = cell(1,D);
    % Compute them
    for d = 1:D,
        uniformIntervals = (options.uniformBounds{d}(2)-options.uniformBounds{d}(1))/options.numberOfIntervals(d);
        uniformPoints = options.uniformBounds{d}(1) + uniformIntervals*(1:(options.numberOfIntervals(d)-1));
        if options.useExtrema(d)
            uniformPoints = cat(2,options.uniformBounds{d}(1),uniformPoints,options.uniformBounds{d}(2));
        end
        options.parameterPoints{d} = options.uniformMaps{2,d}(uniformPoints);
    end
end