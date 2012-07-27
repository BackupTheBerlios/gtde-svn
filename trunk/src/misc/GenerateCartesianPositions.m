function Positions = GenerateCartesianPositions(Coordinates)

%Generate positions in cartesian coordinates
%
% USAGE: Positions = GenerateCartesianPositions(Radia,Azimuths,Elevations)
%
% PARAMETERS:
%   Coordinates ~ Cell containing the points of the grid for each
%   dimension. The dimension will be numel(Coordinates).
%
% RETURN VALUE:
%   Positions ~ Matrix of size N-by-D with the positions.
%
% EXAMPLE:
%   Coordinates = {0:0.1:1; 2:0.1:4};
%   Positions = GenerateCartesianPositions(Coordinates);
%   In that case Positions will be a N-by-2 matrix with the cartesian
%   coordinates of the points lying the specified grid.
%
%   see also GeneratePositions

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

    %%% Main variables
    % Dimension
    D = numel(Coordinates);
    % Max block size
    BlockSize = 1;
    for d = 1:D
        BlockSize = BlockSize * numel(Coordinates{d});
    end
    % How many times do we repeat
    RepeatBlock = 1;

    %%% Compute positions
    % Set up the variables
    Positions = zeros(BlockSize,D);
    % Fill each dimension
    for d = 1:D
        % Create the block
        block = zeros(BlockSize,1);
        % Element size
        ElementSize = BlockSize / numel(Coordinates{d});
        % Fill the block
        for e = 1:numel(Coordinates{d}),
            block( (e-1)*ElementSize+1:e*ElementSize ) = ...
                Coordinates{d}(e);
        end
        Positions(:,d) = repmat(block,RepeatBlock,1);
        % Update repeat block and the block size
        RepeatBlock = RepeatBlock * numel(Coordinates{d});
        BlockSize = ElementSize;
    end
    
end
    