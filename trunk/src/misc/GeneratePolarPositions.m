function Positions = GeneratePolarPositions(Radia,Azimuths)

%Generate positions lying in circles expressed in cartesian coordinates
%
% USAGE: Positions = GeneratePolarPositions(Radia,Azimuths)
%
% PARAMETERS:
%   Radia ~ Vector of size R with the radia values in the range (0,infty)
%   Azimuths ~ Vector of size A with the azimuth values in the range
%   [0,2*pi]
%
% RETURN VALUE:
%   Positions ~ Matrix of size R*A-by-2 with the positions
%
% EXAMPLE:
%   Positions = GeneratePolarPositions([3,4],[0,pi]);
%   In that case Positions will be a 4-by-3 matrix with the cartesian
%   coordinates of the points lying the specified radia and azimuths.
%
%   see also GeneratePositions, pol2cart

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
    if nargin < 2
        error('Usage: Positions = GeneratePositions(Radia, Azimuths)');
    end

    if isempty(Radia)
        error('Radia cannot be empty');
    end

    if isempty(Azimuths)
        error('Azimuths cannot be empty');
    end
    
    if sum(Radia < 0)
        error('Radia should not contain negative values.');
    end
    
    %%% Compute
    % Useful numbers
    NRadia = numel(Radia);
    NAzimuths = numel(Azimuths);
    
    Positions = zeros(NRadia*NAzimuths,2);
    % Position index
    iPositions = 1;
    for iRadia = 1:NRadia,
        for iAzimuth = 1:NAzimuths,
            [Positions(iPositions,1) Positions(iPositions,2)] = pol2cart(Azimuths(iAzimuth),Radia(iRadia));
            iPositions = iPositions + 1;
        end
    end

end