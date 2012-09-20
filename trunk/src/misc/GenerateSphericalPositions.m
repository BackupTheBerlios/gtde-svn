function Positions = GenerateSphericalPositions(Radia,Azimuths,Elevations)

%Generate positions lying in spheres expressed in cartesian coordinates
%
% USAGE: Positions = GenerateSphericalPositions(Radia,Azimuths,Elevations)
%
% PARAMETERS:
%   Radia ~ Vector of size R with the radia values in the range (0,infty)
%   Azimuths ~ Vector of size A with the azimuth values in the range
%   [0,2*pi]
%   Elevations ~ Vector of size E with the elevation values in the range [-pi/2,pi/2]
%
% RETURN VALUE:
%   Positions ~ Matrix of size R*A*E-by-3 with the positions
%
% EXAMPLE:
%   Positions = GenerateSphericalPositions([3,4],[0,pi],[-pi/4,pi/4]);
%   In that case Positions will be a 8-by-3 matrix with the cartesian
%   coordinates of the points lying the the radia, azimuth and elevation
%   specified.
%
%   see also GeneratePositions, sph2cart

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
    if nargin < 3
        error('Usage: Positions = GeneratePositions(Radius, Azimuths, Evelations)');
    end

    if isempty(Radia)
        error('Radia cannot be empty');
    end

    if isempty(Azimuths)
        error('Azimuths cannot be empty');
    end

    if isempty(Elevations)
        error('Elevations cannot be empty');
    end
    
    if sum(Radia < 0)
        error('Radia should not contain negative numbers.');
    end

    NRadia = numel(Radia);
    NAzimuths = numel(Azimuths);
    NElevations = numel(Elevations);
    

    %%% Generate positions
    Positions = zeros(NRadia*NAzimuths*NElevations,3);
    iPosition = 1;
    % Looping
    for iRadius = 1:NRadia,
        radius = Radia(iRadius);
        for iAzimuth = 1:NAzimuths,
            azimuth = Azimuths(iAzimuth);
            for iElevation = 1:NElevations,
                elevation = Elevations(iElevation);
                % Generate the position
                [Positions(iPosition,1) Positions(iPosition,2) Positions(iPosition,3)] = ...
                    sph2cart(azimuth, elevation, radius);
                % Increase the counter
                iPosition = iPosition + 1;
            end
        end
    end
    
end