function Sh = GeneratePositionsCasaRedmine(panIndices, tiltIndices, field)
% GeneratePositionsCasaRedmine
%
%   [microphonesPositions, azimuths, elevations, radia] = GeneratePositionsCasaRedmine(panIndices, tiltIndices, field)
%
%   panIndices, tiltIndices  - Indices for the motor states of the robot
%
%     microphonesPositions - microphone positions in the robot's head coordinate system (does
%         not depend on the pan and tilt of course)
%     azimuts, elevations, radia - azimuth, elevation, radius w.r.t. robot's head coordinate
%                      system
%     
% THE DOC HAS NOT BEEN DONE YET

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

%%% Importdata
data = importdata(strcat('/scratch/Andromeda/amat/WorkingDevLocal/CASA_RED_MINE_19122011/snapshots/',field,'/points_3d_new_coord_syst.txt'));
% Microphones' positions and speaker's position
p = data(:,2:end);

% microphonesPositions = zeros(4,3);
% microphonesPositions(1,:) = p(3,:);
% microphonesPositions(2,:) = p(4,:);
% microphonesPositions(3,:) = p(2,:);
% microphonesPositions(4,:) = p(1,:);

Sw = p(5,:)'; %coordinates w.r.t. World

%%% Load pan/tilt values
data = importdata(strcat('/scratch/pictor/deleforg/the_CASA_REDMINE/audio_recordings/',field,'/motorsGT.txt'));
data(1,:) = [];
motorStates = zeros(numel(data)/2,2);
motorStates(:,1) = data(1:2:end);
motorStates(:,2) = data(2:2:end);
% Construct the indices
pairIndices1 = false(size(motorStates,1),1);
pairIndices2 = false(size(motorStates,1),1);
for ii=1:numel(panIndices)
    pairIndices1(panIndices(ii):37:end) = true;
end
for ii=1:numel(tiltIndices),
    pairIndices2( (37*(tiltIndices(ii)-1)+1):37*tiltIndices(ii) ) = true;
end
pairIndices = pairIndices1 & pairIndices2;
% Get pans and tils
phi = motorStates(pairIndices,1)*pi/180;
psi = motorStates(pairIndices,2)*pi/180;

% Sound source position
Sh = zeros(numel(phi),3);
for sPos = 1:numel(phi),
    R_pan = [cos(phi(sPos)), sin(phi(sPos)), 0; ...
            -sin(phi(sPos)), cos(phi(sPos)), 0; ...
             0       , 0       , 1];

    R_tilt = [1,         0,        0; ...
              0,  cos(psi(sPos)), sin(psi(sPos)); ...
              0, -sin(psi(sPos)), cos(psi(sPos))];


    Sh(sPos,:) = R_tilt * R_pan * Sw;  %coordinates w.r.t. Robot head
end

% [azimuths, elevations, radia] = cartesian2spherical(Sh(:,1),Sh(:,2),Sh(:,3));


