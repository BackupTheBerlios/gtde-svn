function experimentOptions = GetMicrophonePositions(experimentOptions)

%Get the microphone positions depending on the options of the experiment.
%
% USAGE: experimentOptions = GetMicrophonePositions(experimentOptions)
%
% PARAMETERS:
%   experimentOptions ~ the options of the experiment. If
%     experimentOptions.microphonePositions does not exist, you should specify
%     experimentOptions.microphonePositionOptions. The
%     microphonePositionOptions structure may have the following fields (M
%     means mandatory)
%
%       type [M] ~ Type of microphone array used. 'tetrahedon' is the only
%         option right now.
%       scale ~ The size of the array (Default = 0.1m).
%       offset ~ The offset of the array with respect to the origin
%         (Default = null vector).
%


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

    if isfield(experimentOptions,'microphonePositions')
        fprintf('Warning: The field "microphonePositions" is set, clearing the field "microphonePositionOptions".\n');
        experimentOptions.microphonePositionOptions = [];
    else
        % Check different cases
        if ~strcmp(experimentOptions.dataOptions.type,'real')
            % Just a shorcut to the scale
            a = experimentOptions.microphonePositionOptions.scale;
            % Switch the type of microphone array
            switch experimentOptions.microphonePositionOptions.type
                case 'tetrahedron'
                    % Dimension is either 2 or 3
                    if experimentOptions.dimension == 2
                        % Scale parameters
                        b = -a/2;
                        c = a*sqrt(9/12);
                        % Microphones in a regular tetrahedron
                        experimentOptions.microphonePositions = ...
                            [ a   0 ;...
                              b   c ;...
                              b  -c];
                    else
                        % Scale parameters
                        c = a / sqrt(2);
                        % Microphones in a regular tetrahedron
                        experimentOptions.microphonePositions = ...
                            [ a   0   -c ;...
                             -a   0   -c ;...
                              0   a    c ;...
                              0  -a    c];

                    end
                case 'icosahedron'
                     mP = GenerateIcosahedron(a);
                     % If the user desires to use a subset of the
                     % microphones, extract this subset.
                    if isfield(experimentOptions.microphonePositionOptions,'subSet')
                        mP = mP(experimentOptions.microphonePositionOptions.subSet,:);
                    end
                    experimentOptions.microphonePositions = mP;
            end
            % If the field offset exists, add it to the microphone
            % positions
            if isfield(experimentOptions.microphonePositionOptions,'offset')
                experimentOptions.microphonePositions = experimentOptions.microphonePositions + ...
                    repmat(experimentOptions.microphonePositionOptions.offset,size(experimentOptions.microphonePositions,1),1);
            end
        elseif strcmp(experimentOptions.dataOptions.type,'real')
            error('Warning: Microphone positions not ready on real data.');
%             %%% Importdata
%             data = importdata(strcat('/scratch/Andromeda/amat/WorkingDevLocal/CASA_RED_MINE_19122011/snapshots/',experimentOptions.field,'/points_3d_new_coord_syst.txt'));
%             % Microphones' positions and speaker's position
%             p = data(:,2:end);
%             % Allocate microphones positions
%             experimentOptions.microphonesPositions = zeros(4,3);
%             experimentOptions.microphonesPositions(1,:) = p(3,:);
%             experimentOptions.microphonesPositions(2,:) = p(4,:);
%             experimentOptions.microphonesPositions(3,:) = p(2,:);
%             experimentOptions.microphonesPositions(4,:) = p(1,:);        
        end
    end

end