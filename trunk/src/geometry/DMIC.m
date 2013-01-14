function D = DMIC(X,MIC)

%Get the maximum possible value of the TDE given the microphones' position
%
% DMIC(X,MIC) returns the Euclidean distance of each row of X to the 
%       microphone's position MIC
%
% PARAMETERS:
%   X ~ Positions of the source
%   MIC ~ Position of the microphones

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

    %%% Check the input
    % The dimension of both vectors should be the same
    if(size(X,2)~=size(MIC,2))
        error('X and MIC1 should have the same number of columns.');
    end
    % There should be just one microphone position
    if(size(MIC,1)~=1)
        fprintf('[DMIC Warning]: Computing the distance to the first row of MIC.');
    end

    %%% Compute the distance of every row on X to the first row of MIC
    D = sqrt(sum(bsxfun(@minus,X,MIC(1,:)).^2,2));
    
return
