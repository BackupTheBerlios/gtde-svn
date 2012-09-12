function ISMFileName = GetISMFileName(experimentOptions,sT60,sPos,sSF)

%Build the ISM file name
%
% USAGE: ISMFileName = GetISMFileName(experimentOptions,sT60,sPos,sSF)
%
% PARAMETERS:
%   experimentOptions ~ we need the ismOption.ismFolder
%   t60 ~ reverberation time index
%   pos ~ sound source position index
%   SF ~ sampling frequency index
%
% RETURN VALUE: is the filename for the ISM filters in the format:
%  rootFolder/ismFolder/ism_rir_t60_SF_pos(1)_pos(2)_pos(3).mat
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

    % Build ISM file name
    ISMFileName = strcat(experimentOptions.ismOptions.folder,...
                     'ismFilter_',...
                     num2str(sT60),'_',...
                     num2str(sSF),'_',...
                     num2str(sPos),'.mat');
end