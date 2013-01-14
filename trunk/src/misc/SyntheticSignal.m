function y = SyntheticSignal(x,choice,params)

%Some synthetic basic signals
%
% DESCRIPTION:
%   This is used for test purposes. Not really worth to document.
%
%   see also testDiscreteCrossCorrelation, testInterpolateCrossCorrelation

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

    switch choice
        case 1
            y = sin(params(1)*x);
        case 2
            y = zeros(size(x));
            y(x>=0) = exp(-abs(params(2)*x(x>=0))).*sin(params(1)*x(x>=0));
%             y = exp(-abs(params(2)*x)).*sin(params(1)*x);
        case 3
            y = zeros(size(x));
            y(x>=0) = exp(-abs(params(1)*x(x>=0)));
%             y = exp(-abs(params(1)*x));
        otherwise
            error('Not implemented');
    end

return
