function noisy = AddWhiteNoise(noiseless,SNR)

%Add white noise to a signal
%
% USAGE: noisy = addWhiteNoise(noiseless,SNR)
%
% PARAMETERS:
%   noiseless ~ noise less signal
%   SNR ~ signal-to-noise ratio (in dB!)
% 
% RETURN VALUE:
%   noisy ~ noisy signal
%
% DESCRIPTION: Adds gaussian white noise to a signal
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

    % Compute signal's energy
    signalEnergy = sum(noiseless.^2)/numel(noiseless);
    % Compute noise variance
    noiseVariance = 10^( log10(signalEnergy) - SNR/10 );
    % Add gaussian noise
    noisy = noiseless + sqrt(noiseVariance)*randn(size(noiseless));
end