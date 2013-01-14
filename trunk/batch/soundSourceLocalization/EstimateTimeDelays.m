function foundTDE = EstimateTimeDelays(signals,experimentOptions)

%Perform the time delay estimation for a particular discrete signal
%
% USAGE: TDE = EstimateTimeDelays(signals,experimentOptions)
%
% PARAMETERS:
%   signals ~ discrete signal received at the microphones
%   experimentOptions ~ options of the experiment
% 
% RETURN VALUE:
%   TDE ~ time delay estimates
%
% DESCRIPTION: Estimates the time delay to the discrete signals following
%   the options specified in the structure.
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

    %%% Compute the interpolation coefficients
    % Number of samples
    NMics = size(signals,1);
    NSamples = size(signals,2);
    % Compute maxlag
    maxLag = TDEmax(experimentOptions.microphonePositions);
    maxLag = ceil(experimentOptions.samplingFrequency*max(maxLag));
%     % Sampling times
%     samplingTimes = 0:1/experimentOptions.samplingFrequency:(NSamples-1)/experimentOptions.samplingFrequency;
    % Allocate
    PC = cell(NMics,1);
    % Compute
    for ss = 1:NMics,
        PC{ss} = PolynomialInterpolationCoefficients(signals(ss,:),1/experimentOptions.samplingFrequency);
    end
    % Compute the Cross-correlation of the interpolation coefficients
    % Allocate
    PCCC = cell(NMics);
    % Compute
    for mic1 = 1:NMics,
        for mic2 = mic1:NMics,
            PCCC{mic1,mic2} = PolynomialCoefficientsCrossCorrelation(PC{mic1},PC{mic2},maxLag);
        end
    end
    % Check if I want to use the constraint or not
    switch experimentOptions.methodOptions.type
        case 'gtde'
            % Compute TDEs
            [foundTDE, foundCriterion, foundConstraint] = gTDE_parallel(PCCC,...
                 experimentOptions.microphonePositions,...
                 1./experimentOptions.samplingFrequency,...
                 (experimentOptions.samplingFrequency*experimentOptions.initializationPositions)');
            % Select the minimum criterion of those who satisfy the constraint
            foundCriterion( foundConstraint > 0 ) = max(foundCriterion);
            [~,minIndex] = min(foundCriterion);
            foundTDE = foundTDE(:,minIndex)/experimentOptions.samplingFrequency;
        case 'tde'
            % Compute TDEs
            [foundTDE, foundCriterion] = ngTDE_parallel(PCCC,...
                 experimentOptions.microphonePositions,...
                 1./experimentOptions.samplingFrequency,...
                 (experimentOptions.samplingFrequency*experimentOptions.initializationPositions)');
            % Select the minimum criterion
            [~,minIndex] = min(foundCriterion);
            foundTDE = foundTDE(:,minIndex)/experimentOptions.samplingFrequency;
        case {'init','truth'}
            % Compute the criterion for all the initial values
            [foundCriterion] = gTDECriterion(...
                experimentOptions.samplingFrequency*experimentOptions.initializationPositions',...
                PCCC,...
                experimentOptions.microphonePositions,...
                1/experimentOptions.samplingFrequency);
            % Select the minimum criterion
            [~,minIndex] = min(foundCriterion);
            foundTDE = experimentOptions.initializationPositions(minIndex,:);
        case 'bypairs'
            [foundTDE] = gTDEByPairs(PCCC,...
                 experimentOptions.microphonePositions,...
                 1./experimentOptions.samplingFrequency,...
                 experimentOptions.initializationPositions);
        case 'sqplab'
            [foundTDE] = gTDE_SQPLAB(PCCC,...
                 experimentOptions.microphonePositions,...
                 1./experimentOptions.samplingFrequency,...
                 (experimentOptions.samplingFrequency*experimentOptions.initializationPositions)');
    end
end
