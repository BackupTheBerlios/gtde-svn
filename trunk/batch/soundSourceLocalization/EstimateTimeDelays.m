function [foundTDE, solverOutput] = EstimateTimeDelays(signals,experimentOptions)

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
    %%% From signals to Polynomial Coefficient Cross Correlation
%     PCCC = Signals2PCCC(signals,experimentOptions.microphonePositions,experimentOptions.samplingFrequency);
    XCPC = Signals2XCPC(signals,experimentOptions.microphonePositions,experimentOptions.samplingFrequency);
    
    % Check if I want to use the constraint or not
    switch experimentOptions.methodOptions.type
        case 'dip'
            % Compute TDEs
            [solverOutput, experimentOptions.methodOptions] = gTDE_DIP(XCPC,...
                 experimentOptions.microphonePositions,...
                 1./experimentOptions.samplingFrequency,...
                 (experimentOptions.samplingFrequency*experimentOptions.initializationPositions),...
                 experimentOptions.methodOptions);
        case 'tde'
            % Compute TDEs
            [solverOutput, experimentOptions.methodOptions] = ngTDE_DIP(PCCC,...
                 experimentOptions.microphonePositions,...
                 1./experimentOptions.samplingFrequency,...
                 (experimentOptions.samplingFrequency*experimentOptions.initializationPositions),...
                 experimentOptions.methodOptions);
        case {'init','truth'}
            % Compute the criterion for all the initial values
            [foundCriterion] = gTDECriterion(...
                experimentOptions.samplingFrequency*experimentOptions.initializationPositions,...
                PCCC,...
                experimentOptions.microphonePositions,...
                1/experimentOptions.samplingFrequency);
            % SolverOutput building
            solverOutput.x = experimentOptions.initializationPositions;
            solverOutput.f = foundCriterion;
        case 'bypairs'
            [solverOutput] = gTDEByPairs(PCCC,...
                 experimentOptions.microphonePositions,...
                 1./experimentOptions.samplingFrequency,...
                 experimentOptions.initializationPositions);
        case 'sqplab'
            [solverOutput] = gTDE_SQPLAB(PCCC,...
                 experimentOptions.microphonePositions,...
                 1./experimentOptions.samplingFrequency,...
                 (experimentOptions.samplingFrequency*experimentOptions.initializationPositions));
    end
    
    if ~strcmp(experimentOptions.methodOptions.type,'bypairs')
        foundTDE = getOptimalTDE(solverOutput);
    else
        foundTDE = solverOutput.xstar;
    end
end

function xstar = getOptimalTDE(solverOutput)
    % Get the optimal x
    [~,index] = min(solverOutput.f);
    xstar = solverOutput.x(:,index);
end
