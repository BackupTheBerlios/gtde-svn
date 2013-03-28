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
    [XCPC, Energies] = Signals2XCPC(signals,experimentOptions.microphonePositions,experimentOptions.samplingFrequency);
    
    foundTDE = nan(experimentOptions.dimension,1);
    
    % Check if I want to use the constraint or not
    switch experimentOptions.methodOptions.type
        case 'bab'
            % Compute TDEs
            [solverOutput, experimentOptions.methodOptions] = gTDE_BAB(XCPC,...
                 experimentOptions.microphonePositions,...
                 1./experimentOptions.samplingFrequency,...
                 (experimentOptions.samplingFrequency*experimentOptions.initializationPositions),...
                 experimentOptions.methodOptions,Energies);
             % Feasible points
            indices = solverOutput.flag >= 0 & solverOutput.c >= 0;
            % Retrieve the optimal point
            if sum(indices)
%                 if sum(solverOutput.f == min(solverOutput.f(indices))) > 1, 
%                     fprintf('Voilu'); 
%                 end
                foundTDE = solverOutput.x(:,solverOutput.f == min(solverOutput.f(indices)))/experimentOptions.samplingFrequency;
                % If two points are equivalent, they'd probably be the
                % same, so take the first one. THIS IS HEURISTIC!!!!!
                if size(foundTDE,2) > 1
                    foundTDE = foundTDE(:,1);
                end
            end
        case 'eusipco'
            % Compute TDEs
            [solverOutput] = gTDE_EUSIPCO(XCPC,...
                 experimentOptions.microphonePositions,...
                 1./experimentOptions.samplingFrequency,...
                 (experimentOptions.samplingFrequency*experimentOptions.initializationPositions));
             % Feasible points
            indices = solverOutput.c >= 0;
            % Retrieve the optimal point
            if sum(indices)
                foundTDE = solverOutput.x(:,solverOutput.f == min(solverOutput.f(indices)))/experimentOptions.samplingFrequency;
            end
        case 'dip'
            % Compute TDEs
            [solverOutput, experimentOptions.methodOptions] = gTDE_DIP(XCPC,...
                 experimentOptions.microphonePositions,...
                 1./experimentOptions.samplingFrequency,...
                 (experimentOptions.samplingFrequency*experimentOptions.initializationPositions),...
                 experimentOptions.methodOptions,Energies);
             % Feasible points
            indices = solverOutput.flag' >= 0 & solverOutput.c >= 0;
            % Retrieve the optimal point
            if sum(indices)
                foundTDE = solverOutput.x(:,solverOutput.f == min(solverOutput.f(indices)))/experimentOptions.samplingFrequency;
            end
%         case 'tde'
%             % Compute TDEs
%             [solverOutput, experimentOptions.methodOptions] = ngTDE_DIP(PCCC,...
%                  experimentOptions.microphonePositions,...
%                  1./experimentOptions.samplingFrequency,...
%                  (experimentOptions.samplingFrequency*experimentOptions.initializationPositions),...
%                  experimentOptions.methodOptions);
        case {'init','truth'}
            % Compute the criterion for all the initial values
            [foundCriterion] = gTDECriterionMod(...
                experimentOptions.samplingFrequency*experimentOptions.initializationPositions,...
                XCPC,...
                experimentOptions.microphonePositions,...
                1/experimentOptions.samplingFrequency,...
                Energies);
            % SolverOutput building
            solverOutput.x = experimentOptions.initializationPositions;
            solverOutput.f = foundCriterion;
            solverOutput.c = TDEDiscriminant(solverOutput.x,experimentOptions.microphonePositions);
            % Constraints
            indices = solverOutput.c >= 0 & solverOutput.f >= 0;
            % Retrieve optimal point if any
            if sum(indices)
                foundTDE = experimentOptions.initializationPositions(:, solverOutput.f == min(solverOutput.f(indices)));
            end
        case 'bypairs'
            % Compute bypairs
            [solverOutput] = gTDEByPairs(XCPC,...
                 experimentOptions.microphonePositions,...
                 1./experimentOptions.samplingFrequency,...
                 experimentOptions.initializationPositions,...
                 experimentOptions.methodOptions);
             % Retrieve optimal poing
             foundTDE = solverOutput.xstar;
%         case 'sqplab'
%             [solverOutput] = gTDE_SQPLAB(PCCC,...
%                  experimentOptions.microphonePositions,...
%                  1./experimentOptions.samplingFrequency,...
%                  (experimentOptions.samplingFrequency*experimentOptions.initializationPositions));
        case 'bp2dip'
            experimentOptions.methodOptions.numMaxima = 3;
            % Compute bypairs
            [solverOutput] = gTDEByPairs(XCPC,...
                experimentOptions.microphonePositions,...
                1./experimentOptions.samplingFrequency,...
                experimentOptions.initializationPositions,...
                experimentOptions.methodOptions);
            % Get middle optimal points
            middleInitPoints = solverOutput.xstar;
            % Compute the new initialization points
            % That is all the possible combinations of all the values
            % output by the bypairs algorithms
            NPointsDimension = sum(~isnan(middleInitPoints),1);
            NInitPoints = prod(NPointsDimension);
            newInitializationPositions = zeros(experimentOptions.dimension,NInitPoints);
            % Initialize variables
            SubChunkSize = NInitPoints;
            NRepeat = 1;
            % Loop
            for d=1:experimentOptions.dimension,
                % SubChunkSize update
                SubChunkSize = SubChunkSize/NPointsDimension(d);
                Chunk = [];
                for ii=1:NPointsDimension(d),
                    Chunk = cat(2,Chunk,repmat(middleInitPoints(ii,d),1,SubChunkSize));
                end
                newInitializationPositions(d,:) = repmat(Chunk,1,NRepeat);
                % NRepeat update
                NRepeat = NRepeat * NPointsDimension(d);
            end
            % Launch DIP
            [solverOutput, experimentOptions.methodOptions] = gTDE_DIP(XCPC,...
                experimentOptions.microphonePositions,...
                1./experimentOptions.samplingFrequency,...
                (experimentOptions.samplingFrequency*newInitializationPositions),...
                experimentOptions.methodOptions,Energies);
            % Feasible points
            indices = solverOutput.flag' >= 0 & solverOutput.c >= 0;
            % Retrieve the optimal point
            if sum(indices)
                foundTDE = solverOutput.x(:,solverOutput.f == min(solverOutput.f(indices)))/experimentOptions.samplingFrequency;
            end
        case 'bab2dip'
             tolerance = TDEmax(experimentOptions.microphonePositions)*...
                experimentOptions.samplingFrequency/2^2;
            experimentOptions.methodOptions.tolerance = min(tolerance);
            experimentOptions.methodOptions.init = true;
            % Compute TDEs
            [solverOutput, experimentOptions.methodOptions] = gTDE_BAB(XCPC,...
                 experimentOptions.microphonePositions,...
                 1./experimentOptions.samplingFrequency,...
                 (experimentOptions.samplingFrequency*experimentOptions.initializationPositions),...
                 experimentOptions.methodOptions,Energies);
            % Launch DIP
            [solverOutput, experimentOptions.methodOptions] = gTDE_DIP(XCPC,...
                experimentOptions.microphonePositions,...
                1./experimentOptions.samplingFrequency,...
                solverOutput.x,...
                experimentOptions.methodOptions,Energies);
            % Feasible points
            indices = solverOutput.flag' >= 0 & solverOutput.c >= 0;
            % Retrieve the optimal point
            if sum(indices)
                foundTDE = solverOutput.x(:,solverOutput.f == min(solverOutput.f(indices)))/experimentOptions.samplingFrequency;
            end
        case 'fminunc'
            % Compute the minization
            [solverOutput, experimentOptions.methodOptions] = gTDE_FMINUNC(XCPC,...
                 experimentOptions.microphonePositions,...
                 1./experimentOptions.samplingFrequency,...
                 (experimentOptions.samplingFrequency*experimentOptions.initializationPositions),...
                 experimentOptions.methodOptions,Energies);
            % Feasibility
            indices = solverOutput.f >= 0 & solverOutput.flag > 0 & solverOutput.c >= 0;
            % Retrieve the optimum
            if sum(indices)
                foundTDE = solverOutput.x(:,solverOutput.f == min(solverOutput.f(indices)))/experimentOptions.samplingFrequency;
            end
    end
end
