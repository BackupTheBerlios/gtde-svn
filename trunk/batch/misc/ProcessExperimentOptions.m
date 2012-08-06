function experimentOptions = ProcessExperimentOptions(experimentOptions)

%Check the options of the experiment
%
% USAGE: experimentOptions = ProcessExperimentOptions(experimentOptions)
%
% PARAMETERS:
%   experimentOptions ~ structure containing the options of the experiment
%
% DESCRIPTION: the options are checked, if something fails an error is
%   reported.
%
%   The structure MUST contain the following fields:
%
%     globalOptions ~ stucture with the global options:
%       dimension -- dimension of the experiment (either 2 or 3)
%       rootFolder -- rootFolder of the code, all the results,
%         precomputations, wav signals etc, are supposed to be in there (in
%         some particular, further on specified, folder structure).   
%
%     dataOptions ~ structure describing the data options:
%       type -- synthetic, simulated or real
%       wavFolder -- in case of simulated data, we need the emitted signals folder
%       snrValues -- in case of simulated data, we need the SNR to add noise
%
%     microphonePositionOptions ~ structure with the microphones' options:
%       scale -- scale of the microphone array (TODO, see what?)
%
%     method ~ the method used to estimate the time delays: 'gtde' for
%       geometrically constrained tde, 'tde' for unconstrained tde, 'init'
%       for unconstrained tde without local minimization and 'bypairs' for
%       pairwise estimation.
%
%     sourcePositionOptions ~ structure: if you use simulated or synthetic 
%       data, check GeneratePosition to fill it.
%
%     snrValues ~ the values for the snr
%
%     ism ~ boolean: true if you want to use the ISM model [1,2]. In that
%       case you need the fields: absorption_weights, ISMFolder and T60.
%       The first corresponds to the relative absorption weights, the
%       second to the folder in which the ISM filetrs are precomputed (see
%       the generateISMData option) and T60 is the parameter controlling
%       the level of reverberation.
%
%     rootFolder ~ 
%
%     room ~ this specifies the room field of the ISM structure (see [1])/
%
%     microphoneOffset ~ specifies the offset of the microphones'
%       positions. This vector will be added to the positions
%
%     length ~ 
%
%     generateISMData ~
%
%     generateSourcePositions ~
% 
%   REFERENCES:
%     [1] E. A. Lehmann. Matlab code for image-source model in room acoustics.
%         http://www.eric-lehmann.com/ism code.html.
%     [2] E. A. Lehmann and A. M. Johansson. Prediction of energy decay in 
%         room impulse responses simulated with an image-source model. The 
%         Journal of the Acoustical Society of, 124(1):269–277, 2008.
%
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

    % Check that there is the kind of data used and the type of sensor
    mandatoryFields = {'source','dimension','method',...
        'sourcePositionOptions','snrValues','ism','rootFolder','room',...
        'microphoneOffset','length','generateISMData','generateSourcePositions'};
    for f = 1:numel(mandatoryFields),
        if ~isfield(experimentOptions,mandatoryFields{f})
            error(['Field ' mandatoryFields{f} ' is mandatory.']);
        end
    end

    data = {'synthetic','simulated','real'};
    % Check the data specified
    dataFlag = 0;
    for d = 1:numel(data)
        dataFlag = dataFlag + strcmp(experimentOptions.dataUsed,data{d});
    end
    if dataFlag == 0
        error('Data specified not known.');
    end
    
    method =  {'gtde','tde','init','bypairs'};
    % Check the method specified
    methodFlag = 0;
    for d = 1:numel(method)
        methodFlag = methodFlag + strcmp(experimentOptions.method,method{d});
    end
    if methodFlag == 0
        error('Method specified not known.');
    end
    
    % Particular fields needed
    if strcmp(experimentOptions.dataUsed,'synthetic')
        mandatoryFields = {'samplingFrequency'};
        for f = 1:numel(mandatoryFields),
            if ~isfield(experimentOptions,mandatoryFields{f})
                error(['Field ' mandatoryFields{f} ' is mandatory when using synthetic data.']);
            end
        end
    end
    if strcmp(experimentOptions.dataUsed,'simulated')
        mandatoryFields = {'wavFolder'};
        for f = 1:numel(mandatoryFields),
            if ~isfield(experimentOptions,mandatoryFields{f})
                error(['Field ' mandatoryFields{f} ' is mandatory when using simulated data.']);
            end
        end
    end
    if strcmp(experimentOptions.dataUsed,'real')
        mandatoryFields = {'field','panIndices','tiltIndices'};
        for f = 1:numel(mandatoryFields),
            if ~isfield(experimentOptions,mandatoryFields{f})
                error(['Field ' mandatoryFields{f} ' is mandatory when using real data.']);
            end
        end
        if ~experimentOptions.generateSourcePositions
            error('In the real data experiments the source positions are always generated (read from file).');
        end
    end
    
    if strcmp(experimentOptions.dataUsed,'real') && experimentOptions.ism
        fprintf('Warning: No needing for simulate the ISM with real data.');
    end
    
    % Check the dimension
    if experimentOptions.dimension ~= 2 && experimentOptions.dimension ~= 3
        error('Dimension should be 2 or 3.');
    end
    if experimentOptions.dimension == 2 && experimentOptions.ism
        error('Not simulating the room in the 2D case.');
    end
    if experimentOptions.dimension == 2 && strcmp(experimentOptions.dataUsed,'real'),
        error('We do not have real data with dimension 2.');
    end
    
    % ISM
    if experimentOptions.ism
        if ~isfield(experimentOptions,'absorption_weights')
            error('If you want to use ISM, you need to specify the absortion weights.');
        end
        if ~isfield(experimentOptions,'ismFolder')
            error('If you want to use ISM, you need to specify an ISM folder.');
        end
        if ~isfield(experimentOptions,'T60')
            error('If you want to use ISM, you need to specify value(s) for T60.');
        end
    end

end