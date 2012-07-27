function batchTDEExperiments(experimentOptions)

%Batch the experiments on geometrically constrained time delay estimation
%
% USAGE: batchTDEExperiments(experimentOptions)
%
% PARAMETERS:
%   experimentOptions ~ options of the experiment
%
% TODO: describe the options
%
%
%   see also 

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

    if ~isfield(experimentOptions,'rootDir')
        experimentOptions.rootDir = '/local_scratch/alamedap/Software/gTDESSL/trunk/';
    end

    %%% General variables
    % Check options
    checkExperimentOptions(experimentOptions);
    
    % Experiment comments
%    comment = inputdlg('Please enter the experiment comments.');
%    experimentOptions.comment = comment{1};

    % Get Source positions
    if experimentOptions.generateSourcePositions
        experimentOptions.sourcePositions = getSourcePositions(experimentOptions);
    else
        if ~isfield(experimentOptions,'sourcePositions')
            error('If you do not want the source positions to be generated, you have to specify them.');
        end
    end
    nSourcePositions = size(experimentOptions.sourcePositions,1);
    
    % Get microphones positions
    experimentOptions.microphonesPositions = getMicrophonesPositions(experimentOptions);
    if ~strcmp(experimentOptions.dataUsed,'real'),
        experimentOptions.microphonesPositions = experimentOptions.microphonesPositions + ...
            repmat(experimentOptions.microphoneOffset,size(experimentOptions.microphonesPositions,1),1);
        experimentOptions.ismStruct.mic_pos = experimentOptions.microphonesPositions;
        % Generate ISM data (RIR)
        if experimentOptions.generateISMData
            generateISMData(experimentOptions);
            return
        end
    end
    
    % Generate Initialization positions
    experimentOptions.initializationPositions = getInitializationPositions(experimentOptions);
    
    % Get signals
    experimentOptions.signals = getSignals(experimentOptions);
    nSignals = numel(experimentOptions.signals);
    if strcmp(experimentOptions.dataUsed,'simulated')
        nSignals = 3;
    elseif strcmp(experimentOptions.dataUsed,'real')
        nSignals = 1;
    end 

    % SNR values
    nSNRValues = length(experimentOptions.snrValues);
    
    % T60 values
    nT60 = length(experimentOptions.T60);

    % Results
    foundTDEs = cell(nSignals,...
                      nSourcePositions,...
                      nSNRValues,...
                      nT60);
    
    % Control the time
    donePart = 0;
    tStart = clock;
    
   
    % For each position
    for sPosition = 1:nSourcePositions,
        fprintf('Position %d/%d:\n',sPosition,nSourcePositions);
        % For each signal
        for sSignal = 1:nSignals
            fprintf('    Signal %d/%d   <',sSignal,nSignals);
            % For each SNR value
            for sSNR = 1:nSNRValues,
                % For each T60 value,
                for sT60 = 1:nT60,
                    % Get the signals
                    [signals fs] = getSignalRealization(experimentOptions,sSignal,sPosition,sSNR,sT60);
                    % Reassign the sampling frequency
                    experimentOptions.samplingFrequency = fs;
                    
                    % Allocate results
                    foundTDEs{sSignal,sPosition,sSNR,sT60} = zeros(numel(signals),experimentOptions.dimension);
                    % Loop for each partial signal
                    for subs = 1:numel(signals),
                        % Get partial signal
                        signalPiece = signals{subs};

%                         partStart = tic;
                        % Apply the method
                        foundTDE = applyTDE(signalPiece,experimentOptions);
%                         partElapsed = toc(partStart);
%                         fprintf('%1.1f ',partElapsed);

                        % Save the result
                        foundTDEs{sSignal,sPosition,sSNR,sT60}(subs,:) = foundTDE; 
                    end
                    % Notify :)
                    donePart = donePart + 100/(nSignals*nSourcePositions*nSNRValues*nT60);
                    fprintf('.');
                end % sT60
            end % sSNR
            tEnd = clock;
            eTime = etime(tEnd,tStart);
            rTime = eTime*(100-donePart)/donePart;
            fprintf('>   [Done by: %s]\n',datestr(now+(rTime/86400),0));
        end % sPosition
    end % sSignal
    
    % Save
    % NOTE, we modified the filename. We assume that the computations on
    % the same type of sensor, the same data used, the same t60 and the
    % same snrValue are not running at the same time.
    if strcmp(experimentOptions.dataUsed,'real')
        fileName = strcat(experimentOptions.rootFolder,...
                          '../results/',...
                          experimentOptions.method,'_',...
                          experimentOptions.dataUsed,'_',...
                          experimentOptions.sensorType,'_',...
                          datestr(clock,'yyyy-mm-dd-HH-MM-SS'),'_',...
                          experimentOptions.field,'.mat');
    else
        fileName = strcat(experimentOptions.rootFolder,...
                          '../results/',...
                          experimentOptions.method,'_',...
                          experimentOptions.dataUsed,'_',...
                          experimentOptions.sensorType,'_',...
                          datestr(clock,'yyyy-mm-dd-HH-MM-SS'),'_',...
                          num2str(experimentOptions.T60(1)),'_',...
                          num2str(experimentOptions.snrValues(1)),'.mat');
    end
    save(fileName,'experimentOptions','foundTDEs');
    
end % function

function checkExperimentOptions(experimentOptions)
    
    % Check that there is the kind of data used and the type of sensor
    mandatoryFields = {'dataUsed','sensorType','dimension','method',...
        'sourcePositionOptions','snrValues','ism','rootFolder','room',...
        'microphoneOffset','sourceMotion','length','generateISMData','generateSourcePositions'};
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
    
    sensor = {'linearArray','tetrahedronArray'};
    % Check the sensor specified
    sensorFlag = 0;
    for s = 1:numel(sensor),
        sensorFlag = sensorFlag + strcmp(experimentOptions.sensorType,sensor{s});
    end
    if sensorFlag == 0
        error('Sensor type specified not known.');
    end    
    
    % The combination of real data and linear array is forbidden.
    if strcmp(experimentOptions.dataUsed,'real') && strcmp(experimentOptions.sensorType,'linearArray')
        error('We do not have real data for linear array.');
    end
    if strcmp(experimentOptions.dataUsed,'real') && experimentOptions.ism
        fprintf('Warning: No needing for simulate the ISM with real data.');
    end
    if strcmp(experimentOptions.sensorType,'linearArray') && strcmp(experimentOptions.constraint,true)
        error('The use of the geometric constraint in the linear array case makes no sense.');
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

function sourcePositions = getSourcePositions(experimentOptions)
    % Initializing to the empty set
    sourcePositions = [];
    % Choose depending on the data used for the experiment
    switch experimentOptions.dataUsed
        case 'real'
            sourcePositions = GeneratePositionsCasaRedmine(experimentOptions.panIndices,experimentOptions.tiltIndices,experimentOptions.field);
        case 'simulated'
            sourcePositions = GeneratePositions(experimentOptions.dimension,experimentOptions.sourcePositionOptions);
        case 'synthetic'
            sourcePositions = GeneratePositions(experimentOptions.dimension,experimentOptions.sourcePositionOptions);
    end
end

function microphonesPositions = getMicrophonesPositions(experimentOptions)
    % Initialization
    microphonesPositions = [];
    
    % Check different cases
    if ~strcmp(experimentOptions.dataUsed,'real')
        % Dimension is either 2 or 3
        if experimentOptions.dimension == 2
            % Scale parameters
            a = 0.1;
            b = -a/2;
            c = a*sqrt(9/12);
            % Microphones in a regular tetrahedron
            microphonesPositions = ...
                [ a   0 ;...
                  b   c ;...
                  b  -c];
        else
            % Scale parameters
            a = 0.1;
            c = a / sqrt(2);
            % Microphones in a regular tetrahedron
            microphonesPositions = ...
                [ a   0   -c ;...
                 -a   0   -c ;...
                  0   a    c ;...
                  0  -a    c];
            
        end
    elseif strcmp(experimentOptions.dataUsed,'real')
        %%% Importdata
        data = importdata(strcat('/scratch/Andromeda/amat/WorkingDevLocal/CASA_RED_MINE_19122011/snapshots/',experimentOptions.field,'/points_3d_new_coord_syst.txt'));
        % Microphones' positions and speaker's position
        p = data(:,2:end);
        % Allocate microphones positions
        microphonesPositions = zeros(4,3);
        microphonesPositions(1,:) = p(3,:);
        microphonesPositions(2,:) = p(4,:);
        microphonesPositions(3,:) = p(2,:);
        microphonesPositions(4,:) = p(1,:);        
    end
end

function initializationPositions = getInitializationPositions(experimentOptions)

    % Dimension
    initializationPositionOptions.dimension = size(experimentOptions.microphonesPositions,1)-1;
    % Coordinate system
    initializationPositionOptions.coordinateSystem = 'cartesian';
    % Bounds
    maxTDEs = TDEmax(experimentOptions.microphonesPositions);
    initializationPositionOptions.bounds = cell(initializationPositionOptions.dimension,1);
    for d = 1:initializationPositionOptions.dimension
        initializationPositionOptions.bounds{d} = [-maxTDEs(d),maxTDEs(d)];
    end
    % Number of intervals
    initializationPositionOptions.numberOfIntervals = 15;
    % Generate positions
    initializationPositions = GeneratePositions(initializationPositionOptions.dimension,initializationPositionOptions);

    % Modify the initial positions depending on the method used
    switch experimentOptions.method
        case 'gtde'
            Constraint = TDEDiscriminant(initializationPositions',experimentOptions.microphonesPositions);
            initializationPositions = initializationPositions(Constraint>0,:);
        case 'tde'
        case 'init'
        case 'bypairs'
            auxInit = initializationPositions;
            % Reallocate them
            initializationPositions = cell(1,experimentOptions.dimension);
            % Chose the unique values for each column
            for d = 1:experimentOptions.dimension,
                initializationPositions{d} = unique(auxInit(:,d));
            end
    end
end

function signals = getSignals(experimentOptions)
    % Initialize
    signals = [];
    % Depending on the data
    % In case of synthetic data, they are functions
    if strcmp(experimentOptions.dataUsed,'synthetic')
        NSignals = 1;
        signals = cell(NSignals,1);
        % Synthetic signals
        for ns = 1:NSignals,
            [sinusCoefficients sinusPhases sinusFrequencies] = GenerateRandomSpectrumParameters(1/experimentOptions.samplingFrequency,experimentOptions.length);
            signals{ns} = @(x) GenerateRandomSpectrumSignal(x,sinusCoefficients,sinusPhases,sinusFrequencies);
        end  
    elseif strcmp(experimentOptions.dataUsed,'simulated')
        % Retrieve folder
        [~, signalNames] = system(['find ' experimentOptions.wavFolder ' | grep .wav']);
        % Cut in names
        indices = strfind(signalNames,'.wav');
        % Number of signals, and cell declaration
        NSignals = numel(indices);
        signals = cell(NSignals,1);
        % Complete to ease the loop
        indices = cat(2,-4,indices);
        for kk = 1:numel(indices)-1,
            signals{kk} = signalNames(indices(kk)+5:indices(kk+1)+3);
        end
    elseif strcmp(experimentOptions.dataUsed,'real')
        signals = strcat('/scratch/pictor/deleforg/the_CASA_REDMINE/audio_recordings/',experimentOptions.field,'/sound12/Recorded/');
    else
        error('Not ready yet!');
    end
end

function [signals fs] = getSignalRealization(experimentOptions,sSignal,sPosition,sSNR,sT60)
    % Depending on the data
    % In case of synthetic data, they are functions
    if strcmp(experimentOptions.dataUsed,'synthetic')
        fs = experimentOptions.samplingFrequency;
        % ISM case
        if experimentOptions.ism
            % Generate signal
            signal = GenerateDiscreteSignals(1./experimentOptions.samplingFrequency,...
                                             experimentOptions.length,...
                                             experimentOptions.signals{sSignal});
            % ISM File Name
            ISMFileName = getISMFileName(experimentOptions,...
                experimentOptions.T60(sT60),...
                experimentOptions.sourcePositions(sPosition,:),...
                fs);
            if ~(exist(ISMFileName,'file') == 2)
                error(strcat('ISM filter does not exist (',ISMFileName,').'));
            end
            % Apply RIR filter
            partialSignals = ISM_AudioData(ISMFileName,signal,'SilentFlag',1);
            partialSignals = partialSignals'; 
        % Non ISM case
        else
            % Generate the signals
            partialSignals = GenerateDiscreteSignals(...
                experimentOptions.sourcePositions(sPosition,:),...
                experimentOptions.microphonesPositions,...
                1./experimentOptions.samplingFrequency,...
                experimentOptions.length,...
                experimentOptions.signals{sSignal}...
                );
        end
    elseif strcmp(experimentOptions.dataUsed,'simulated')
        %%%%%%%%%%%%% WARNING %%%%%%%%%%%%%%%
        % THE NON REVERBERANT CASE IS T60 = 0
        %%%%%%%%%%%%% WARNING %%%%%%%%%%%%%%%
        % Load the signal
        switch sSignal,
            case 1,
                [signal fs] = wavread(experimentOptions.signals{2});
            case 2,
                [signal fs] = wavread(experimentOptions.signals{9});
            case 3,
                index = randi(numel(experimentOptions.signals),1);
                [signal fs] = wavread(experimentOptions.signals{index});
            otherwise
                error('Signal should be either 1 (caugh), 2 (female speech) or 3 (random sound).');
        end
        % ISM File Name
        ISMFileName = getISMFileName(experimentOptions,...
            experimentOptions.T60(sT60),...
            experimentOptions.sourcePositions(sPosition,:),...
            fs);
        if ~(exist(ISMFileName,'file') == 2)
            error(strcat('ISM filter does not exist (',ISMFileName,').'));
        end
        % Apply RIR filter
        partialSignals = ISM_AudioData(ISMFileName,signal,'SilentFlag',1);
        partialSignals = partialSignals'; 
    elseif strcmp(experimentOptions.dataUsed,'real'),
        % Compute the recording index, from sPos
%         numTilts = numel(experimentOptions.tiltIndices);
%         numPans = numel(experimentOptions.panIndices);
        sPos = sPosition-1;
        tiltIndex = 3*floor(sPos/6)+1;
        panIndex = 6*(sPos-6*fix(sPos/6))+1;
        wavIndex = 37*(tiltIndex-1)+panIndex;
%         fprintf('%s\n',strcat(experimentOptions.signals,num2str(wavIndex),'.wav'));
        [partialSignals fs] = wavread(strcat(experimentOptions.signals,num2str(wavIndex),'.wav'));
        partialSignals = partialSignals';
    end
        
    %%% Cut the signals
    pieceLength = experimentOptions.length*fs; 
    NPieces = floor(size(partialSignals,2)/pieceLength);
    signals = cell(NPieces,1);

    for np = 1:NPieces,
        signals{np} = partialSignals( :, ((np-1)*pieceLength+1):(np*pieceLength) );
        if ~strcmp(experimentOptions.dataUsed,'real'),
            %%% Add noise
            for channel = 1:size(partialSignals,1),
                signals{np}(channel,:) = addWhiteNoise(signals{np}(channel,:),experimentOptions.snrValues(sSNR));
            end
        end
    end
end

function noisy = addWhiteNoise(noiseless,SNR)
    % Compute signal's energy
    signalEnergy = sum(noiseless.^2)/numel(noiseless);
    % Compute noise variance
    noiseVariance = 10^( log10(signalEnergy) - SNR/10 );
    % Add gaussian noise
    noisy = noiseless + sqrt(noiseVariance)*randn(size(noiseless));
end

function foundTDE = applyTDE(signals,experimentOptions)
    %%% Compute the interpolation coefficients
    % Number of samples
    NMics = size(signals,1);
    NSamples = size(signals,2);
    % Sampling times
    samplingTimes = 0:1/experimentOptions.samplingFrequency:(NSamples-1)/experimentOptions.samplingFrequency;
    % Allocate
    PC = cell(NMics,1);
    % Compute
    for ss = 1:NMics,
        PC{ss} = PolynomialInterpolationCoefficients(signals(ss,:),samplingTimes);
    end
    % Compute the Cross-correlation of the interpolation coefficients
    % Allocate
    PCCC = cell(NMics);
    % Compute
    for mic1 = 1:NMics,
        for mic2 = mic1:NMics,
            PCCC{mic1,mic2} = PolynomialCoefficientsCrossCorrelation(PC{mic1},PC{mic2});%,maxLAG);
        end
    end
    % Check if I want to use the constraint or not
    switch experimentOptions.method
        case 'gtde'
            % Compute TDEs
            [foundTDE foundCriterion foundConstraint] = gTDE_parallel(PCCC,...
                 experimentOptions.microphonesPositions,...
                 1./experimentOptions.samplingFrequency,...
                 (experimentOptions.samplingFrequency*experimentOptions.initializationPositions)');
            % Select the minimum criterion of those who satisfy the constraint
            foundCriterion( foundConstraint > 0 ) = max(foundCriterion);
            [~,minIndex] = min(foundCriterion);
            foundTDE = foundTDE(:,minIndex)/experimentOptions.samplingFrequency;
        case 'tde'
            % Compute TDEs
            [foundTDE foundCriterion] = ngTDE_parallel(PCCC,...
                 experimentOptions.microphonesPositions,...
                 1./experimentOptions.samplingFrequency,...
                 (experimentOptions.samplingFrequency*experimentOptions.initializationPositions)');
            % Select the minimum criterion
            [~,minIndex] = min(foundCriterion);
            foundTDE = foundTDE(:,minIndex)/experimentOptions.samplingFrequency;
        case 'init'
            % Compute the criterion for all the initial values
            [foundCriterion] = gTDECriterion(...
                experimentOptions.samplingFrequency*experimentOptions.initializationPositions',...
                PCCC,...
                experimentOptions.microphonesPositions,...
                1/experimentOptions.samplingFrequency);
            % Select the minimum criterion
            [~,minIndex] = min(foundCriterion);
            foundTDE = experimentOptions.initializationPositions(minIndex,:);
        case 'bypairs'
            [foundTDE] = gTDEByPairs(PCCC,...
                 experimentOptions.microphonesPositions,...
                 1./experimentOptions.samplingFrequency,...
                 experimentOptions.initializationPositions);
    end
end

function generateISMData(experimentOptions)
    % Sampling Frequencies
    samplingFrequencies = [16000 44100 48000];
    % Static source case
    if ~experimentOptions.sourceMotion
        donePart = 0;
        % Done part done
        % For each value of T60
        for nt = 1:numel(experimentOptions.T60)
            t60 = experimentOptions.T60(nt);
            % For each sampling frequency
            for nSF = 1:numel(samplingFrequencies),
            SF = samplingFrequencies(nSF);
                % For each sound source positions
                for np = 1:size(experimentOptions.sourcePositions,1)
                    % Extract position
                    position = experimentOptions.sourcePositions(np,:);
                    % ISM file name
                    ISMFileName = getISMFileName(experimentOptions,t60,position,SF);
                    % If the file does not exist, create it
                    if ~(exist(ISMFileName,'file') == 2)
                        ISMStruct = gTDE_ISM_setup();
                        % Modify it
                        ISMStruct.Fs = SF;
                        ISMStruct.room = experimentOptions.room;
                        ISMStruct.mic_pos = experimentOptions.microphonesPositions;
                        ISMStruct.src_traj = position;
                        ISMStruct.T60 = t60;
                        ISMStruct.abs_weights = experimentOptions.absorption_weights;
                        % Create the RIR filter
                        ISM_RIR_bank(ISMStruct,ISMFileName);
                    end
                    donePart = donePart + 100/(numel(experimentOptions.T60)*size(experimentOptions.sourcePositions,1)*numel(samplingFrequencies));
                    fprintf('T60: %1.1f %d [x= %1.3f, y= %1.3f, z= %1.3f] (%1.1f %%)\n',t60,SF,position(1),position(2),position(3),donePart);

                end
            end
        end
    else
        error('Not ready for dynamic sources.');
    end
end

function ISMFileName = getISMFileName(experimentOptions,t60,position,SF)
    % Build ISM file name
    ISMFileName = strcat(experimentOptions.rootFolder,...
                     '../',...
                     experimentOptions.ismFolder,...
                     'ism_rir_',...
                     num2str(t60),'_',...
                     num2str(SF),'_',...
                     num2str(position(1)),'_',...
                     num2str(position(2)),'_',...
                     num2str(position(3)),'.mat');
end
