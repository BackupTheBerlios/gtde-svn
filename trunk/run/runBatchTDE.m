%%% Clearing
close all; clc;

%%% Basic parameters
[~, hostname] = system('hostname');
if strcmp(hostname(1:end-1),'cassiopeia') || strcmp(hostname(1:end-1),'peirasmos')
    experimentOptions.rootFolder = '/local_scratch/alamedap/Software/gTDESSL/trunk/';
else
    cd('/scratch/cassiopeia');
    experimentOptions.rootFolder = '/scratch/cassiopeia/alamedap/Software/gTDESSL/trunk/';
end

%%% Add path
addpath(svngenpath(experimentOptions.rootFolder));

%%% Signal properties
% Type of signal
experimentOptions.dataUsed = 'simulated';
experimentOptions.field = 'far-field';
% experimentOptions.sourceMotion = false;
% Signal properties: sampling frequency and length
experimentOptions.samplingFrequency = 48000;
experimentOptions.length = 0.1;

%%% Sensor properties
% Type of sensor
% experimentOptions.sensorType = 'tetrahedronArray';
% Method to apply
% 'gtde' with geometric constraint
% 'tde' withouth geometric constraint
% 'init' just the minimum over the initialization set, no local
%           minimization
% 'bypairs' estimating by pairs
experimentOptions.method = 'bypairs';
% SNR values (this could also be a vector)
experimentOptions.snrValues = 5;

%%% ISM
experimentOptions.ism = true;

% If ISM
experimentOptions.microphoneOffset = [2.25 1.25 1.25];
experimentOptions.room = [3  4  2.5];
experimentOptions.absorption_weights =  [1 1 1 1 1 1];
experimentOptions.T60 = 0.2;
% ISM structure
experimentOptions.ismFolder= 'ism_rir/';
experimentOptions.generateISMData = false;

%%% Wav folder
experimentOptions.wavFolder= '/scratch/pictor/deleforg/the_CASA_REDMINE/emitted_sounds/';

%%% Sound source properties
% Dimension of the ambient space
experimentOptions.dimension = 3;
% Source position options
experimentOptions.sourcePositionOptions.coordinateSystem = 'cartesian';
experimentOptions.sourcePositionOptions.bounds = cell(3,1);
experimentOptions.sourcePositionOptions.bounds{1} = [experimentOptions.room(1)/20 19*experimentOptions.room(1)/20];
experimentOptions.sourcePositionOptions.bounds{2} = [experimentOptions.room(2)/20 19*experimentOptions.room(2)/20];
experimentOptions.sourcePositionOptions.bounds{3} = [experimentOptions.room(3)/20 19*experimentOptions.room(3)/20];
experimentOptions.sourcePositionOptions.numberOfIntervals = [4,4,4];
% Source positions precomputed or generated
experimentOptions.generateSourcePositions = false;
if ~experimentOptions.generateSourcePositions
%     tmp = load(strcat(experimentOptions.rootFolder,'data/sourcePositions.mat'));
    tmp = load(strcat(experimentOptions.rootFolder,'data/sourcePositionsCenter.mat'));
    experimentOptions.sourcePositions = tmp.sourcePositions;
%     experimentOptions.sourcePositions = experimentOptions.sourcePositions(1,:);
end
experimentOptions.tiltIndices = [1 4 7 10 13];
experimentOptions.panIndices = [1 7 13 19 25 31];

%%% Run the batch
batchTDEExperiments(experimentOptions);
