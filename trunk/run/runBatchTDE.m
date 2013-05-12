
%Running TDE batch example
%
% USAGE: runBatchTDE
%
% DESCRIPTION:
%   This is an example on how to set the options for the
%   batchTDEExperiments script. A lot of parameters and comments are made
%   such that it is somehow intuitive.
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

%%% Clearing
close all; clc; clear;

%%% Used variables
room = [4  4  4];

%%% Global options
% Root folder
[~, hostname] = system('hostname');
if strcmp(hostname(1:end-1),'cassiopeia') || strcmp(hostname(1:end-1),'peirasmos')
    theRootFolder = '/local_scratch/alamedap/';
else
    cd('/scratch/cassiopeia');
    theRootFolder = '/scratch/cassiopeia/alamedap/';
end
experimentOptions.rootFolder = strcat(theRootFolder,'Software/gtde/trunk/');

% Dimension of the ambient space
experimentOptions.dimension = 3;
% Use of ism
experimentOptions.ism = true;
% Add the root folder to the path (the software)
addpath(svngenpath(experimentOptions.rootFolder));

%%% Source position options
% Source position options
sourcePositionOptions.file = strcat(theRootFolder,'audispace/gTDEJournalData/sourceGT.txt');
sourcePositionOptions.coordinateSystem = 'cartesian';
sourcePositionOptions.dimension = experimentOptions.dimension;
sourcePositionOptions.bounds = cell(3,1);
sourcePositionOptions.bounds{1} = [room(1)/20 19*room(1)/20];
sourcePositionOptions.bounds{2} = [room(2)/20 19*room(2)/20];
sourcePositionOptions.bounds{3} = [room(3)/20 19*room(3)/20];
sourcePositionOptions.numberOfIntervals = [4,4,4];
% Save it
experimentOptions.sourcePositionOptions = sourcePositionOptions;

%%% Data options
% Type of signal
dataOptions.type = 'real';
% Length of the signal cuts
dataOptions.cutLength = 0.1;
% Depending
if strcmp(dataOptions.type,'simulated') || ...
   strcmp(dataOptions.type,'synthetic')
    dataOptions.snrValues = 0;
end
if strcmp(dataOptions.type,'synthetic')
    dataOptions.samplingFrequency = 48000;
end
if strcmp(dataOptions.type,'simulated')
    dataOptions.wavFolder = '/scratch/pictor/deleforg/Sounds/TIMIT_normalized/';
    dataOptions.subIndices = [];
end
if strcmp(dataOptions.type,'real')
    dataOptions.wavRootFolder = strcat(theRootFolder,'audispace/gTDEJournalData/');
    dataOptions.subIndices = [];
end
dataOptions.t60 = 0.6;
% Save it
experimentOptions.dataOptions = dataOptions;

%%% Microphone properties
microphonePositionOptions.type = 'tetrahedron';
microphonePositionOptions.scale = 0.1;
microphonePositionOptions.offset = [1.9 2.1 1.9];
% Subset of the microphone set
% 4 Microphones
% microphonePositionOptions.subSet = [1,2,11,12];
% 6 Microphones
% microphonePositionOptions.subSet = [1,3,7,8,11,12];
% 8 Microphones
% microphonePositionOptions.subSet = [5,6,7,8,9,10,11,12];
% Increasing # of microphones
% microphonePositionOptions.subSet = [1,2,11,12,6,7,4,10,3,9,5,6];

% Save itf
experimentOptions.microphonePositionOptions = microphonePositionOptions;

%%% ISM Options
ismOptions.room = room;
ismOptions.absorptionWeights =  [1 1 1 1 1 1];
ismOptions.t60 = [0 0.1 0.2 0.4 0.6];
ismOptions.samplingFrequencies = [16000 44100 48000];
ismOptions.folder= strcat(experimentOptions.rootFolder,'SPHERICAL-2013/');
% Save it
experimentOptions.ismOptions = ismOptions;

%%% Method options
methodOptions.type = 'dip';
methodOptions.gridSize = 10;
methodOptions.regConstant = 10;
% Save it
experimentOptions.methodOptions = methodOptions;

%%% Result options

%%% Run the batch
batchTDEExperiments(experimentOptions);
