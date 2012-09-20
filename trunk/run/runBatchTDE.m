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
room = [3  4  2.5];

%%% Global options
% Root folder
[~, hostname] = system('hostname');
if strcmp(hostname(1:end-1),'cassiopeia') || strcmp(hostname(1:end-1),'peirasmos')
    experimentOptions.rootFolder = '/local_scratch/alamedap/Software/gtde/trunk/';
else
    cd('/scratch/cassiopeia');
    experimentOptions.rootFolder = '/scratch/cassiopeia/alamedap/Software/gtde/trunk/';
end
% Dimension of the ambient space
experimentOptions.dimension = 3;
% Use of ism
experimentOptions.ism = true;
% Add the root folder to the path (the software)
addpath(svngenpath(experimentOptions.rootFolder));

%%% Source position options
% Source position options
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
dataOptions.type = 'simulated';
% Length of the signal cuts
dataOptions.cutLength = 0.1;
% Depending
if strcmp(dataOptions.type,'simulated') || ...
   strcmp(dataOptions.type,'synthetic')
    dataOptions.snrValues = 10;
end
if strcmp(dataOptions.type,'synthetic')
    dataOptions.samplingFrequency = 48000;
end
if strcmp(dataOptions.type,'simulated')
    dataOptions.wavFolder = '/scratch/pictor/deleforg/the_CASA_REDMINE/emitted_sounds/cocktail_sounds_normalized/';
    dataOptions.subIndices = [];
end
% Save it
experimentOptions.dataOptions = dataOptions;

%%% Microphone properties
microphonePositionOptions.type = 'tetrahedron';
microphonePositionOptions.scale = 0.1;
microphonePositionOptions.offset = [2.25 1.25 1.25];
% Save it
experimentOptions.microphonePositionOptions = microphonePositionOptions;

%%% ISM Options
ismOptions.room = room;
ismOptions.absorptionWeights =  [1 1 1 1 1 1];
ismOptions.t60 = [0 0.1 0.2 0.4 0.6];
ismOptions.samplingFrequencies = [16000 44100 48000];
ismOptions.folder= strcat(experimentOptions.rootFolder,'EUSIPCO-2012/');
% Save it
experimentOptions.ismOptions = ismOptions;

%%% Method options
methodOptions.type = 'gtde';
% Save it
experimentOptions.methodOptions = methodOptions;

%%% Result options


%%% Run the batch
batchTDEExperiments(experimentOptions);
