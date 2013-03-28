function showResultsTDEExperiments(experimentOptions,foundTDEs)

    if ~isfield(experimentOptions,'rootDir')
        experimentOptions.rootDir = '/local_scratch/alamedap/Software/gTDESSL/';
    end

    %%% General variables
    % Check options
%     checkExperimentOptions(experimentOptions);

    % Get Source positions
%     experimentOptions.sourcePositions = getSourcePositions(experimentOptions);
%     nSourcePositions = size(experimentOptions.sourcePositions,1);
    
    % Get microphones positions
%     experimentOptions.microphonesPositions = getMicrophonesPositions(experimentOptions);
    
    % Generate Initialization positions
%     experimentOptions.initializationPositions = getInitializationPositions(experimentOptions);
    if strcmp(experimentOptions.dataOptions.type,'real')
        experimentOptions.T60 = 0.2;
        experimentOptions.snrValues = 0;
    end
    
    % Get signals
    nSignals = numel(experimentOptions.signals);

    % SNR values
    nSNRValues = length(experimentOptions.snrValues);
    
    % T60 values
    nT60 = numel(experimentOptions.T60);
    
    %%% First, bias and standard deviation
    % Compute true TDE
    trueTDEs = TDEGeometricDirect(experimentOptions.sourcePositions,...
                                 experimentOptions.microphonePositions);    
%     % PLOT BY SIGNALS
%     plotBySignal(foundTDEs,experimentOptions,trueTDEs);

%     % PLOT LOCALIZATION
%     plotLocalizationError(foundTDEs,experimentOptions);
    
%     % Plot by T60
%     plotByT60(foundTDEs,experimentOptions,trueTDEs);

%     % Plot by method
%     plotByMethod(foundTDEs,experimentOptions,trueTDEs);

%     % Plot by TDE
%     plotByTDE(foundTDEs,experimentOptions,trueTDEs);
    
end % function

% function checkExperimentOptions(experimentOptions)
%     
%     % Check that there is the kind of data used and the type of sensor
%     mandatoryFields = {'dataUsed','sensorType','dimension','sourcePositionOptions','snrValues'};
%     for f = 1:numel(mandatoryFields),
%         if ~isfield(experimentOptions,mandatoryFields{f})
%             error(['Field ' mandatoryFields{f} ' is mandatory.']);
%         end
%     end
% 
%     data = {'synthetic','simulated','real'};
%     % Check the data specified
%     dataFlag = 0;
%     for d = 1:numel(data)
%         dataFlag = dataFlag + strcmp(experimentOptions.dataUsed,data{d});
%     end
%     if dataFlag == 0
%         error('Data specified not known.');
%     end
% %     if strcmp(experimentOptions.dataUsed,data{3})
% %         error('Not ready to work with real data.');
% %     end
% end

% The SNR values will be in the x-axis, the T60 values will produce one
% curve each. The plots will be error bars
function plotBySignal(foundTDEs,experimentOptions,trueTDEs)
    % Variables
    nSNR = length(experimentOptions.snrValues);
    nT60 = length(experimentOptions.T60);
    if strcmp(experimentOptions.dataUsed,'simulated')
        nSignals = 3;
    elseif strcmp(experimentOptions.dataUsed,'real')
        nSignals = 1;
    else
        nSignals = numel(experimentOptions.signals);
    end
    % Plot variables
%     lineSpec = {'-o',':x','-.s','--d'};
    lineSpec = {'-o','-d','-s','-p'};
    tdeColor = {'g','b','r','k'};
%     % Random comparison
%     randomTDE = generateRandomEstimates(experimentOptions.microphonesPositions,experimentOptions.constraint,10*size(experimentOptions.sourcePositions,1));
%     randomErrorTDEs = [];
%     for sPos = 1:size(experimentOptions.sourcePositions,1),
%         randomErrorTDEs = cat(1,randomErrorTDEs,...
%             abs(randomTDE( (10*(sPos-1)+1):10*sPos,:)-repmat(trueTDEs(sPos,1:experimentOptions.dimension),10,1)));
%     end
%     rErrorTDE = cat(1,mean(randomErrorTDEs,1),std(randomErrorTDEs,0,1));
    % For each signal
    for sSignal = nSignals:nSignals,
        % Signal results
        signalErrorTDEs = cell(nSNR,nT60);
        errorTDEs = zeros(nSNR,nT60,experimentOptions.dimension+1,2);
        % Recover the found TDEs for this signal
        for sSNR = 1:nSNR,
            for sT60 = 1:nT60,
                % Compute the difference with each instance
                for sPos = 1:size(experimentOptions.sourcePositions,1),
                    signalErrorTDEs{sSNR,sT60} = cat(1,signalErrorTDEs{sSNR,sT60},...
                        abs(foundTDEs{sSignal,sPos,sSNR,sT60}-repmat(trueTDEs(sPos,1:experimentOptions.dimension),size(foundTDEs{sSignal,sPos,sSNR,sT60},1),1)));
                end
                % Compute error statistics
                errorTDEs(sSNR,sT60,1:end-1,1) = mean(signalErrorTDEs{sSNR,sT60},1);
                errorTDEs(sSNR,sT60,1:end-1,2) = std(signalErrorTDEs{sSNR,sT60},0,1);
                errorTDEs(sSNR,sT60,end,1) = mean(signalErrorTDEs{sSNR,sT60}(:),1);
                errorTDEs(sSNR,sT60,end,2) = std(signalErrorTDEs{sSNR,sT60}(:),0,1);
            end
        end
        fprintf('Signal %d:\n',sSignal);
        fprintf('\t SNR \t'); for sSNR=1:nSNR, fprintf('%f\t',experimentOptions.snrValues(sSNR)); end; fprintf('\n');
        fprintf('\t TDE1 \t'); for sSNR=1:nSNR, fprintf('%1.2e (%1.2e)\t',errorTDEs(sSNR,end,1,1),errorTDEs(sSNR,end,1,2)); end; fprintf('\n');
        fprintf('\t TDE2 \t'); for sSNR=1:nSNR, fprintf('%1.2e (%1.2e)\t',errorTDEs(sSNR,end,2,1),errorTDEs(sSNR,end,2,2)); end; fprintf('\n');
        fprintf('\t TDE3 \t'); for sSNR=1:nSNR, fprintf('%1.2e (%1.2e)\t',errorTDEs(sSNR,end,3,1),errorTDEs(sSNR,end,3,2)); end; fprintf('\n');
        % For each dimension perform a plot
        figure
        hold on
        for d = 1:experimentOptions.dimension,
%             errorbar(experimentOptions.snrValues,...
%                      repmat(rErrorTDE(1,d),nSNR,1),...
%                      repmat(rErrorTDE(2,d),nSNR,1),...
%                      strcat(lineSpec{end},tdeColor{d}));
            for sT60 = 1:nT60,
                plot(experimentOptions.snrValues,...
                    errorTDEs(:,sT60,end,1),...
                    strcat(lineSpec{sT60},tdeColor{sT60}),...
                    'LineWidth',3,'MarkerFaceColor',tdeColor{sT60},'MarkerSize',15);
%                 errorbar(experimentOptions.snrValues,...
%                     errorTDEs(:,sT60,d,1),...
%                     errorTDEs(:,sT60,d,2),...
%                     strcat(lineSpec{sT60},tdeColor{d}));
            end
            axis([-16 11 0 4e-4]);
            set(gca,'FontSize',30);
            xlabel('SNR (dB)','fontsize',30);
            ylabel('Absolute error mean (s)','fontsize',30);
            legend('T60 = 0 s','T60 = 0.1 s','T60 = 0.2 s');%,'T60 = 0.4 s');
        end
        grid on
    end
    
end

function plotLocalizationError(foundTDEs,experimentOptions)
    % Variables
    nSNR = length(experimentOptions.snrValues);
    nT60 = length(experimentOptions.T60);
    if strcmp(experimentOptions.dataUsed,'simulated')
        nSignals = 3;
    elseif strcmp(experimentOptions.dataUsed,'real')
        nSignals = 1;
    else
        nSignals = numel(experimentOptions.signals);
    end
    % Microphones' origin
    origin = mean(experimentOptions.microphonesPositions,1);
    % Compute spherical coordinates
    sphericalCoordinates = computeSphericalCoordinates(experimentOptions.sourcePositions,origin);
    % Plot variables
    lineSpec = {'-o','-d','-s','-p'};
    tdeColor = {'g','b','r','k'};
    % For each signal
    for sSignal = nSignals:nSignals,
        errorLocalization = zeros(nSNR,nT60,3,2);
        % Recover the found TDEs for this signal
        for sSNR = 1:nSNR,
            for sT60 = 1:nT60,
                % Loc error
                locError = [];
                if sSNR == nSNR && sSignal == 1 %&& sT60 == 1 %&&  sSignal == 1
                    figure
                    hold on
                    % Plot spherical coordinates
%                     plot(sphericalCoordinates(:,1),...
%                          sphericalCoordinates(:,2),...
%                          'ko','MarkerFaceColor','k');
%                     plot3(experimentOptions.microphonesPositions(:,1),...
%                           experimentOptions.microphonesPositions(:,2),...
%                           experimentOptions.microphonesPositions(:,3),...
%                           'bs','MarkerFaceColor','b');
%                     view(-20,60);
%                     if strcmp(experimentOptions.dataUsed,'simulated')
%                         axis([0 experimentOptions.room(1) 0 experimentOptions.room(2) 0 experimentOptions.room(3)]);
%                     end
                end
                % Compute the difference with each instance
                for sPos = 1:size(experimentOptions.sourcePositions,1),
%                     if strcmp(experimentOptions.dataUsed,'real')
%                         figure
%                         plot(180*sphericalCoordinates(sPos,1)/pi,180*sphericalCoordinates(sPos,2)/pi,'gd');
%                         axis([-180 180 -90 90]);
%                         hold on
%                     end
                    % Number of sub-signals
                    nPartialTrials = size(foundTDEs{sSignal,sPos,sSNR,sT60},1);
                    % Localization partial error
                    locPartError = [];
                    for partialTrial = 1:nPartialTrials,
                        % Comptue position and its spherical coordinates
                        pos = TDEGeometricInverse(experimentOptions.microphonesPositions,foundTDEs{sSignal,sPos,sSNR,sT60}(partialTrial,:));
                        if(numel(pos))
                            % Centered position
                            pos = pos - origin;
                            % Spherical coordinates
                            [pos(1) pos(2) pos(3)] = cart2sph(pos(1), pos(2), pos(3));
                            % Angle error
                            partError = [angdist(pos(1),sphericalCoordinates(sPos,1)) abs(pos(2)-sphericalCoordinates(sPos,2))];
                            plotError = sqrt((partError(1)/pi)^2 + (partError(2)/(pi/2))^2);
                            % Store
                            locPartError = cat(1,locPartError,180*partError/pi);
                            % If any, plot
%                             if sSNR == nSNR && sT60 == 1
%                                 plot(pos(:,1),pos(:,2),'ks','MarkerFaceColor',[1-exp(-plotError);exp(-plotError);0]);
%                             end
%                             if strcmp(experimentOptions.dataUsed,'real')
%                                 plot(180*pos(1)/pi,180*pos(2)/pi,'ko');
%                             end
                        end
                    end
%                     If the max SNR and the last signal
                    if sSNR == nSNR && sSignal == 1 %&& sT60 == 1 %
%                         az = 180*sphericalCoordinates(sPos,1)/pi;
%                         el = 180*sphericalCoordinates(sPos,2)/pi;
%                         plot(az,el,'ko');
%                         plot([az-std(locPartError(:,1)),az+std(locPartError(:,1))],...
%                              [el,el],'bd-');
%                         plot([az,az],...
%                              [el-std(locPartError(:,2)),el+std(locPartError(:,2))],'gx-');
                        color = mean(locPartError(:,1))/180;
                        msize = mean(locPartError(:,2))+5;
                        %%%% 3D balls
%                         plot3(experimentOptions.sourcePositions(sPos,1),experimentOptions.sourcePositions(sPos,2),experimentOptions.sourcePositions(sPos,3),...
%                             'ko','MarkerFaceColor',[ 1-color, 1-color, 1-color],'MarkerSize',msize);
                        %%%% Plot agains elevation & azimuth
%                         plot(180*sphericalCoordinates(sPos,1)/pi,180*sphericalCoordinates(sPos,2)/pi,...
%                             'ko','MarkerFaceColor',[ 1-color, 1-color, 1-color],'MarkerSize',msize);
%                         plot(180*sphericalCoordinates(sPos,1)/pi,180*sphericalCoordinates(sPos,2)/pi,...
%                              'ko','MarkerSize',mean(locPartError(locPartError(:,1)<30,1))+5);
                        plot(sphericalCoordinates(sPos,3),mean(locPartError(locPartError(:,1)<30,1)),...
                             'ko');
%                         axis([-180 180 -90 90 0 30]);
                    end
                    % Cumulate the results
                    locError = cat(1,locError,locPartError);
                end
%                 if sSNR == nSNR
%                     figure
%                     hist(locError,100);
%                 end
                % Compute error statistics
                nonAnomalies = locError(:,1) < 30 & locError(:,2) < 30;
                errorLocalization(sSNR,sT60,1,:) = 100*sum(~nonAnomalies)/numel(nonAnomalies);
                errorLocalization(sSNR,sT60,2,:) = mean(locError(nonAnomalies,:),1);
                errorLocalization(sSNR,sT60,3,:) = std(locError(nonAnomalies,:),1);
            end
        end
        fprintf('Signal %d:\n',sSignal);
        fprintf('\t SNR \t'); for sSNR=1:nSNR, fprintf('%d\t\t',experimentOptions.snrValues(sSNR)); end; fprintf('\n');
        fprintf('\t SNR \t'); for sSNR=1:nSNR, fprintf('%1.2f%%\t',errorLocalization(sSNR,end,1,1)); end; fprintf('\n');
        fprintf('\t Az \t'); for sSNR=1:nSNR, fprintf('%1.2f (%1.2f)\t',errorLocalization(sSNR,end,2,1),errorLocalization(sSNR,end,3,1)); end; fprintf('\n');
        fprintf('\t El \t'); for sSNR=1:nSNR, fprintf('%1.2f (%1.2f)\t',errorLocalization(sSNR,end,2,2),errorLocalization(sSNR,end,3,2)); end; fprintf('\n');
        
        % SNR-plot
        % For each localization variables       
        %%% INLIERS
        figure
        hold on
        for sT60 = 1:nT60,
            plot(experimentOptions.snrValues,100-squeeze(errorLocalization(:,sT60,1,1)),...
                strcat(lineSpec{sT60},tdeColor{sT60}),'LineWidth',3,...
                'MarkerFaceColor',tdeColor{sT60},'MarkerSize',15);
            grid on
        end
        axis([-16 11 0 100]);
        set(gca,'FontSize',30);
        xlabel('SNR (dB)','fontsize',30);
        ylabel('Inliers (%)','fontsize',30);
        %%% AZIMUTH
        figure
        hold on
        for sT60 = 1:nT60,
            plot(experimentOptions.snrValues,squeeze(errorLocalization(:,sT60,2,1)),...
                strcat(lineSpec{sT60},tdeColor{sT60}),'LineWidth',3,...
                'MarkerFaceColor',tdeColor{sT60},'MarkerSize',15);
%             errorbar(experimentOptions.snrValues,...
%                 errorLocalization(:,sT60,2,1),...
%                 errorLocalization(:,sT60,3,1),...
%                 lineSpec{sT60});
            grid on
        end
        axis([-16 11 0 15]);
        set(gca,'FontSize',30);
        xlabel('SNR (dB)','fontsize',30);
        ylabel('Azimuth MAE','fontsize',30);
        %%% ELEVATION
        figure
        hold on
        for sT60 = 1:nT60,
            plot(experimentOptions.snrValues,squeeze(errorLocalization(:,sT60,3,1)),...
                strcat(lineSpec{sT60},tdeColor{sT60}),'LineWidth',3,...
                'MarkerFaceColor',tdeColor{sT60},'MarkerSize',15);
%             errorbar(experimentOptions.snrValues,...
%                 errorLocalization(:,sT60,2,2),...
%                 errorLocalization(:,sT60,3,2),...
%                 lineSpec{sT60});
            grid on
        end
        axis([-16 11 0 15]);
        set(gca,'FontSize',30);
        xlabel('SNR (dB)','fontsize',30);
        ylabel('Elevation MAE','fontsize',30);
    end
end

% The SNR values will be in the x-axis, the T60 values will produce one
% curve each. The plots will be error bars
function plotByT60(foundTDEs,experimentOptions,trueTDEs)
    % Variables
    nSNR = length(experimentOptions.snrValues);
    nT60 = length(experimentOptions.T60);
    
    if strcmp(experimentOptions.dataUsed,'simulated')
        nSignals = 3;
    elseif strcmp(experimentOptions.dataUsed,'real')
        nSignals = 1;
    else
        nSignals = numel(experimentOptions.signals);
    end
    % Microphones' origin
    origin = mean(experimentOptions.microphonesPositions,1);
    % MAXTDE
    maxTDE = TDEmax(experimentOptions.microphonesPositions);
    % Compute spherical coordinates
    sphericalCoordinates = computeSphericalCoordinates(experimentOptions.sourcePositions,origin);
    % Plot variables
%     lineSpec = {'-o',':x','-.s','--d'};
    lineSpec = {'-d','-s','-o','^-'};
    tdeColor = {'g','b','r','k'};
    % For each signal
    for sSignal = nSignals:nSignals,
        % Error in localization
        errorLocalization = zeros(nSNR,nT60,3,2);
        % Recover the found TDEs for this signal
        % Error in tdes
        signalErrorTDEs = cell(nSNR,nT60);
        satisfyConstraint = cell(nSNR,nT60);
        errorTDEs = zeros(nSNR,nT60,experimentOptions.dimension+1,4);
        % Recover the found TDEs for this signal
        for sSNR = 1:nSNR,
            for sT60 = 1:nT60,
                %%% TDE ERROR
                % Compute the difference with each instance
                for sPos = 1:size(experimentOptions.sourcePositions,1),
                    % Store errors in this positions
                    partialSignalErrorTDE = foundTDEs{sSignal,sPos,sSNR,sT60}-repmat(trueTDEs(sPos,1:experimentOptions.dimension),size(foundTDEs{sSignal,sPos,sSNR,sT60},1),1);
                    % Concatenate to the rest
                    signalErrorTDEs{sSNR,sT60} = cat(1,signalErrorTDEs{sSNR,sT60},partialSignalErrorTDE);
                    constraint = foundTDEs{sSignal,sPos,sSNR,sT60} > -repmat(maxTDE(1:3)',size(foundTDEs{sSignal,sPos,sSNR,sT60},1),1) &...
                        foundTDEs{sSignal,sPos,sSNR,sT60} < repmat(maxTDE(1:3)',size(foundTDEs{sSignal,sPos,sSNR,sT60},1),1);
                    for d=2:experimentOptions.dimension,
                        constraint(:,1) = constraint(:,1) & constraint(:,d);
                    end
                    constraint = constraint(:,1);% & (TDEDiscriminant(foundTDEs{sSignal,sPos,sSNR,sT60}',experimentOptions.microphonesPositions) > 0)';
                    for d=1:experimentOptions.dimension,
                        constraint = constraint & abs(partialSignalErrorTDE(:,d)) < 1e-4;
                    end
                    satisfyConstraint{sSNR,sT60} = cat(1,satisfyConstraint{sSNR,sT60},constraint);
                end
                signalErrorTDEs{sSNR,sT60} = signalErrorTDEs{sSNR,sT60}(logical(satisfyConstraint{sSNR,sT60}),:);
                % Compute error statistics
                errorTDEs(sSNR,sT60,1:end-1,1) = abs(mean(signalErrorTDEs{sSNR,sT60},1));
                errorTDEs(sSNR,sT60,1:end-1,2) = std(signalErrorTDEs{sSNR,sT60},0,1);
                errorTDEs(sSNR,sT60,1:end-1,3) = sqrt(mean(signalErrorTDEs{sSNR,sT60}.^2,1));
                errorTDEs(sSNR,sT60,end,1) = abs(mean(signalErrorTDEs{sSNR,sT60}(:)));
                errorTDEs(sSNR,sT60,end,2) = std(signalErrorTDEs{sSNR,sT60}(:),0);
                errorTDEs(sSNR,sT60,end,3) = sqrt(mean((signalErrorTDEs{sSNR,sT60}(:)).^2));

                errorTDEs(sSNR,sT60,1:end,4) = 100*mean(satisfyConstraint{sSNR,sT60}(:));
                %%%% LOCALIZATION
                % Loc error
                locError = [];
                nonAnomalies = [];
                for sPos = 1:size(experimentOptions.sourcePositions,1),
                    % Number of sub-signals
                    nPartialTrials = size(foundTDEs{sSignal,sPos,sSNR,sT60},1);
                    % Localization partial error
                    locPartError = [];
                    for partialTrial = 1:nPartialTrials,
                        % Comptue position and its spherical coordinates
                        pos = TDEGeometricInverse(experimentOptions.microphonesPositions,foundTDEs{sSignal,sPos,sSNR,sT60}(partialTrial,:));
                        if(numel(pos))
                            % Centered position
                            pos = pos - origin;
                            % Spherical coordinates
                            [pos(1) pos(2) pos(3)] = cart2sph(pos(1), pos(2), pos(3));
                            % Angle error
                            partError = [angdist(pos(1),sphericalCoordinates(sPos,1)) abs(pos(2)-sphericalCoordinates(sPos,2))];
                            plotError = sqrt((partError(1)/pi)^2 + (partError(2)/(pi/2))^2);
                            % Store
                            locPartError = cat(1,locPartError,180*partError/pi);
                            nonAnomalies = cat(1,nonAnomalies,true);
                        else
                            nonAnomalies = cat(1,nonAnomalies,false);
                            locPartError = cat(1,locPartError,[0,0]);
                        end
                    end
                    % Cumulate the results
                    locError = cat(1,locError,locPartError);
                end
                % Compute error statistics
                nonAnomalies = nonAnomalies & locError(:,1) < 30 & locError(:,2) < 30 & satisfyConstraint{sSNR,sT60};
                errorLocalization(sSNR,sT60,1,:) = 100*sum(~nonAnomalies)/numel(nonAnomalies);
                errorLocalization(sSNR,sT60,2,:) = mean(locError(nonAnomalies,:),1);
                errorLocalization(sSNR,sT60,3,:) = std(locError(nonAnomalies,:),1);
            end
        end
        
%         histTimes = (-1e-3):1e-6:1e-3;
%         for sMethod=1:nMethod,
%             figure
%             hist(signalErrorTDEs{nSNR,1}(:),histTimes);
%             title(experimentOptions.method{sMethod});
%             fprintf('Found in %s: %d.\n',experimentOptions.method{sMethod},size(signalErrorTDEs{nSNR,1},1));
%         end
        
        %%%%% BIAS
        % For each dimension perform a plot
%         for d = 1:experimentOptions.dimension,
            figure
            hold on
%             errorbar(experimentOptions.snrValues,...
%                      repmat(rErrorTDE(1,d),nSNR,1),...
%                      repmat(rErrorTDE(2,d),nSNR,1),...
%                      strcat(lineSpec{end},tdeColor{d}));
            for sT60 = 1:nT60,
                plot(experimentOptions.snrValues,...
                    errorTDEs(:,sT60,end,1),...
                    strcat(lineSpec{sT60},tdeColor{sT60}),...
                    'LineWidth',3,'MarkerFaceColor',tdeColor{sT60},'MarkerSize',15);
%                 errorbar(experimentOptions.snrValues+sMethod*0.05,...
%                     errorTDEs(:,sT60,d,1),...
%                     errorTDEs(:,sT60,d,2),...
%                     strcat(lineSpec{sMethod},tdeColor{end}),...
%                     'LineWidth',3,'MarkerFaceColor',tdeColor{end},'MarkerSize',15);
            end
%             axis([-16 11 0 4e-4]);
            set(gca,'FontSize',30);
            xlabel('SNR (dB)','fontsize',30);
            ylabel('TDE Bias (s)','fontsize',30);
            legend('T60=0','T60=0.1','T60=0.2');
%         end
%         grid on
        
        %%%%% STD
        % For each dimension perform a plot
%         for d = 1:experimentOptions.dimension,
            figure
            hold on
%             errorbar(experimentOptions.snrValues,...
%                      repmat(rErrorTDE(1,d),nSNR,1),...
%                      repmat(rErrorTDE(2,d),nSNR,1),...
%                      strcat(lineSpec{end},tdeColor{d}));
            for sT60 = 1:nT60,
                plot(experimentOptions.snrValues,...
                    errorTDEs(:,sT60,end,2),...
                    strcat(lineSpec{sT60},tdeColor{sT60}),...
                    'LineWidth',3,'MarkerFaceColor',tdeColor{sT60},'MarkerSize',15);
%                 errorbar(experimentOptions.snrValues+sMethod*0.05,...
%                     errorTDEs(:,sT60,d,1),...
%                     errorTDEs(:,sT60,d,2),...
%                     strcat(lineSpec{sMethod},tdeColor{end}),...
%                     'LineWidth',3,'MarkerFaceColor',tdeColor{end},'MarkerSize',15);
            end
%             axis([-16 11 0 4e-4]);
            set(gca,'FontSize',30);
            xlabel('SNR (dB)','fontsize',30);
            ylabel('TDE Standard Deviation (s)','fontsize',30);
            legend('T60=0','T60=0.1','T60=0.2');
%         end
%         grid on
        
        %%%%% Constraint
        % For each dimension perform a plot
%         for d = 1:experimentOptions.dimension,
            figure
            hold on
%             errorbar(experimentOptions.snrValues,...
%                      repmat(rErrorTDE(1,d),nSNR,1),...
%                      repmat(rErrorTDE(2,d),nSNR,1),...
%                      strcat(lineSpec{end},tdeColor{d}));
            for sT60 = 1:nT60,
                plot(experimentOptions.snrValues,...
                    errorTDEs(:,sT60,end,4),...
                    strcat(lineSpec{sT60},tdeColor{sT60}),...
                    'LineWidth',3,'MarkerFaceColor',tdeColor{sT60},'MarkerSize',15);
%                 errorbar(experimentOptions.snrValues+sMethod*0.05,...
%                     errorTDEs(:,sT60,d,1),...
%                     errorTDEs(:,sT60,d,2),...
%                     strcat(lineSpec{sMethod},tdeColor{end}),...
%                     'LineWidth',3,'MarkerFaceColor',tdeColor{end},'MarkerSize',15);
            end
            axis([-10 10 0 100]);
            set(gca,'FontSize',30);
            xlabel('SNR (dB)','fontsize',30);
            ylabel('Non-anomalous points (%)','fontsize',30);
            legend('T60=0','T60=0.1','T60=0.2');
%         end
%         grid on
        
        fprintf('Signal %d:\n',sSignal);
        for sT60 = 1:nT60
            fprintf('  T60 %f:\n',experimentOptions.T60(sT60));
            fprintf('\t SNR \t'); for sSNR=1:nSNR, fprintf('%d\t\t',experimentOptions.snrValues(sSNR)); end; fprintf('\n');
            fprintf('\t INL \t'); for sSNR=1:nSNR, fprintf('%1.2f%%\t\t',100-errorLocalization(sSNR,sT60,1,1)); end; fprintf('\n');
            fprintf('\t Az \t'); for sSNR=1:nSNR, fprintf('%1.2f (%1.2f)\t',errorLocalization(sSNR,sT60,2,1),errorLocalization(sSNR,sT60,3,1)); end; fprintf('\n');
            fprintf('\t El \t'); for sSNR=1:nSNR, fprintf('%1.2f (%1.2f)\t',errorLocalization(sSNR,sT60,2,2),errorLocalization(sSNR,sT60,3,2)); end; fprintf('\n');
        end
        
        % SNR-plot
        % For each localization variables       
        %%% INLIERS
        figure
        hold on
        for sT60 = 1:nT60,
            plot(experimentOptions.snrValues,100-squeeze(errorLocalization(:,sT60,1,1)),...
                strcat(lineSpec{sT60},tdeColor{sT60}),'LineWidth',3,...
                'MarkerFaceColor',tdeColor{sT60},'MarkerSize',15);
%             grid on
        end
        legend('T60=0','T60=0.1','T60=0.2');
        axis([-10 10 0 100]);
        set(gca,'FontSize',30);
        xlabel('SNR (dB)','fontsize',30);
        ylabel('Inliers (%)','fontsize',30);
        %%% AZIMUTH
        figure
        hold on
        for sT60 = 1:nT60,
            plot(experimentOptions.snrValues,squeeze(errorLocalization(:,sT60,2,1)),...
                strcat(lineSpec{sT60},tdeColor{sT60}),'LineWidth',3,...
                'MarkerFaceColor',tdeColor{sT60},'MarkerSize',15);
%             errorbar(experimentOptions.snrValues,...
%                 errorLocalization(:,sT60,2,1),...
%                 errorLocalization(:,sT60,3,1),...
%                 lineSpec{sT60});
%             grid on
        end
        legend('T60=0','T60=0.1','T60=0.2');
        axis([-10 10 0 15]);
        set(gca,'FontSize',30);
        xlabel('SNR (dB)','fontsize',30);
        ylabel('Azimuth MAE','fontsize',30);
        %%% ELEVATION
        figure
        hold on
        for sT60 = 1:nT60,
            plot(experimentOptions.snrValues,squeeze(errorLocalization(:,sT60,3,1)),...
                strcat(lineSpec{sT60},tdeColor{sT60}),'LineWidth',3,...
                'MarkerFaceColor',tdeColor{sT60},'MarkerSize',15);
%             errorbar(experimentOptions.snrValues,...
%                 errorLocalization(:,sT60,2,2),...
%                 errorLocalization(:,sT60,3,2),...
%                 lineSpec{sT60});
%             grid on
        end
        legend('T60=0','T60=0.1','T60=0.2');
        axis([-10 10 0 10]);
        set(gca,'FontSize',30);
        xlabel('SNR (dB)','fontsize',30);
        ylabel('Elevation MAE','fontsize',30);
    end
    
end

% The SNR values will be in the x-axis, the T60 values will produce one
% curve each. The plots will be error bars
function plotByMethod(foundTDEs,experimentOptions,trueTDEs)
    % Variables
    nSNR = length(experimentOptions.snrValues);
    nT60 = length(experimentOptions.T60);
    nMethod = length(experimentOptions.method);
    
    if strcmp(experimentOptions.dataUsed,'simulated')
        nSignals = 3;
    elseif strcmp(experimentOptions.dataUsed,'real')
        nSignals = 1;
    else
        nSignals = numel(experimentOptions.signals);
    end
    % Microphones' origin
    origin = mean(experimentOptions.microphonesPositions,1);
    % MAXTDE
    maxTDE = TDEmax(experimentOptions.microphonesPositions);
    % Compute spherical coordinates
    sphericalCoordinates = computeSphericalCoordinates(experimentOptions.sourcePositions,origin);
    % Plot variables
%     lineSpec = {'-o',':x','-.s','--d'};
    lineSpec = {'-o','--s','-d','^-'};
    tdeColor = {'g','b','r','k'};
    % For each signal
    for sSignal = nSignals:nSignals,
        % Error in localization
        errorLocalization = zeros(nSNR,nT60,nMethod,3,2);
        % Recover the found TDEs for this signal
        % Error in tdes
        signalErrorTDEs = cell(nSNR,nT60,nMethod);
        satisfyConstraint = cell(nSNR,nT60,nMethod);
        errorTDEs = zeros(nSNR,nT60,nMethod,experimentOptions.dimension+1,4);
        % Recover the found TDEs for this signal
        for sSNR = 1:nSNR,
            for sT60 = 1:nT60,
                for sMethod = 1:nMethod,
                    %%% TDE ERROR
                    % Compute the difference with each instance
                    for sPos = 1:size(experimentOptions.sourcePositions,1),
                        % Store errors in this positions
                        partialSignalErrorTDE = foundTDEs{sSignal,sPos,sSNR,sT60,sMethod}-repmat(trueTDEs(sPos,1:experimentOptions.dimension),size(foundTDEs{sSignal,sPos,sSNR,sT60,sMethod},1),1);
                        % Concatenate to the rest
                        signalErrorTDEs{sSNR,sT60,sMethod} = cat(1,signalErrorTDEs{sSNR,sT60,sMethod},partialSignalErrorTDE);
                        constraint = foundTDEs{sSignal,sPos,sSNR,sT60,sMethod} > -repmat(maxTDE(1:3)',size(foundTDEs{sSignal,sPos,sSNR,sT60,sMethod},1),1) &...
                            foundTDEs{sSignal,sPos,sSNR,sT60,sMethod} < repmat(maxTDE(1:3)',size(foundTDEs{sSignal,sPos,sSNR,sT60,sMethod},1),1);
                        for d=2:experimentOptions.dimension,
                            constraint(:,1) = constraint(:,1) & constraint(:,d);
                        end
                        constraint = constraint(:,1);% & (TDEDiscriminant(foundTDEs{sSignal,sPos,sSNR,sT60,sMethod}',experimentOptions.microphonesPositions) > 0)';
                        for d=1:experimentOptions.dimension,
                            constraint = constraint & abs(partialSignalErrorTDE(:,d)) < 1e-4;
                        end
                        satisfyConstraint{sSNR,sT60,sMethod} = cat(1,satisfyConstraint{sSNR,sT60,sMethod},constraint);
                    end
                    signalErrorTDEs{sSNR,sT60,sMethod} = signalErrorTDEs{sSNR,sT60,sMethod}(logical(satisfyConstraint{sSNR,sT60,sMethod}),:);
                    % Compute error statistics
                    errorTDEs(sSNR,sT60,sMethod,1:end-1,1) = abs(mean(signalErrorTDEs{sSNR,sT60,sMethod},1));
                    errorTDEs(sSNR,sT60,sMethod,1:end-1,2) = std(signalErrorTDEs{sSNR,sT60,sMethod},0,1);
                    errorTDEs(sSNR,sT60,sMethod,1:end-1,3) = sqrt(mean(signalErrorTDEs{sSNR,sT60,sMethod}.^2,1));
                    errorTDEs(sSNR,sT60,sMethod,end,1) = abs(mean(signalErrorTDEs{sSNR,sT60,sMethod}(:)));
                    errorTDEs(sSNR,sT60,sMethod,end,2) = std(signalErrorTDEs{sSNR,sT60,sMethod}(:),0);
                    errorTDEs(sSNR,sT60,sMethod,end,3) = sqrt(mean((signalErrorTDEs{sSNR,sT60,sMethod}(:)).^2));
                    
                    errorTDEs(sSNR,sT60,sMethod,1:end,4) = 100*mean(satisfyConstraint{sSNR,sT60,sMethod}(:));
                    %%%% LOCALIZATION
                    % Loc error
                    locError = [];
                    nonAnomalies = [];
                    for sPos = 1:size(experimentOptions.sourcePositions,1),
                        % Number of sub-signals
                        nPartialTrials = size(foundTDEs{sSignal,sPos,sSNR,sT60,sMethod},1);
                        % Localization partial error
                        locPartError = [];
                        for partialTrial = 1:nPartialTrials,
                            % Comptue position and its spherical coordinates
                            pos = TDEGeometricInverse(experimentOptions.microphonesPositions,foundTDEs{sSignal,sPos,sSNR,sT60,sMethod}(partialTrial,:));
                            if(numel(pos))
                                % Centered position
                                pos = pos - origin;
                                % Spherical coordinates
%                                 [pos(1) pos(2) pos(3)] = cart2sph(pos(1), pos(2), pos(3));
                                % Angle error
%                                 partError = [angdist(pos(1),sphericalCoordinates(sPos,1)) abs(pos(2)-sphericalCoordinates(sPos,2))];
%                                 partError = [pos(1)-sphericalCoordinates(sPos,1) pos(2)-sphericalCoordinates(sPos,2)];
%                                 plotError = sqrt((partError(1)/pi)^2 + (partError(2)/(pi/2))^2);
                                partError = acosd( sum( pos.*(experimentOptions.sourcePositions(sPos,:)-origin)) / (norm(pos)*norm(experimentOptions.sourcePositions(sPos,:)-origin)) );
                                % Store
                                locPartError = cat(1,locPartError,partError);%,180*partError/pi);
                                nonAnomalies = cat(1,nonAnomalies,true);
                            else
                                nonAnomalies = cat(1,nonAnomalies,false);
                                locPartError = cat(1,locPartError,0);
                            end
                        end
                        % Cumulate the results
                        locError = cat(1,locError,locPartError);
                    end
%                     % Correc the error
%                     locError(locError>pi) = locError(locError>pi)-2*pi;
%                     locError(locError<-pi) = locError(locError<-pi)+2*pi;
                    % Compute error statistics
                    nonAnomalies = nonAnomalies & locError(:,1) < 30 & satisfyConstraint{sSNR,sT60,sMethod};
                    errorLocalization(sSNR,sT60,sMethod,1,:) = 100*sum(~nonAnomalies)/numel(nonAnomalies);
                    errorLocalization(sSNR,sT60,sMethod,2,:) = mean(locError(nonAnomalies,:),1);
                    errorLocalization(sSNR,sT60,sMethod,3,:) = std(locError(nonAnomalies,:),1);
                end
            end
        end
        
%         histTimes = (-1e-3):1e-6:1e-3;
%         for sMethod=1:nMethod,
%             figure
%             hist(signalErrorTDEs{nSNR,1,sMethod}(:),histTimes);
%             title(experimentOptions.method{sMethod});
%             fprintf('Found in %s: %d.\n',experimentOptions.method{sMethod},size(signalErrorTDEs{nSNR,1,sMethod},1));
%         end

        myLegend = cell(0);
        for sT60 = 1:nT60,
            for sMethod = 1:sMethod,
                myLegend = cat(1,myLegend,strcat(experimentOptions.method{sMethod},' ','T60=',num2str(experimentOptions.T60(sT60))));
            end
        end
        
        %%%%% BIAS
        % For each dimension perform a plot
%         for d = 1:experimentOptions.dimension,
            figure
            hold on
%             errorbar(experimentOptions.snrValues,...
%                      repmat(rErrorTDE(1,d),nSNR,1),...
%                      repmat(rErrorTDE(2,d),nSNR,1),...
%                      strcat(lineSpec{end},tdeColor{d}));
            for sT60 = 1:nT60
                for sMethod = 1:nMethod,
                    plot(experimentOptions.snrValues,...
                        errorTDEs(:,sT60,sMethod,end,1),...
                        strcat(lineSpec{sT60},tdeColor{sMethod}),...
                        'LineWidth',3,'MarkerFaceColor',tdeColor{sMethod},'MarkerSize',15);
%                 errorbar(experimentOptions.snrValues+sMethod*0.05,...
%                     errorTDEs(:,sT60,sMethod,d,1),...
%                     errorTDEs(:,sT60,sMethod,d,2),...
%                     strcat(lineSpec{sMethod},tdeColor{end}),...
%                     'LineWidth',3,'MarkerFaceColor',tdeColor{end},'MarkerSize',15);
                end
            end
%             axis([-16 11 0 4e-4]);
            set(gca,'FontSize',30);
            xlabel('SNR (dB)','fontsize',30);
            ylabel('TDE Bias (s)','fontsize',30);
            legend(myLegend);
%         end
%         grid on
        
        %%%%% STD
        % For each dimension perform a plot
%         for d = 1:experimentOptions.dimension,
            figure
            hold on
%             errorbar(experimentOptions.snrValues,...
%                      repmat(rErrorTDE(1,d),nSNR,1),...
%                      repmat(rErrorTDE(2,d),nSNR,1),...
%                      strcat(lineSpec{end},tdeColor{d}));
            for sT60 = 1:nT60
                for sMethod = 1:nMethod,
                    plot(experimentOptions.snrValues,...
                        errorTDEs(:,sT60,sMethod,end,2),...
                        strcat(lineSpec{sT60},tdeColor{sMethod}),...
                        'LineWidth',3,'MarkerFaceColor',tdeColor{sMethod},'MarkerSize',15);
%                 errorbar(experimentOptions.snrValues+sMethod*0.05,...
%                     errorTDEs(:,sT60,sMethod,d,1),...
%                     errorTDEs(:,sT60,sMethod,d,2),...
%                     strcat(lineSpec{sMethod},tdeColor{end}),...
%                     'LineWidth',3,'MarkerFaceColor',tdeColor{end},'MarkerSize',15);
                end
            end
%             axis([-16 11 0 4e-4]);
            set(gca,'FontSize',30);
            xlabel('SNR (dB)','fontsize',30);
            ylabel('TDE Standard Deviation (s)','fontsize',30);
            legend(myLegend);
%         end
%         grid on
        
        %%%%% Constraint
        % For each dimension perform a plot
%         for d = 1:experimentOptions.dimension,
            figure
            hold on
%             errorbar(experimentOptions.snrValues,...
%                      repmat(rErrorTDE(1,d),nSNR,1),...
%                      repmat(rErrorTDE(2,d),nSNR,1),...
%                      strcat(lineSpec{end},tdeColor{d}));
            for sT60 = 1:nT60
                for sMethod = 1:nMethod,
                    plot(experimentOptions.snrValues,...
                        errorTDEs(:,sT60,sMethod,end,4),...
                        strcat(lineSpec{sT60},tdeColor{sMethod}),...
                        'LineWidth',3,'MarkerFaceColor',tdeColor{sMethod},'MarkerSize',15);
%                 errorbar(experimentOptions.snrValues+sMethod*0.05,...
%                     errorTDEs(:,sT60,sMethod,d,1),...
%                     errorTDEs(:,sT60,sMethod,d,2),...
%                     strcat(lineSpec{sMethod},tdeColor{end}),...
%                     'LineWidth',3,'MarkerFaceColor',tdeColor{end},'MarkerSize',15);
                end
            end
            axis([-10 10 0 100]);
            set(gca,'FontSize',30);
            xlabel('SNR (dB)','fontsize',30);
            ylabel('Non-anomalous points (%)','fontsize',30);
            legend(myLegend);
%         end
%         grid on
        
        fprintf('Signal %d:\n',sSignal);
        for sMethod = 1:nMethod
            fprintf('  Method %s:\n',experimentOptions.method{sMethod});
            fprintf('\t SNR \t'); for sSNR=1:nSNR, fprintf('%d\t\t',experimentOptions.snrValues(sSNR)); end; fprintf('\n');
            fprintf('\t INL \t'); for sSNR=1:nSNR, fprintf('%1.2f%%\t\t',100-errorLocalization(sSNR,end,sMethod,1,1)); end; fprintf('\n');
            fprintf('\t Az \t'); for sSNR=1:nSNR, fprintf('%1.2f (%1.2f)\t',errorLocalization(sSNR,end,sMethod,2,1),errorLocalization(sSNR,end,sMethod,3,1)); end; fprintf('\n');
            fprintf('\t El \t'); for sSNR=1:nSNR, fprintf('%1.2f (%1.2f)\t',errorLocalization(sSNR,end,sMethod,2,2),errorLocalization(sSNR,end,sMethod,3,2)); end; fprintf('\n');
        end
        
        % SNR-plot
        % For each localization variables       
        %%% INLIERS
        figure
        hold on
        for sT60 = 1:nT60
            for sMethod = 1:nMethod,
                plot(experimentOptions.snrValues,100-squeeze(errorLocalization(:,sT60,sMethod,1,1)),...
                    strcat(lineSpec{sT60},tdeColor{sMethod}),'LineWidth',3,...
                    'MarkerFaceColor',tdeColor{sMethod},'MarkerSize',15);
%             grid on
            end
        end
        legend(myLegend);
%         axis([-16 11 0 100]);
        set(gca,'FontSize',30);
        xlabel('SNR (dB)','fontsize',30);
        ylabel('Inliers (%)','fontsize',30);
        %%% AZIMUTH
        figure
        hold on
        for sT60 = 1:nT60
            for sMethod = 1:nMethod,
                plot(experimentOptions.snrValues,squeeze(errorLocalization(:,sT60,sMethod,3,1)),...
                    strcat(lineSpec{sT60},tdeColor{sMethod}),'LineWidth',3,...
                    'MarkerFaceColor',tdeColor{sMethod},'MarkerSize',15);
%             errorbar(experimentOptions.snrValues,...
%                 errorLocalization(:,sT60,2,1),...
%                 errorLocalization(:,sT60,3,1),...
%                 lineSpec{sT60});
%             grid on
            end
        end
        legend(myLegend);
%         axis([-16 11 0 15]);
        set(gca,'FontSize',30);
        xlabel('SNR (dB)','fontsize',30);
        ylabel('Direction MSE','fontsize',30);
%         %%% ELEVATION
%         figure
%         hold on
%         for sT60 = 1:nT60
%             for sMethod = 1:nMethod,
%                 plot(experimentOptions.snrValues,squeeze(errorLocalization(:,sT60,sMethod,3,2)),...
%                     strcat(lineSpec{sT60},tdeColor{sMethod}),'LineWidth',3,...
%                     'MarkerFaceColor',tdeColor{sMethod},'MarkerSize',15);
% %             errorbar(experimentOptions.snrValues,...
% %                 errorLocalization(:,sT60,2,2),...
% %                 errorLocalization(:,sT60,3,2),...
% %                 lineSpec{sT60});
% %             grid on
%             end
%         end
%         legend(myLegend);
% %         axis([-16 11 0 15]);
%         set(gca,'FontSize',30);
%         xlabel('SNR (dB)','fontsize',30);
%         ylabel('Elevation MAE','fontsize',30);
    end
    
end

function randomTDE = generateRandomEstimates(MICS,constraint,MinimumNumberEstimates)
    % Dimension
    Dimension = size(MICS,1)-1;
    % TDE Limits
    maxTDE = TDEmax(MICS);
    maxTDE = maxTDE(1:Dimension)';
    % If using the geometric constraint
    if constraint
        randomTDE = repmat(maxTDE,MinimumNumberEstimates,1).*(rand(MinimumNumberEstimates,Dimension)-0.5);
    else
        randomTDE = [];
        while size(randomTDE,1) < MinimumNumberEstimates
            newRandomTDE = repmat(maxTDE,MinimumNumberEstimates,1).*(rand(MinimumNumberEstimates,Dimension)-0.5);
            constraint = TDEDiscriminant(newRandomTDE',MICS);
            randomTDE = cat(1,randomTDE,newRandomTDE(constraint>0,:));
        end
    end
end