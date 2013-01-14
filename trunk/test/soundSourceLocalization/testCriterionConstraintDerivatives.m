function testCriterionConstraintDerivatives

    clc; close all;

    % Scale parameters
    a = 0.1;
    c = a / sqrt(2);
    % Microphones in a regular tetrahedron
    MICS = [ a   0   -c ;...
            -a   0   -c ;...
             0    a   c ;...
             0   -a   c];
    % Random rotation
    Vx = rand(3,1)-0.5;
    Vy = rand(3,1)-0.5;
    alpha = 2*pi*rand(1);
    r = rotation(Vx,Vy,alpha);
    MICS = MICS * r;

    % Generate potision2
    positionOptions = struct('coordinateSystem','cartesian');
    % Declare the bounds in the three variables
    positionOptions.bounds = cell(1,3);
    % Set the bounds for the radius, the azimuth and elevation will de the
    % defaults
    positionOptions.bounds{1} = [-2 2];
    positionOptions.bounds{2} = [-2 2];
    positionOptions.bounds{3} = [-2 2];
    % Number of intervals
    positionOptions.numberOfIntervals = [8, 8, 8];
    % Do not use the evelation's extrema
    positionOptions.useExtrema = true(1,3);
    positionOptions.useExtrema(3) = false;
    % Generate positions
    global Positions;
    Positions = GeneratePositions(3, positionOptions);
    Positions(284,:) = [];
    

    % Choose the position
    truePosition = 56;

    % Sinusoidal
    % 120Hz
    % myFun = @(x) SyntheticSignal(x,1,120);
    % Sinusoidal with exponential decay
    % 120Hz and decay constant = 3
    % myFun = @(x) SyntheticSignal(x,2,[120,3]);
    % Exponential decay
    myFun = @(x) SyntheticSignal(x,2,[2000,200]);

    % Sampling frequency
    samplingPeriod = 1/48000;

    % Length
    Length = 0.01;

    % Generate the signals
    [signals] = GenerateDiscreteSignals(Positions(5,:),...
                                      MICS,...
                                      samplingPeriod,...
                                      Length,...
                                      myFun);
%     signals = 0.000001*randn(size(signals));
    %%% General variables
    % Number of microphones
    NMics = size(MICS,1);
    % Number of samples
    NSamples = size(signals,2);
    % Sampling times
    samplingTimes = 0:samplingPeriod:(NSamples-1)*samplingPeriod;
    % Maximum TDE values from the microphones
    maxTDE = TDEmax(MICS);
    maxTDESamples = ceil(1.1*maxTDE/samplingPeriod);
    maxLAG = max(maxTDESamples);
    
    %%% Compute the interpolation coefficients
    % Allocate
    PC = cell(NMics,1);
    % Compute
    for ss = 1:NMics,
        PC{ss} = PolynomialInterpolationCoefficients(signals(ss,:),samplingPeriod);
    end
    
    % Compute the Cross-correlation of the interpolation coefficients
    % Allocate
    PCCC = cell(NMics);
    % Compute
    for mic1 = 1:NMics,
        for mic2 = mic1:NMics,
            PCCC{mic1,mic2} = PolynomialCoefficientsCrossCorrelation(PC{mic1},PC{mic2},maxLAG);
        end
    end
                   
    X0 = TDEGeometricDirect(Positions,MICS);
    X0 = (X0(:,1:3)/samplingPeriod)';
    step = 1e-10;
        
    % Test gradient Objective function
    gradObjTest = @(x) gTDECriterion(x,PCCC,MICS,samplingPeriod);
    checkDerivatives(gradObjTest,X0,step,'Gradient');
    
    % Test hessian objective function
    for ii = 1:3,
        funcTest = @(x) hessObjTest(x,ii,PCCC,MICS,samplingPeriod);
        checkDerivatives(funcTest,X0,step,strcat('Hessian (',num2str(ii),')'));
    end

%     % Test gradient log-Objective function
%     gradObjTest = @(x) gTDELogCriterion(x,PCCC,MICS,samplingPeriod);
%     checkDerivatives(gradObjTest,X0,step,'Log-Gradient');
%     
%     % Test hessian log-objective function
%     for ii = 1:3,
%         funcTest = @(x) hessObjLogTest(x,ii,PCCC,MICS,samplingPeriod);
%         checkDerivatives(funcTest,X0,step,strcat('Log-Hessian (',num2str(ii),')'));
%     end

%     
%     % Test gradient's Constraint
%     X0 = TDEGeometricDirect(Positions(58:59,:),MICS);
%     X0 = (X0(:,1:3)/samplingPeriod)';
%     
%     gradConsTest = @(x) TDEDiscriminant(x,MICS,samplingPeriod);
%     checkManyDerivatives(gradConsTest,X0,step,'Gradient');
%     
% %     Test hessian's constraint
%     for ii = 1:3,
%         funcTest = @(x) hessConsTest(x,ii,MICS,samplingPeriod);
%         checkManyDerivatives(funcTest,X0,step,strcat('Hessian (',num2str(ii),')'));
%     end
    
%     % Test gradient chen's objective function
%     gradChenTest = @(x) chenTDECriterion(x,PCCC,[0.1,0.2,0.1],samplingPeriod);
%     checkDerivatives(gradChenTest,0,step,'Chen Gradient');
%     
%     % Test chen's hessian
%     funcTest = @(x) hessChenTest(x,PCCC,[0.1,0.2,0.1],samplingPeriod);
%     checkDerivatives(funcTest,0,step,'Chen Hessian');

%     % Test gradient Objective function
%     checkManyDerivatives(gradObjTest,X0,step,'Gradient');
%     
%     % Test hessian objective function
%     for ii = 1:3,
%         funcTest = @(x) hessObjTest(x,ii,PCCC,MICS,samplingPeriod);
%         checkManyDerivatives(funcTest,X0,step,strcat('Hessian (',num2str(ii),')'));
%     end
end

% function [C GC] = GradTDEDiscriminant(x,MICS,samplingPeriod)
%     [C GC] = TDEDiscriminant(x,MICS,samplingPeriod);
%     C = C';
%     GC = GC';
% end

function [F, GF] = hessObjTest(x,component,PCCC,MICS,samplingPeriod)
    [ddd, GO, HO] = gTDECriterion(x,PCCC,MICS,samplingPeriod);
    if size(x,2) == 1
        F = GO(component);
        GF = HO(component,:)';
    else
        F = GO(component,:);
        GF = squeeze(HO(component,:,:));
    end
end

function [F, GF] = hessObjLogTest(x,component,PCCC,MICS,samplingPeriod)
    [ddd, GO, HO] = gTDELogCriterion(x,PCCC,MICS,samplingPeriod);
    if size(x,2) == 1
        F = GO(component);
        GF = HO(component,:);
    else
        F = GO(component,:);
        GF = squeeze(HO(component,:,:));
    end
end

function [F, GF] = hessConsTest(x,component,MICS,samplingPeriod)
    [ddd, GO, HO] = TDEDiscriminant(x,MICS,samplingPeriod);
    F = GO(component);
    GF = squeeze(HO(component,:,:));
end

function [F, GF] = hessChenTest(x,PCCC,micDistances,samplingPeriod)
    [ddd,F,GF] = chenTDECriterion(x,PCCC,micDistances,samplingPeriod);
end

function checkDerivatives(fun, X0, step, name)
    % Check out the dimension
    Dimension = size(X0,1);
    % Compute the values provided by the function
    [~, Gradient] = fun(X0);
    % Check the gradient
    NGradient = zeros(size(Gradient));
    for ii=1:Dimension
        % Direction
        Direction = zeros(size(Gradient));
        Direction(ii,:) = 1;
        % Forward point
        XF = X0 + step*Direction;
        % Backward point
        XB = X0 - step*Direction;
        % Compute the gradient
        NGradient(ii,:) = (fun(XF)-fun(XB));
    end
    % Normalize
    NGradient = NGradient / (2*step);
    % Printout
    fprintf('========\n');
    fprintf('%s\n',name);
    fprintf('========\n');
    for ii=1:Dimension,
        AbsDiff = abs(Gradient(ii,:)-NGradient(ii,:));
        RelDiff = AbsDiff./abs(Gradient(ii,:));
%         fprintf('(%d) %1.5e \t%1.5e \t%1.5e \t%1.5e\n',ii,Gradient(ii),NGradient(ii),AbsDiff,RelDiff);
        fprintf('(%d) %1.5e \t%1.5e\n',ii,mean(AbsDiff),mean(RelDiff));
        figure;
%         scatter3(X0(1,:),X0(2,:),X0(3,:),RelDiff);
        global Positions;
        scatter3(Positions(:,1),Positions(:,2),Positions(:,3),RelDiff);
        title(name);
    end
end

function checkManyDerivatives(fun, X0, step, name)
    % Check out the dimension
    [Dimension, NPoints] = size(X0);
    % Compute the values provided by the function
    [ddd, Gradient] = fun(X0);
    % Check the gradient
    NGradient = zeros(size(Gradient));
%     for np=1:NPoints,
        for ii=1:Dimension
            % Direction
            Direction = zeros(size(X0));
            Direction(ii,:) = 1;
            % Forward point
            XF = X0 + step*Direction;
            % Backward point
            XB = X0 - step*Direction;
            % Compute the gradient
            NGradient(ii,:) = (fun(XF)-fun(XB));
        end
%     end
    % Normalize
    NGradient = NGradient / (2*step);
    % Printout
    fprintf('========\n');
    fprintf('%s\n',name);
    fprintf('========\n');
    for np=1:NPoints,
        for ii=1:Dimension,
            AbsDiff = abs(Gradient(ii,np)-NGradient(ii,np));
            RelDiff = AbsDiff/abs(Gradient(ii,np));
            fprintf('(%d) %1.5e \t%1.5e \t%1.5e \t%1.5e\n',ii,Gradient(ii,np),NGradient(ii,np),AbsDiff,RelDiff);
        end
    end
end