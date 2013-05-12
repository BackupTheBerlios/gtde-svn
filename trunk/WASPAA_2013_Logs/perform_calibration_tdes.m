function perform_calibration_tdes

% Load data
tmp = load('tdes_calib');
tdes = tmp.tdes;
clear tmp;

% Perform calibration
costFunction = @(x) calibration_cost(x,tdes.positions,tdes);
x0 = tdes.mics(:);
options = optimset('DerivativeCheck','off','Display','iter','GradObj','on');
xF = fminunc(costFunction,x0,options);
% fprintf('%1.4f  ',xF);
% fprintf('\n');
plot3(tdes.mics(:,1),tdes.mics(:,2),tdes.mics(:,3),'bo');
hold on
plot3(xF(1),xF(2),xF(3),'rx');
plot3(xF(4),xF(5),xF(6),'rx');
plot3(xF(7),xF(8),xF(9),'rx');
plot3(xF(10),xF(11),xF(12),'rx');
hold off

end

function [C,J] = calibration_cost(mics,positions,tdes)
% Fixed values
v = 334.2;

% Output variables
C = 0;
J = zeros(size(mics));

% Microphones
mic1 = mics(1:3);
mic2 = mics(4:6);
mic3 = mics(7:9);
mic4 = mics(10:12);

% Loop on the positions
for p = 1:size(positions,1),
    % Compute distances to microphones and normalized vectors
    v1 = positions(p,:) - mic1'; n1 = norm(v1,2); v1 = v1/n1;
    v2 = positions(p,:) - mic2'; n2 = norm(v2,2); v2 = v2/n2;
    v3 = positions(p,:) - mic3'; n3 = norm(v3,2); v3 = v3/n3;
    v4 = positions(p,:) - mic4'; n4 = norm(v4,2); v4 = v4/n4;

    % Build the true TDE
    ttde = zeros(1,size(tdes.bab{1,p},2));
    ttde(1) = (n1-n2)/v;
    ttde(2) = (n1-n3)/v;
    ttde(3) = (n1-n4)/v;
    
    % Build the jacobian matrix
    pJ = zeros(3,12);
    % V1
    pJ(1,1:3) = - v1;
    pJ(2,1:3) = - v1;
    pJ(3,1:3) = - v1;
    % V2
    pJ(1,4:6) = v2;
    % V3
    pJ(2,7:9) = v3;
    % V4
    pJ(3,10:12) = v4;
    
    % Provisional partial gradient
    auxJ = zeros(1,size(tdes.bab{1,p},2));
    % On the first group of tdes
    for ii = 1:size(tdes.bab{1,p},1),
        if isnan(tdes.bab{1,p}(ii,:)),
            continue;
        end
        C = C + norm(tdes.bab{1,p}(ii,:)-ttde,2).^2;
        auxJ = auxJ + tdes.bab{1,p}(ii,:) - v*ttde;
    end
%     % On the second group of tdes
%     for ii = 1:size(tdes.bab{1,p},1),
%         if isnan(tdes.dip{1,p}(ii,:)),
%             continue;
%         end
%         C = C + norm(tdes.dip{1,p}(ii,:)-ttde,2).^2;
%         auxJ = auxJ + tdes.dip{1,p}(ii,:) - v*ttde;
%     end
    
    % Multiply by pJ
    J = J + (auxJ * pJ)';
end

J = -2*J/v;

end