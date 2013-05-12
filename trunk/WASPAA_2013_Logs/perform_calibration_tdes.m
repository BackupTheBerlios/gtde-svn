function perform_calibration_tdes

% Load data
tmp = load('tdes_calib');
tdes = tmp.tdes;
clear tmp;

% Mics distance
d12 = norm(tdes.mics(1,:)-tdes.mics(2,:),2);
d13 = norm(tdes.mics(1,:)-tdes.mics(3,:),2);
d14 = norm(tdes.mics(1,:)-tdes.mics(4,:),2);
d23 = norm(tdes.mics(2,:)-tdes.mics(3,:),2);
d24 = norm(tdes.mics(2,:)-tdes.mics(4,:),2);
d34 = norm(tdes.mics(3,:)-tdes.mics(4,:),2);

% Perform calibration
costFunction = @(x) calibration_cost(x,tdes.positions,tdes,d12,d13,d14,d23,d24,d34);
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

function [C,J] = calibration_cost(mics,positions,tdes,d12,d13,d14,d23,d24,d34)
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

%%% Constraint

% Precomputing
v12 = mic1-mic2; n12 = norm(v12,2); v12 = v12/n12;
v13 = mic1-mic3; n13 = norm(v13,2); v13 = v13/n13;
v14 = mic1-mic4; n14 = norm(v14,2); v14 = v14/n14;
v23 = mic2-mic3; n23 = norm(v23,2); v23 = v23/n23;
v24 = mic2-mic4; n24 = norm(v24,2); v24 = v24/n24;
v34 = mic3-mic4; n34 = norm(v34,2); v34 = v34/n34;

CO = (n12-d12)^2+...
     (n13-d13)^2+...
     (n14-d14)^2+...
     (n23-d23)^2+...
     (n24-d24)^2+...
     (n34-d34)^2;
 
GO = [ (n12-d12)*v12 + (n13-d13)*v13 + (n14-d14)*v14;...
      -(n12-d12)*v12 + (n23-d23)*v23 + (n24-d24)*v24;...
      -(n13-d13)*v13 - (n23-d23)*v23 + (n34-d34)*v34;...
      -(n14-d14)*v14 - (n24-d24)*v24 - (n34-d34)*v34];

GO = 2*GO;
    
% Lambda
lambda = 10;

% Regularize
C = C + lambda*CO;
J = J + lambda*GO;

end