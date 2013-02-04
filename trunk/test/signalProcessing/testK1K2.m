function testK1K2
%test the K1 and K2 functions

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

    close all; clear;
    
    step = 1e-6;

    % Polynomial degrees
    for D1 = 0:3,
        for D2 = 0:3,

            % Taus
            tauinc = 0.01;
            % Sampling period
            T = 1/48000;
%             tau = -T:tauinc*T:T;
%             tau = tau + 1e-7*rand(size(tau));
            tau = -T*10.^(-50:0);
%             tau = cat(2,-tau,tau);

            % Compute K1, K2 and their derivatives with respect to tau
            [k1, k1d, k1dd] = K1(D1,D2,tau,T);
            [k2, k2d, k2dd] = K2(D1,D2,tau,T);

            tauPlot = tau;
            
            % Derivative step
            h = 0.9*tauinc;
            
%             if  D1==1 & D2==1,
%                 disp('stop');
%             end
                        
            nameCheck = ['K1 widht D1=' num2str(D1) ' D2=' num2str(D2)];

            % Plot k1
            figure
    %         subplot(2,3,1)
            plot(tauPlot,k1,'bo');
            title(nameCheck);
            hold on
%             for ii = 1:numel(tau),
%                 plot([tau(ii);tau(ii)+h],[k1(ii);k1(ii) + k1d(ii)*h],'r');
%             end
%             for ii = 1:numel(tau),
%                 plot([tau(ii);tau(ii)+h],[k1(ii);k1(ii) + k1d(ii)*h + k1dd(ii)*h^2/2],'g');
%             end
%             hold off
%             legend('k1','k1d');
% 
%             subplot(2,3,2)
%             hold on
            plot(tauPlot,k1d,'rx');
%             k1de = (k1(2:end)-k1(1:end-1))/tauinc;
%             plot(tau(1:end-1),k1de,'bo');
%             hold off
% 
%             subplot(2,3,3)
            plot(tauPlot,k1dd,'gd');
%             hold on
%             k1dde = (k1d(2:end)-k1d(1:end-1))/tauinc;
%             plot(tau(1:end-1),k1dde,'bo');
            hold off
            
%             nameCheck = ['Derivative of K1 widht D1=' num2str(D1) ' D2=' num2str(D2)];
%             myFun = @(x) K1(D1,D2,x,T);
%             checkDerivatives(myFun,tau,step,nameCheck);
%             nameCheck = ['Second derivative of K1 widht D1=' num2str(D1) ' D2=' num2str(D2)];
%             myFun = @(x) K1D(D1,D2,x,T);
%             checkDerivatives(myFun,tau,step,nameCheck);


            % Plot k2
%             subplot(2,3,4)
            nameCheck = ['K2 widht D1=' num2str(D1) ' D2=' num2str(D2)];

            figure
            plot(tauPlot,k2,'bo');
            title(nameCheck);
            hold on
%             for ii = 1:numel(tau),
%                 plot([tau(ii);tau(ii)+h],[k2(ii);k2(ii) + k2d(ii)*h],'r');
%             end
%             for ii = 1:numel(tau),
%                 plot([tau(ii);tau(ii)+h],[k2(ii);k2(ii) + k2d(ii)*h + k2dd(ii)*h^2/2],'g');
%             end
%             hold off
%             legend('k2','k2d');
% 
%             subplot(2,3,5)
%             hold on
            plot(tauPlot,k2d,'rx');
%             k2de = (k2(2:end)-k2(1:end-1))/tauinc;
%             plot(tau(1:end-1),k2de,'bo');
%             hold off
% 
%             subplot(2,3,6)
            plot(tauPlot,k2dd,'gd');
%             hold on
%             k2dde = (k2d(2:end)-k2d(1:end-1))/tauinc;
%             plot(tau(1:end-1),k2dde,'bo');
            hold off
            
%             nameCheck = ['Derivative of K2 widht D1=' num2str(D1) ' D2=' num2str(D2)];
%             myFun = @(x) K2(D1,D2,x,T);
%             checkDerivatives(myFun,tau,step,nameCheck);
%             nameCheck = ['Second derivative of K2 widht D1=' num2str(D1) ' D2=' num2str(D2)];
%             myFun = @(x) K2D(D1,D2,x,T);
%             checkDerivatives(myFun,tau,step,nameCheck);
        end
    end

end

function [k1d, k1dd] = K1D(D1,D2,tau,T)
    [~, k1d, k1dd] = K1(D1,D2,tau,T);
end

function [k2d, k2dd] = K2D(D1,D2,tau,T)
    [~, k2d, k2dd] = K2(D1,D2,tau,T);
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
    fprintf('%s: ',name);
    fprintf('========\n');
    for ii=1:Dimension,
        AbsDiff = abs(Gradient(ii,:)-NGradient(ii,:));
        RelDiff = zeros(size(AbsDiff));
        indices = Gradient~=0 ;%& AbsDiff > 10e-10;
        RelDiff(indices) = AbsDiff(indices)./abs(Gradient(ii,indices));
%         fprintf('(%d) %1.5e \t%1.5e \t%1.5e \t%1.5e\n',ii,Gradient(ii),NGradient(ii),AbsDiff,RelDiff);
        fprintf('Abs ~ (%1.5e,%1.5e)   Rel: (%d;%1.5e)\n',min(AbsDiff),max(AbsDiff),sum(indices),max(RelDiff));
%         figure;
%         plot(log10(AbsDiff));
%         scatter3(X0(1,:),X0(2,:),X0(3,:),log(2000*(RelDiff)+2));
%         global Positions;
%         scatter3(Positions(:,1),Positions(:,2),Positions(:,3),20000*(RelDiff+eps));
%         title(name);
    end
    % New figure
%     figure;
%     scatter3(X0(1,:),X0(2,:),X0(3,:),3,Value);
%     hold on;
%     for np=1:size(X0,2),
%     plot3([X0(1,np),X0(1,np)+5000*Gradient(1,np)],...
%           [X0(2,np),X0(2,np)+5000*Gradient(2,np)],...
%           [X0(3,np),X0(3,np)+5000*Gradient(3,np)],'km-');
%     end
%     hold off
end