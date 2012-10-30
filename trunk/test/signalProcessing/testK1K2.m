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



% Polynomial degrees
for D1 = 0:3,
    for D2 = 0:3,

        % Taus
        tauinc = 0.01;
        % Sampling period
        T = 0.5;
        tau = -5*T:tauinc:5*T;
        tau = tau + 1e-7*rand(size(tau));

        % Compute K1, K2 and their derivatives with respect to tau
        [k1 k1d k1dd] = K1(D1,D2,tau,T);
        [k2 k2d k2dd] = K2(D1,D2,tau,T);

        % Derivative step
        h = 0.9*tauinc;

        % Plot k1
        figure
%         subplot(2,3,1)
        plot(tau,k1,'bo');
        title(['K1 widht D1=' num2str(D1) ' D2=' num2str(D2)]);
        hold on
%         for ii = 1:numel(tau),
%             plot([tau(ii);tau(ii)+h],[k1(ii);k1(ii) + k1d(ii)*h],'r');
%         end
%         for ii = 1:numel(tau),
%             plot([tau(ii);tau(ii)+h],[k1(ii);k1(ii) + k1d(ii)*h + k1dd(ii)*h^2/2],'g');
%         end
%         hold off
%         legend('k1','k1d');
        
%         subplot(2,3,2)
%         hold on
        plot(tau,k1d,'rx');
%         k1de = (k1(2:end)-k1(1:end-1))/tauinc;
%         plot(tau(1:end-1),k1de,'bo');
%         hold off
        
%         subplot(2,3,3)
        plot(tau,k1dd,'gd');
%         hold on
%         k1dde = (k1d(2:end)-k1d(1:end-1))/tauinc;
%         plot(tau(1:end-1),k1dde,'bo');
        hold off
        

        % Plot k2
%         subplot(2,3,4)
        figure
        plot(tau,k2,'bo');
        title(['K2 widht D1=' num2str(D1) ' D2=' num2str(D2)]);
        hold on
%         for ii = 1:numel(tau),
%             plot([tau(ii);tau(ii)+h],[k2(ii);k2(ii) + k2d(ii)*h],'r');
%         end
%         for ii = 1:numel(tau),
%             plot([tau(ii);tau(ii)+h],[k2(ii);k2(ii) + k2d(ii)*h + k2dd(ii)*h^2/2],'g');
%         end
%         hold off
%         legend('k2','k2d');
        
%         subplot(2,3,5)
%         hold on
        plot(tau,k2d,'rx');
%         k2de = (k2(2:end)-k2(1:end-1))/tauinc;
%         plot(tau(1:end-1),k2de,'bo');
%         hold off
        
%         subplot(2,3,6)
        plot(tau,k2dd,'gd');
%         hold on
%         k2dde = (k2d(2:end)-k2d(1:end-1))/tauinc;
%         plot(tau(1:end-1),k2dde,'bo');
        hold off
    end
end