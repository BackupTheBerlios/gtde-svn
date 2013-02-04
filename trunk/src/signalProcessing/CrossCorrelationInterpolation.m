function [XC, XCd, XCdd, flag] = CrossCorrelationInterpolation(PCCC,Delay,SamplingPeriod,ZeroIndex)

%ccInterpolation Cross-correlation function interpolation
%
% USAGE: [XC, XCd, XCdd, flag] = CrossCorrelationInterpolation(PCCC,Delay,T[,ZI])
%
% PARAMETERS:
%  PCCC ~ cross-correlation signals of the polynomial coefficients.
%  Delay ~ delay is the point(s) in which we want to estimate the signals'
%     cross-correlation function.
%  T ~ the signals' sampling period.
%  ZI ~ the index on the PCCC representing no delay (the default value is ceil(L/2), where L is the
%     size of the coefficients cross-correlation function.)
% 
% RETURN VALUE:
%  cc,ccd,ccdd ~ The value of the cross-correlation function, its first and
%     its second derivative.
% 
% DESCRIPTION:
%     This function estimates the value of the cross-correlation function
%     at Delay from from the polynomial coefficients' cross-correlation 
%     signals PCCC.
%    
% REFERENCES:
%     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
%     estimation-based sound source localisation (gTDESSL). Research Report 
%     RR-7988, INRIA, June 2012.
%
%     see also K1, K2, InterpolateCrossCorrelation

% Copyright 2012, Xavier Alameda-Pineda
% INRIA Grenoble Rhone-Alpes
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

    % Input check
    if nargin < 3
        error('Usage: [XC, XCd, XCdd, flag] = CrossCorrelationInterpolation(PCCC,Delay,SamplingPeriod[,ZeroIndex])');
    end
    if nargin < 4
        ZeroIndex = ceil(size(PCCC,3)/2);
    end

    % Output value
    XC = zeros(size(Delay));
    XCd = zeros(size(Delay));
    XCdd = zeros(size(Delay));
    if nargout > 3,
        flag = zeros(size(Delay));
    end
    
    % Index in the cross-correlation
    ccIndex = floor(Delay/SamplingPeriod);
    % Remaining of the division (0 <= tau < SamplingPeriod)
    tau = Delay - ccIndex*SamplingPeriod;
    
    % Interpolation polynomial degrees
    D1 = size(PCCC,1)-1;
    D2 = size(PCCC,2)-1;
    
    if ndims(PCCC) == 3,
        % Which indices are inside?
        ccIndex = ccIndex + ZeroIndex;
        inside = ccIndex > 1 & ccIndex <= size(PCCC,3);
%         if sum(~inside) > 0
%             disp('inside');
%         end
        % Set the correspondent flag
        flag(~inside) = 1;
        % No points to evaluate
        if sum(inside) == 0
            return;
        end
        % Reduce the taus
        tau = tau(inside);
        % Allocate space
        cc = zeros(size(tau));
        ccd = zeros(size(tau));
        ccdd = zeros(size(tau));
%         h = figure;
%         hd = figure;
%         fk1 = figure;
%         fk2 = figure;
%         fr1 = figure;
%         fr2 = figure;
%         set(h,'Name','XC');
%         set(hd,'Name','XCD');
        for d1 = 0:D1,
            for d2 = 0:D2,
                % Value of the first interpolating function
                [k1, k1d, k1dd] = L1(d1,d2,tau,SamplingPeriod);
%                 k1d = -k1d;
%                 k1dd = k1dd*SamplingPeriod^2;
                % Value of the second interpolation function
                [k2, k2d, k2dd] = L2(d1,d2,tau,SamplingPeriod);
%                 k2d = -k2d;
%                 k2dd = k2dd*SamplingPeriod^2;
%                 if sum(isnan([k1 k2 k1d k2d k1dd k2dd])) > 0
%                    disp('isnan');
%                 end
                % First cross-correlation value
                R1 = (squeeze(PCCC(d1+1,d2+1,ccIndex(inside)))');
                % Second cross-correlation value
                R2 = (squeeze(PCCC(d1+1,d2+1,ccIndex(inside)+1))');
%                 if length(cc) ~= length(k1) || length(cc) ~= length(k2) || length(cc) ~= length(R1) || length(cc) ~= length(R2)
%                     fprintf('CC %d, K1 %d, K2 %d, R1 %d, R2 %d\n',length(cc), length(k1), length(k2), length(R1), length(R2));
%                 end
%                 fprintf('%1.10e %1.10e || %1.10e %1.10e\n',k1d,R1,k2d,R2);
                % Cross-correlation's value
                cc = cc + R1.*k1 + R2.*k2;
                % Cross-correlation's derivatives
                ccd = ccd + R1.*k1d + R2.*k2d;
                % Cross-correlation's second derivative
                ccdd = ccdd + R1.*k1dd + R2.*k2dd;
                % Some partial plotting
%                 figure(h);
%                 subplot(D1+1,D2+1,(D2+1)*d1+d2+1), plot(R1.*k1 + R2.*k2);
%                 figure(hd);
%                 subplot(D1+1,D2+1,(D2+1)*d1+d2+1), plot(R1.*k1d + R2.*k2d);
%                 figure(fk1);
%                 subplot(D1+1,D2+1,(D2+1)*d1+d2+1), plot(k1);
%                 figure(fk2);
%                 subplot(D1+1,D2+1,(D2+1)*d1+d2+1), plot(R1);
%                 figure(fr1);
%                 subplot(D1+1,D2+1,(D2+1)*d1+d2+1), plot(k2);
%                 figure(fr2);
%                 subplot(D1+1,D2+1,(D2+1)*d1+d2+1), plot(R2);
%                 h = figure;
%                 set(h,'Name',['D1 = ' num2str(d1) ' & D2 = ' num2str(d2)]);
%                 subplot(2,5,1), plot(Delay,k1,'bo'); title('K1'); 
%                 subplot(2,5,2), plot(Delay,k2,'bo'); title('K2'); 
%                 subplot(2,5,3), plot(Delay,R1.*k1,'bo'); title('R1*K1'); 
%                 subplot(2,5,4), plot(Delay,R2.*k2,'bo'); title('R2*K2'); 
%                 subplot(2,5,5), plot(Delay,cc,'bo'); title('CC');
%                 subplot(2,5,6), plot(Delay,k1d,'bo'); title('K1D'); 
%                 subplot(2,5,7), plot(Delay,k2d,'bo'); title('K2D'); 
%                 subplot(2,5,8), plot(Delay,R1.*k1d,'bo'); title('R1*K1D'); 
%                 subplot(2,5,9), plot(Delay,R2.*k2d,'bo'); title('R2*K2D'); 
%                 subplot(2,5,10), plot(Delay,ccd,'bo'); title('CCD');
            end
        end
        % Store good values
        XC(inside) = cc;
        XCd(inside) = ccd;
        XCdd(inside) = ccdd;
        
    elseif ismatrix(PCCC)
        for d1 = 0:D1,
            for d2 = 0:D2,
                % Value of the first interpolating function
                [k1, k1d, k1dd] = L1(d1,d2,tau,SamplingPeriod);
                % Value of the second interpolation function
                [k2, k2d, k2dd] = L2(d1,d2,tau,SamplingPeriod);
                % First cross-correlation value
                R1 = PCCC(d1+1,d2+1);
                % Second cross-correlation value
                R2 = PCCC(d1+1,d2+1);
                % Cross-correlation's value
                XC = XC + R1.*k1 + R2.*k2;
                % Cross-correlation's derivative
                XCd = XCd + R1.*k1d + R2.*k2d;
                % Cross-correlation's second derivative
                XCdd = XCdd + R1.*k1dd + R2.*k2dd;
            end
        end
    else
        error('Dimension of PCCC should be 2 or 3.');
    end
return

