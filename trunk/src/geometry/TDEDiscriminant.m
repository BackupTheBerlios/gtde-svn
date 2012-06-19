function [Delta GDelta HDelta] = TDEDiscriminant(TDEs,MICS,samplingPeriod,C)

%TDEDiscriminant implements the geometric constraint for geometric TDE.
%
% USAGE: [Delta GDelta HDelta] = TDEDiscriminant(TDEs,MICS,[samplingPeriod,C])
%
% PARAMETERS:
%  TDEs ~ Point at which the discriminant should be evaluated.
%  MICS ~ Positions of the microphones.
%  samplingPeriod ~ Sampling period of the discrete signals.
%  C ~ sound propagation speed.
%
% RETURN VALUE:
%  Delta,GDelta,HDelta ~ Value of the function, the first and the second 
%     derivatives.
% 
% DESCRIPTION:
%     [Delta GDelta HDelta] = TDEDiscriminant(TDEs,MICS) returns the value 
%     of such discriminant. This value is positive if the set of TDEs is 
%     consistent (i.e. correspond to a position in the space) and negative 
%     if it is not consistent. The TDEs are expected in seconds and the 
%     microphones' positions in meters. It also returns the gradient and
%     the Hessian.
% 
%     NOTE: The reference microphone for computing the TDEs is the first one
%     in the variables MICS.
% 
%     Delta = TDEDiscriminant(TDEs,MICS,samplingPeriod,C) allows the specification 
%     of the sampling Period and of the sound speed. The default is 
%     samplingPeriod = 0 and C = 343.2. Here the units have to agree. The
%     samplingPeriod is used to scale the optimization problem. If not
%     specified, the problem is not scaled.
%
% REFERENCES:
%     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
%     estimation-based sound source localisation (gTDESSL). Research Report 
%     RR-7988, INRIA, June 2012.
%
%   see TDEGeometricDirect TDEGeometricInverse

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


    %%% Input check
    % Default sound speed
    if nargin < 4
         C = 343.2;
        % No sampling period
        if nargin < 3
           samplingPeriod = 0;
            % Error, bad usage
            if nargin < 2
                error('Usage: [Delta JDelta] = TDEDiscriminant(TDEs,MICS[, samplingPeriod, C])');
            end
        else
            if samplingPeriod > 0
                TDEs = TDEs*samplingPeriod;
            end
        end
    end
    
    % Number of MICS and TDEs
    NMics = size(MICS,1);
    [Dimension NTDEs] = size(TDEs);
    if Dimension ~= NMics-1
        error('The number of TDEs should be NMics-1.');
    end
    
    %%% Useful variables
    % Compute the system matrix
    M = 2*(MICS-repmat(MICS(1,:),NMics,1));
    % The first row does not belong to the matrix
    M(1,:) = [];
    % Its inverse
    MInv = pinv(M);
    % A vector
    A = zeros(NMics-1,NTDEs);
    for ii=1:NMics-1
        A(ii,:) = -2*C*sum(repmat(MInv(ii,:),NTDEs,1).*(TDEs'),2);
    end
    % K and B vectors
    K = sum( MICS.^2, 2) - sum(MICS(1,:).^2);
    K(1) = [];
    % First B
    B = repmat(MInv*K,1,NTDEs);
    for ii=1:NMics-1
        B(ii,:) = B(ii,:) - C^2*(sum(repmat(MInv(ii,:),NTDEs,1).*((TDEs.^2)'),2)');
    end
    
    BM1 = B - repmat(MICS(1,:)',1,NTDEs);
    
    %%% Computation
    % Discriminant value
    ABM1 = sum( A.*BM1 ,1);
    BM12 = sum(BM1.^2,1);
    A2 = sum(A.^2,1);
    Delta = ABM1.^2 - BM12.*(A2-1);
    
    % If asked, compute the gradient
    if nargout > 1
        % Jacobian of the A and B vectors
        JA = -2*C*repmat(MInv,[1,1,NTDEs]);
        JB = zeros(NMics-1,NMics-1,NTDEs);
        for ii=1:NMics-1
            JB(:,ii,:) = -2*C^2*repmat(MInv(:,ii),[1,NTDEs]).*repmat(TDEs(ii,:),[NMics-1,1]);
        end
        % Gradient
        GDelta = zeros(NMics-1,NTDEs);
        % Loop in the dimension
        for ii=1:NMics-1
            GDelta(ii,:) = 2*ABM1.*( sum(squeeze(JA(:,ii,:)).*BM1,1) + sum(squeeze(JB(:,ii,:)).*A,1) ) -...
                     2*(A2-1).*sum(squeeze(JB(:,ii,:)).*BM1,1) -...
                     2*BM12.*sum(squeeze(JA(:,ii,:)).*A,1);
        end
        % Scale if needed
        if samplingPeriod > 0
            GDelta = GDelta * samplingPeriod;
        end
    end
    
    % If asked, compute the Hessian
    if nargout > 2
        % Declare and build matrices C, D and E
        CMat = -2*C^2*repmat(MInv,[1,1,NTDEs]);
        DMat = zeros(Dimension,Dimension,NTDEs);
        EMat = zeros(Dimension,Dimension,NTDEs);
        for mic = 1:Dimension,
            DMat(mic,mic,:) = sum( squeeze(CMat(:,mic,:)).*A, 1);
            EMat(mic,mic,:) = sum( squeeze(CMat(:,mic,:)).*BM1, 1);
        end
        HDelta = zeros(Dimension,Dimension,NTDEs);
        % Compute the hessian
        for d1 = 1:Dimension,
            for d2 = 1:Dimension,
                PartialHess = - (A2-1) .* squeeze(EMat(d1,d2,:))';
                PartialHess = PartialHess - 2*sum( squeeze(JB(:,d1,:)).*BM1 ,1).*...
                    sum( squeeze(JA(:,d2,:)).*A ,1);
                PartialHess = PartialHess - ...
                    (A2-1).*squeeze( sum( JB(:,d1,:).*JB(:,d2,:) ,1) )';
                PartialHess = PartialHess - ...
                    2*sum( squeeze(JA(:,d1,:)).*A ,1).*...
                    sum( squeeze(JB(:,d2,:)).*BM1 ,1);
                PartialHess = PartialHess - ...
                    BM12.*squeeze( sum( JA(:,d1,:).*JA(:,d2,:) ,1) )';
                PartialHess = PartialHess + ...
                    sum( squeeze(JA(:,d1,:)).*BM1 + squeeze(JB(:,d1,:)).*A ,1).*...
                    sum( squeeze(JA(:,d2,:)).*BM1 + squeeze(JB(:,d2,:)).*A ,1);
                PartialHess = PartialHess + ...
                    ABM1.*squeeze( sum( JA(:,d1,:).*JB(:,d2,:),1) + DMat(d1,d2,:) +...
                    sum( JB(:,d1,:).*JA(:,d2,:),1) )';
                HDelta(d1,d2,:) = PartialHess;
            end
        end
        HDelta = 2*HDelta;  
        % Scale if needed
        if samplingPeriod > 0
            HDelta = HDelta * (samplingPeriod^2);
        end
    end
end