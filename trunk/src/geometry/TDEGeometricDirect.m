function [TDEs Jacobians] = TDEGeometricDirect(X, MICS, C)

% TDEGeometricDirect geometric model for TDE
% 
%   TDEGeometricDirect(X, MICS) return the ITD values that a source on X 
%   generates with the microphones place at MICS. X and MICS should have 
%   the same number of columns. MICS should have at least two rows.
%   
%   TDEGeometricDirect(X, MICS, C) returns the same as ITDFunc(X, MICS) 
%   with the specified sound propagation speed. The default is C = 343.2
%   
%   [TDEs Jacobians] = TDEGeometricDirect(X, MICS) returns the Jacobian 
%   of each of the ITD functions at each potision X.
%   
%   see also DMIC

    %%% Check Input
    % If less than two arguments
    if(nargin < 2)
        error('Usage: [TDEs Jacobians] = TDEGeometricDirect(X, MICS[, C])');
    end

    % Default argument for sound propagation speed
    if(nargin < 3)
        C = 343.2;
    end

    % Check the dimension
    if size(X,2)~=size(MICS,2)
        error('X and MICS should have the same number of columns.');
    end

    % Number of microphones and points
    NMics = size(MICS,1);
    NPoints = size(X,1);

    % Number of microphones
    if NMics < 2
        error('At least two microphones are needed to compute TDEs.');
    end

    % Number of points
    if NPoints < 1
        error('At least one point is needed to compute TDEs.');
    end

    % Useful parameters
    Dimension = size(X,2);% = size(MICS,2)
    NTDEs = NMics*(NMics-1)/2;

    %%% Compute the ITD values

    % Compute the distance to each microphone
    DMICS = zeros(NPoints,NMics);
    for mic=1:NMics,
        DMICS(:,mic) = DMIC(X,MICS(mic,:));
    end

    TDEs = zeros(NPoints,NTDEs);

    % Compute the time difference of each microphone pair
    for mic1 = 1:NMics,
        for mic2 = mic1+1:NMics,
            % Compute the index corresponding to the mic1,mic2
            si = TDEis(mic1,mic2,NMics);
            TDEs(:,si) = (DMICS(:,mic2) - DMICS(:,mic1))/C;
        end
    end

    %%% Compute the jacobian if desired

    % If desired, compute the Jacobian at each point
    if nargout > 1
        % Initialize the Jacobians varaible
        Jacobians = zeros(NPoints,NTDEs,Dimension);

        % Compute the directions from the sources to the microphones
        DIRMICS = zeros(NPoints,NMics,Dimension);
        for mic = 1:NMics
            DIRMICS(:,mic,:) = bsxfun(@minus,X,MICS(mic,:))./repmat(DMICS(:,mic),1,Dimension);
        end

        % Compute the Jacobians
        for mic1 = 1:NMics,
            for mic2 = mic1+1:NMics,
                % Compute the index corresponding to the mic1,mic2
                si = TDEis(mic1,mic2,NMics);
                Jacobians(:,si,:) = (DIRMICS(:,mic2,:) - DIRMICS(:,mic1,:))/C;
            end
        end
    end
    
end % function [TDEs Jacobians] = TDEGeometricDirect(X, MICS, C)