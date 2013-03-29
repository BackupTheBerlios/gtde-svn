function [SQPLABOutput, SQPLABOptions] = gTDE_SQPLAB(PCCC,MICS,T,X0,Energies)

%Interface for the gtde software with SQPLAB
% 
% USAGE [foundTDE] = gTDE_SQPLAB(PCCC,MICS,T,X0)
% 
% DESCRIPTION
%     This function is just an interface for the gtde program to call the 
%     sqplab minimizer.
%
% SEE ALSO
%     gTDE_simul, sqplab
%     
% REFERENCES:
%     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
%     estimation-based sound source localisation (gTDESSL). Research Report 
%     RR-7988, INRIA, June 2012.

% Copyright 2013, Xavier Alameda-Pineda
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

    %%% Simulation function
    simulationFunction = @(indic,x,varargin) gTDE_simul(PCCC,MICS,T,Energies,indic,x,varargin);

    %%% Hard-coded options
    options.algo_method        = 'quasi-Newton';
    options.algo_globalization = 'unit stepsize';
    
    options.tol(1)  = 1.e-3;  % tolerance on the gradient of the Lagrangian
    options.tol(2)  = 1.e-3;  % tolerance on the feasibility
    options.tol(3)  = 1.e-3;  % tolerance on the complementarity

    options.dxmin   = sqrt(eps);  % minimum size of a step
    options.inf     = inf;   % infinity value for the Moduloptmatlab collection

    options.miter   = 10;   % max iterations
    options.msimul  = 1000;   % max simulations

    options.fout    = 1;      % print file identifier
    options.verbose = 0;      % verbosity level

    %%% Lower and upper bounds
    ub = TDEmax(MICS)/T;
    ub = ub(1:3);
    lb = -ub;
    % Also for the constraint
    lb = cat(1,lb,-inf);
    ub = cat(1,ub,0);
    
    %%% Output
    XF = zeros(size(X0));
    fs = zeros(size(XF,1),1);
    cs = zeros(size(XF,1),1);
    flag = zeros(size(XF,1),1);
    
    %%% Loop for all X0's
    for ii = 1:size(X0,2),
        % Call sqplab
        [XF(:,ii), ~, info] = sqplab(simulationFunction,X0(:,ii),[],lb,ub,options);
        % Evaluate at XF(ii,:)
        fs(ii) = gTDECriterionMod(XF(:,ii),PCCC,MICS,T,Energies);
        cs(ii) = TDEDiscriminant(XF(:,ii),MICS,T);
        % Flag
        flag(ii) = info.flag;
    end
    
    %%% Save
    SQPLABOutput.x = XF;
    SQPLABOutput.f = fs;
    SQPLABOutput.c = cs;
    SQPLABOutput.flag = flag;
        
    SQPLABOptions = 0;
end