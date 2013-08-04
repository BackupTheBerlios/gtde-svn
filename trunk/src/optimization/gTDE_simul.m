    function [outdic, varargout] = gTDE_simul(PCCC,MICS,T,Energies,indic,x,varargin)

%Interface for the gtde software with SQPLAB
% 
% USAGE [outdic, varargout] = gTDE_simul(indic,x,varargin)
% 
% DESCRIPTION
%     This function is just an interface for the gtde program to call the 
%     sqplab minimizer. Please download sqplab and type 'help sqplab' to 
%     understand the syntax of the present function.
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

    %%% INPUT CHECK
    % Not enough arguments
    if nargin < 2 || nargout < 1
        error('Usage: [outdic, varargout] = gTDE_simul(indic,x,varargin)');
    end
    % Check correspondence between indic and varargin
    errorFlag = false;
    switch indic
        case 1,
            if numel(varargin) < 1, errorFlag = true; end
        case {5,6},
            if numel(varargin) < 1, errorFlag = true; end
        case {2,3,4}
        case {7,11,12,13,14,15,16}
            error('gTDE_SQPLAB was asked to compute state constraints (or derivatives) which does not have.');
        otherwise,
            errorFlag = true;
    end
    if errorFlag,
        error('Indic variable does not correspond to input/output arguments. Check SQPlab documentation.');
    end
    
    %%% Define the functions
    costFunction = @(x) gTDECriterionMod(x,PCCC,MICS,T,Energies);
    inequalityConstraint = @(x) TDEDiscriminant(x,MICS,T);
    
    % Output variables
    outdic = 0;

    %%% CASE-DEPENDENT CALLS
    switch indic
        case 1,
            % Do nothing
        case 2,
            % Compute the cost function and all the constraints at x
            % Cost function
            f = costFunction(x);
            % Inequality constraints
            ci = inequalityConstraint(x);
            % The rest is empty
            varargout = {f,-ci,[],[],[],[],[]};
            % Checking f
            if f < 0
                outdic = 2;
            end
        case 3,
            % Compute the gradient of cost function and of all the constraints at x
            % Cost function
            [~,g,~,flag] = costFunction(x);
            % Inequality constraints
            [~,ai] = inequalityConstraint(x);
            % The rest is empty
            varargout = {[],[],[],[],g,-ai',[]};
            % Check
            if flag > 0
                outdic = 2;
            end
        case 4,
            % Compute the cost function and all the constraints at x
            % Cost function
            [f,g,~,flag] = costFunction(x);
            % Inequality constraints
            [ci,ai] = inequalityConstraint(x);
            % The rest is empty
            varargout = {f,-ci,[],[],g,-ai',[]};
            % Check
            if flag > 0
                outdic = 2;
            end
        case 5,
            % Compute the hessian of the lagrangian
            % Hessian of the cost function
            [~,~,H,flag] = costFunction(x);
            % Hessian of the inequality constraint
            [~,~,W] = inequalityConstraint(x);
            % Total
            varargout = {H - varargin{1}{1}(4)*W};
            % Check
            if flag > 0
                outdic = 2;
            end
    end
    
end