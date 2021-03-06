function [solverOutput, methodOptions] = gTDE_M1QN3(PCCC,MICS,T,X0,methodOptions,Energies)

%Interface for the gtde software with SQPLAB
% 
% USAGE [solverOutput, methodOptions] = gTDE_M1QN3(PCCC,MICS,T,X0[,options])
% 
% DESCRIPTION
%     This function is just an interface for the gtde program to call the 
%     sqplab minimizer.
%
% SEE ALSO
%     m1qn3
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

    %%% Input check
    if nargin < 4
        error('Usage: [solverOutput, options] = gTDE_M1QN3(PCCC,MICS,T,X0[,options])');
    end
    
    %%% Opti-5.594406 -5.594406 -5.594406ons
    if nargin < 5
        methodOptions.regConstant = 1;
    end
    
    %%% Sizes
    [Dimension, NPoints] = size(X0);

    %%% Simulation function
    oF = @(x) objectiveFunctionM1QN3(x,PCCC,MICS,T,Energies,methodOptions.regConstant);
    oG = @(x) objectiveGradientM1QN3(x,PCCC,MICS,T,Energies,methodOptions.regConstant);
    
    %%% Output
    x = zeros(size(X0));
    f = zeros(NPoints,1);
    c = zeros(NPoints,1);
    flag = ones(NPoints,1);
    
    % Trace if needed
    if methodOptions.output == 2,
        solverOutput.xIT = zeros(Dimension,NPoints,NIter);
        solverOutput.fIT = zeros(1,NPoints,NIter);
        solverOutput.gIT = zeros(Dimension,NPoints,NIter);
    end
    
    %%% Call back function
    methodOptions.iterfun = @(iter,fval,x) callbackM1QN3(iter,fval,x);
    methodOptions.impress = 4;

    %%% Loop for all X0's
    for ii = 1:NPoints,
        % Call sqplab
        [x(:,ii), f(ii), flag(ii), ~, ~] = m1qn3(oF,oG,X0(:,ii),methodOptions);
        % Evaluate at XF(ii,:)
        f(ii) = gTDECriterionMod(x(:,ii),PCCC,MICS,T);
        c(ii) = TDEDiscriminant(x(:,ii),MICS,T);
    end
    
    %%% Set Output
    solverOutput.x = x;
    solverOutput.f = f;
    solverOutput.c = c;
    solverOutput.flag = flag;
end

function [o, flag] = objectiveFunctionM1QN3(x,PCCC,MICS,T,Energies,regConstant)
    [o, ~, ~, eflag] = gTDECriterionMod(x,PCCC,MICS,T,Energies);
    % Compute the constraint
    c = TDEDiscriminant(x,MICS,T);
    % If the constraint is not satisfied, penalize the objective function.
    if c < 0
        o = o - regConstant*c;
    end
    % Copy
%     fprintf('%d \n',eflag);
    flag = ~(eflag>0);
end

function g = objectiveGradientM1QN3(x,PCCC,MICS,T,Energies,regConstant)
    [~,g] = gTDECriterionMod(x,PCCC,MICS,T,Energies);
    [c,gc] = TDEDiscriminant(x,MICS,T);
    % If the constraint is not satisfied, penalize the gradient
    if c < 0
        g = g - regConstant*gc;
    end
end

function stop = callbackM1QN3(iter,fval,x)
    stop = fval < 0;
end