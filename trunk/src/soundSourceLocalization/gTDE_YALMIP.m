function [YALMIPBOutput, YALMIPOptions] = gTDE_YALMIP(PCCC,MICS,T,Energies)

%Interface for the gtde software with SQPLAB
% 
% USAGE [YALMIPBOutput, YALMIPOptions] = gTDE_YALMIP(PCCC,MICS,T,X0,Energies)
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

    warning('off','YALMIP:strict');

%     x1 = sdpvar(1,1);
%     x2 = sdpvar(1,1);
%     x3 = sdpvar(1,1);
    x = sdpvar(3,1);

%     p = -2*x1+x2-x3;
    objFunction = gTDECriterionMod(x,PCCC,MICS,T,Energies);

%     F = [x1*(4*x1-4*x2+4*x3-20)+x2*(2*x2-2*x3+9)+x3*(2*x3-13)+24>=0,
%      4-(x1+x2+x3)>=0,
%       6-(3*x2+x3)>=0,
%                x1>=0,
%              2-x1>=0,
%                x2>=0,
%                x3>=0,
%              3-x3>=0]

    constraints = TDEDiscriminant(x,MICS,samplingPeriod)>=0;

    options = sdpsettings('verbose',1,'solver','bmibnb');
    solvesdp(constraints,objFunction,options);

    
%     %%% Save
%     SQPLABOutput.x = XF;
%     SQPLABOutput.f = fs;
%     SQPLABOutput.c = cs;
%     SQPLABOutput.flag = flag;
%         
%     SQPLABOptions = 0;
end