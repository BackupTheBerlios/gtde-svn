function testIPSolver

%Test the parallel version of the Interior-Point algorithm.
%
%  see also, ipsolver_parallel

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

    close all;
    % Generating points
    x = 5*(rand(2,10000)-0.5);
%     x = [1 3]';
    f = objectiveFunction(x);
    c = constraintFunction(x);
    % Feasible points
    indices = c < 0;
    % Plot
    figure
    plot3(x(1,indices),x(2,indices),f(indices),'go');
    hold on
    plot3(x(1,~indices),x(2,~indices),f(~indices),'rx');
    % Objective and constraint function
    oF = @(y) objectiveFunction(y);
    cF = @(y,z) constraintFunction(y,z);
    % Performing optimization
    fesibleX = x(:,indices);
    op.verbose = false;
    [xP fP] = ipsolver_parallel(fesibleX,oF,cF,op);
    % Ploting
    plot3(xP(1,:),xP(2,:),fP,'bd');
end

function [f Gf Hf] = objectiveFunction(x)
    % If just one sample
    if size(x,1) == 1
        % Function
        f = sum(x.^2);
        % Gradient
        if nargout > 1
            Gf = 2*x';
        end
        % Hessian
        if nargout > 2
            Hf = 2*eye(2);
        end
    else
        % Objective function
        f = sum(x.^2,1);
        % Fradient
        if nargout > 1
            Gf = 2*x;
        end
        % Hessian
        if nargout > 2
            Hf = zeros(2,2,size(x,2));
            Hf(1,1,:) = 2;
            Hf(2,2,:) = 2;
        end
    end
end

function [c Gc Hc] = constraintFunction(x,z)
    % If just one sample
    if size(x,1) == 1
        % Constraint
        c = x(1)^2 + 1 - x(2);
        % Gradient
        if nargout > 1
            Gc = 2*x';
            Gc(2) = -1;
        end
        % Hessian
        if nargout > 2
            Hc = zeros(2);
            Hc(1,1) = 2*z(1);
        end
    else
        % Constraint
        c = x(1,:).^2 + 1 - x(2,:);
        % Gradient
        if nargout > 1
            Gc = 2*x;
            Gc(2,:) = -1;
        end
        % Hessian
        if nargout > 2
            Hc = zeros(2,2,size(x,2));
            Hc(1,1,:) = 2;
            Hc = z(1)*Hc;
        end
    end
end