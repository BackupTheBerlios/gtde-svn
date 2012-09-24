function [finalX f c op] = ipsolver_parallel (x, objective, constraint, op)

%ipsolver_parallel Interior-Point algorithm for several initializations at
%
% USAGE: [x f c] = ipsolver_parallel (x, objective, constraint, op)
%
% PARAMETERS:
%  x ~ is the dxn matrix of the n initializations of dimesion d.
%  objective ~ it is the handler of the objective function.
%  constraint ~ it is the handler for THE constraint, since in this version
%       of the IP algorithm we accept just one constraint.
%  op ~ structure of options: option name [default value]
%         descentdir      [newton] (only choice)
%         tolerance       [1e-6]   (if the point moves less, does not move)
%         maxiter         [10]     (maximum number of iterations)
%         verbose         [false]
%         condtest        [false]  (invertibility test)
%         eps             [1e-8]   (machine zero)
%       parameters of the IP algorithm & backtrack search see refences
%         sigmamax        [0.5]
%         etamax          [0.25]
%         mumin           [1e-9]
%         alphamax        [0.995]
%         alphamin        [1e-4]
%         beta            [0.75]
%         tau             [0.01]
%
% RETURN VALUE:
%  x ~ the matrix of output values
%  f ~ the vector of objective(x)
%  c ~ the vector of constraint(x)
% 
% DESCRIPTION:
%     This is a parallelized version of the algorithm described and
%     implemented in [Carbonetto 2008]. The code for the non-parallelized 
%     version is also available online.
%
% REFERENCES:
%     X. Alameda-Pineda and R. Horaud. Geometrically-constrained time delay
%     estimation-based sound source localisation (gTDESSL). Research Report 
%     RR-7988, INRIA, June 2012.
% 
%     P. Carbonetto, MATLAB primal-dual interior-point solver for convex 
%     programs with constraints, 2008.
%     http://www.cs.ubc.ca/ pcarbo/convexprog.html.

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

    if nargin < 4
        op = struct();
    end
    op = check_options(op);

    % INITIALIZATION.  
    % Get the number of primal variables (n), the number of constraints (m),
    % and the total number of primal-dual optimization variables (nv).
    % Initialize the Lagrange multipliers. Initialize the second-order
    % information.
    [n ni] = size(x);
    % Which point has already converged?
    finalX = zeros(size(x));
    saved = false(1,ni);
    
    %%%% FORCING to 1 constraint
    m  = 1;
    z  = ones(m,ni);
    nv = n + m;
    B  = eye(n);

    if op.verbose
        fprintf('  i f(x)       lg(mu) sigma   ||rx||  ||rc||  alpha   #ls\n');
    end

    % Repeat while the convergence criterion has not been satisfied, and
    % we haven't reached the maximum number of iterations.
    alpha = zeros(1,ni);
    ls    = zeros(1,ni);
    saved = false(1,ni);
    
    for iter = 1:op.maxiter

        % COMPUTE OBJECTIVE, GRADIENT, CONSTRAINTS, ETC.  
        % Compute the response of the objective function, the gradient of the
        % objective, the response of the inequality constraints, the Jacobian of
        % the inequality constraints, the Hessian of the Lagrangian (minus the
        % Hessian of the objective) and, optionally, the Hessian of the
        % objective.
        [c J W] = constraint(x,z);
        if strcmp(op.descentdir,'newton')
            [f g B eF] = objective(x);
        else
            error('Do not do non-newton.');
        end

        % Compute the responses of the unperturbed Karush-Kuhn-Tucker
        % optimality conditions.
        rx = g + J.*repmat(z,n,1);  % Dual residual.
        rc = c.*z;      % Complementarity.
        r0 = [rx; rc]; 
        
        % Set some parameters that affect convergence of the primal-dual
        % interior-point method.
        auxr0 = sqrt( sum( r0.^2, 1))/nv;
        eta        = min(op.etamax,auxr0);
        sigma      = min(op.sigmamax,sqrt(auxr0));
        dualitygap = -c.*z;
        mu         = max(op.mumin,sigma.*dualitygap/m);

        % CONVERGENCE CHECK.
        % (If the norm of the responses is less than the specified
        % tolerance or we could not compute the oF at the point)
        % and we have not saved the point before, save it!
        newToSave = (auxr0 < op.tolerance | eF) & ~saved;
        % Mark as saved and save
        saved(newToSave) = true;
        finalX(:,newToSave) = x(:,newToSave);
        % If everyone is done, return
        if sum(~saved) == 0
            return
        end

        % SOLUTION TO PERTURBED KKT SYSTEM.
        % Compute the search direction of x and z.
        S  = z./(c-op.eps);
        gb = g - repmat(mu,n,1).*J.*(1./(repmat(c,n,1)-op.eps));
        px = zeros(n,ni);
        for iii = 1:ni,
            if ~saved(iii)
                if op.condtest && det(B(:,:,iii) + W(:,:,iii) - J(:,iii)*S(iii)*J(:,iii)') < op.eps
                    finalX(:,iii) = x(:,iii);
                    saved(iii) = true;
                    continue
                else     
                    px(:,iii) = (B(:,:,iii) + W(:,:,iii) - J(:,iii)*S(iii)*J(:,iii)') \ (-gb(:,iii));
                end
            end
        end
        % If everyone is done, return
        if sum(~saved) == 0
            return
        end
        
        pz = -(z + mu./(c-op.eps) +  S.*sum(J.*px,1));
        

        % BACKTRACKING LINE SEARCH.
        % To ensure global convergence, execute backtracking line search to
        % determine the step length. First, we have to find the largest step
        % size which ensures that z remains feasible. Next, we perform
        % backtracking line search.
        alpha = op.alphamax*ones(1,ni);
        is    = find(z + pz < 0);
        if ~isempty(is)
            alpha(is) = op.alphamax * min(1, z(is) ./ -pz(is) );
        end

        % Compute the response of the merit function and the directional
        % gradient at the current point and search direction.
        psi  = merit(x,z,f,c,mu,op.eps);
        dpsi = gradmerit(x,z,px,pz,g,c,J,mu,op.eps);
        ls   = zeros(1,ni);
        backtrack = true(1,ni);

        while true
            % Backtrack those points which are not saved before
            backtrack = ~saved & backtrack;
            % This condition is necessary in case all the points in the
            % last inner iteration are either stopped begin backtracked or
            % saved!
            if sum(backtrack) == 0
                break;
            end
            % Compute the candidate point, the constraints, and the response of
            % the objective function and merit function at the candidate point.
            
            ls = ls + 1;
            % Initialize
            xnew = zeros(n,ni);
            znew = zeros(1,ni);
            psinew = zeros(1,ni);
            % Compute
            xnew(:,backtrack) = x(:,backtrack) + repmat(alpha(backtrack),n,1) .* px(:,backtrack);
            znew(backtrack) = z(backtrack) + alpha(backtrack) .* pz(backtrack);
            f(backtrack) = objective(xnew(:,backtrack));
            c(backtrack) = constraint(xnew(:,backtrack),[]);
            psinew(backtrack) = merit(xnew(:,backtrack),znew(backtrack),f(backtrack),c(backtrack),mu(backtrack),eps);

            % Stop backtracking search if we've found a candidate point that
            % sufficiently decreases the merit function and satisfies all the
            % constraints (and if we were backtracking it, of course)
            stopBacktracking = (~(c > 0)) & psinew < (psi + op.tau*(eta.*alpha).*dpsi);
            stopBacktracking = stopBacktracking & backtrack;
            % For those which have to stop, save the point
            x(:,stopBacktracking)     = xnew(:,stopBacktracking);
            z(stopBacktracking)     = znew(stopBacktracking);
            
            % Candidate points which do not meet our criteria, will have a
            % smaller alpha next backtracking search
            alpha(~stopBacktracking) = alpha(~stopBacktracking) * op.beta;
            % Mark the candidates which were backtracked and that have
            % a too small alpha
            newToSave = alpha<op.alphamin & backtrack;

            % Update backtrack variable, with those who do not have to be
            % backtracked any more, and check if there are still some
            % poitns to be backtracked
            backtrack = backtrack & ~stopBacktracking;
            if sum(backtrack) == 0
                break;
            end
            
            % If there are points to be backtracked, still, save those who
            % had a too small alpha
            finalX(:,newToSave) = xnew(:,newToSave);
            saved(newToSave) = true;  
            % If everyone is done, return
            if sum(~saved) == 0
                return
            end
                       
        end
    end
    % Save the final results of those initializations which did not
    % converge
    finalX(:,~saved) = x(:,~saved);
end
  
% ------------------------------------------------------------------
% Compute the response of the merit function at (x,z).
function psi = merit (~, z, f, c, mu, eps)
    psi = f - c.*z - mu.*log(c.^2.*z+eps);
end
  
% ------------------------------------------------------------------
% Compute the directional derivative of the merit function at (x,z).
function dpsi = gradmerit (~, z, px, pz, g, c, J, mu, eps)
    dpsi = sum(px.*(g - J.*repmat(z,size(J,1),1) - 2*J.*repmat(mu./(c-eps),size(J,1),1)),1) -...
        pz.*(c + mu./(z+eps));
end

function options = check_options(options)
    % Field names
    defaultValues = {...
        {'descentdir',      'newton'    },...
        {'tolerance',       1e-6        },...
        {'maxiter',         10          },...
        {'verbose',         false       },...
        {'condtest',        false       },...
        {'eps',             1e-8        },...
        {'sigmamax',        0.5         },...
        {'etamax',          0.25        },...
        {'mumin',           1e-9        },...
        {'alphamax',        0.995       },...
        {'alphamin',        1e-4        },...
        {'beta',            0.75        },...
        {'tau',             0.01        }...
        };
    
    % Loop on it
    for no = 1:numel(defaultValues),
        if ~isfield(options,defaultValues{no}{1})
            options.(defaultValues{no}{1}) = defaultValues{no}{2};
        end
    end
end