% X = IPSOLVER(X0,OBJ,GRAD,CONSTR,JACOBIAN,DESCENTDIR,TOL,MAXITER,VERBOSE)
% is a simple yet reasonably robust implementation of a primal-dual
% interior-point solver for convex programs with convex inequality
% constraints (it does not handle equality constraints). Precisely speaking,
% it will compute the solution to the following optimization problem:
%
%     minimize    f(x)
%     subject to  c(x) < 0,
%
% where f(x) is a convex objective and c(x) is a vector-valued function with
% outputs that are convex in x. There are many, many optimization problems
% that can be framed in this form (see the book "Convex Optimization" by
% Boyd and Vandenberghe for a good start). This code is mostly based on
% the descriptions provided in this reference:
%
%    Paul Armand, Jean C. Gilbert, Sophie Jan-Jegou. A Feasible BFGS
%    Interior Point Algorithm for Solving Convex Minimization Problems. 
%    SIAM Journal on Optimization, Vol. 11, No. 1, pp. 199-222.
%
% However, to understand what is going on you will need to read up on
% interior-point methods for constrained optimization. A good starting
% point is the book of Boyd and Vandenberghe.
%
% The input X0 is the initial point for the solver. It must be an n x 1
% matrix, where n is the number of (primal) optimization variables.
% DESCENTDIR must be either: 'newton' for the Newton search direction,
% 'bfgs' for the quasi-Newton search direction with the
% Broyden-Fletcher-Goldfarb-Shanno (BFGS) update, or 'steepest' for the
% steepest descent direction. The steepest descent direction is often quite
% bad, and the solver may fail to converge to the solution if you take this
% option. For the Newton direction, you must be able to compute the the
% Hessian of the objective. Also note that we form a quasi-Newton
% approximation to the objective, not to the Lagrangian (as is usually
% done). This means that you will always have to provide second-order
% information about the inequality constraint functions.
% 
% TOL is the tolerance of the convergence criterion; it determines when the
% solver should stop. MAXITER is the maximum number of iterations. And the
% final input, VERBOSE, must be set to true or false depending on whether
% you would like to see the progress of the solver.
%
% The inputs OBJ, GRAD, CONSTR and JACOBIAN must all be function handles. If
% you don't know what function handles are, type HELP FUNCTION_HANDLE in
% MATLAB.
% 
%    * OBJ must be a handle to a function that takes 1 input, the vector
%      of optimization variables, and returns the value of the function at
%      the given point. The function definition should look something like F
%      = OBJECTIVE(X).
%
%    * GRAD is a pointer to a function of the form G = GRADIENT(X), where
%      G is the n x 1 gradient vector of the objective, or [G H] =
%      GRADIENT(X) if the Newton step is used, in which case H is the n x n
%      Hessian of the objective.
%
%    * CONSTR is a handle to a function of the form C = CONSTRAINTS(X),
%      where C is the m x 1 vector of constraint responses at X.
%
%    * JACOBIAN is a handle to a function of the form [J W] =
%      JACOBIAN(X,Z). The inputs are the primal variables X and the m x 1
%      vector of dual variables Z. The return values are the m x n
%      Jacobian matrix (containing the first-order partial derivatives
%      of the inequality constraint functions), and W is the n x n
%      Hessian of the Lagrangian (minus the Hessian of the objective), 
%      which is basically equal to
%
%          W = z(1)*W1 + z(2)*W2 + ... + z(m)*Wm,
%
%      where Wi is the Hessian of the ith constraint.
%
% If you set VERBOSE to true, then at each iteration the solver will output
% the following information (from left to right): 1. the iteration number,
% 2. the value of the objective, 3. the barrier parameter mu, 4. the
% centering parameter sigma, 4. the residuals of the perturbed
% Karush-Kuhn-Tucker system (rx, rc), 5. the step size, and the number of
% iterations in the line search before we found a suitable descent step.
%
% If your optimization problem is large (i.e. it involves a lot of
% optimization variables or inequality constraints) it might speed up the
% solver to output sparse matrices. (Type HELP SPARSE in the MATLAB console
% for more information on sparse matrices.)
%
% As a final note, the interior-point solver may not work very well if your
% problem is very poorly scaled (i.e. the Hessian of the objective or the
% Hessian of one of the constraint functions is poorly conditioned). It
% is up to you to make sure you look at the conditioning of your problem.
%
%                                         Peter Carbonetto
%                                         Dept. of Computer Science
%                                         University of British Columbia
%                                         Copyright 2008
%
function [finalX f c] = ipsolver_parallelEUSIPCO (x, objective, constraint, ...
		       descentdir, tolerance, maxiter, verbose, condtest)

    if nargin < 8
        condtest = false;
    end
           
    % Some algorithm parameters.
    eps       = 1e-8;   % A number close to zero.
    sigmamax  = 0.5;    % The maximum centering parameter.
    etamax    = 0.25;   % The maximum forcing number.
    mumin     = 1e-9;   % Minimum barrier parameter.
    alphamax  = 0.995;  % Maximum step size.
    alphamin  = 1e-4;   % Minimum step size.
    beta      = 0.75;   % Granularity of backtracking search.
    tau       = 0.01;   % Amount of actual decrease we will accept in 
                      % line search.

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

    if verbose
        fprintf('  i f(x)       lg(mu) sigma   ||rx||  ||rc||  alpha   #ls\n');
    end

    % Repeat while the convergence criterion has not been satisfied, and
    % we haven't reached the maximum number of iterations.
    alpha = zeros(1,ni);
    ls    = zeros(1,ni);
    saved = false(1,ni);
    
    for iter = 1:maxiter

        % COMPUTE OBJECTIVE, GRADIENT, CONSTRAINTS, ETC.  
        % Compute the response of the objective function, the gradient of the
        % objective, the response of the inequality constraints, the Jacobian of
        % the inequality constraints, the Hessian of the Lagrangian (minus the
        % Hessian of the objective) and, optionally, the Hessian of the
        % objective.
        [c J W] = constraint(x,z);
        if strcmp(descentdir,'newton')
            [f g B] = objective(x);
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
        eta        = min(etamax,auxr0);
        sigma      = min(sigmamax,sqrt(auxr0));
        dualitygap = -c.*z;
        mu         = max(mumin,sigma.*dualitygap/m);

        % Print the status of the algorithm.
        if verbose
            fprintf('%3d %+0.3e  %+5.2f %0.1e %0.1e %0.1e %0.1e %3d\n',...
                iter,f(1),log10(mu(1)),sigma(1),norm(rx(:,1)),norm(rc(1)),alpha(1),ls(1));
        end

        % CONVERGENCE CHECK.
        % If the norm of the responses is less than the specified tolerance,
        % and we have not saved the point before, save it!
        newToSave = (auxr0 < tolerance) & ~saved;
        % Mark as saved and save
        saved(newToSave) = true;
        finalX(:,newToSave) = x(:,newToSave);
        % If everyone is done, return
        if sum(~saved) == 0
            return
        end
        
%         % Update the BFGS approximation to the Hessian of the objective.
%         if strcmp(descentdir,'bfgs') && iter > 1
%             B = bfgsupdate(B,alpha*px,g-gprev);
%         end

        % SOLUTION TO PERTURBED KKT SYSTEM.
        % Compute the search direction of x and z.
        S  = z./(c-eps);
        gb = g - repmat(mu,n,1).*J.*(1./(repmat(c,n,1)-eps));
        px = zeros(n,ni);
        for iii = 1:ni,
            if ~saved(iii)
                myMatrix = B(:,:,iii) + W(:,:,iii) - J(:,iii)*S(iii)*J(:,iii)';
                if abs(det(myMatrix)) < eps || isnan(det(myMatrix)) || isinf(det(myMatrix))                   
                    finalX(:,iii) = x(:,iii);
                    saved(iii) = true;
                    continue
                else     
                    px(:,iii) = myMatrix \ (-gb(:,iii));
                end
            end
        end
        % If everyone is done, return
        if sum(~saved) == 0
            return
        end
        
        pz = -(z + mu./(c-eps) +  S.*sum(J.*px,1));

        % BACKTRACKING LINE SEARCH.
        % To ensure global convergence, execute backtracking line search to
        % determine the step length. First, we have to find the largest step
        % size which ensures that z remains feasible. Next, we perform
        % backtracking line search.
        alpha = alphamax*ones(1,ni);
        is    = find(z + pz < 0);
        if ~isempty(is)
            alpha(is) = alphamax * min(1, z(is) ./ -pz(is) );
        end

        % Compute the response of the merit function and the directional
        % gradient at the current point and search direction.
        psi  = merit(x,z,f,c,mu,eps);
        dpsi = gradmerit(x,z,px,pz,g,c,J,mu,eps);
        ls   = zeros(1,ni);
        backtrack = true(1,ni);
% xnew = zeros(size(x));
% znew = zeros(size(z));
% innerIt = zeros(1,ni);
% positions = 1:ni;
% figure
% hold on
        while true
%             fprintf('.');
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
%             ls(backtrack)     = ls(backtrack) + 1;
%             xnew(:,backtrack)   = x(:,backtrack) + repmat(alpha(backtrack),n,1) .* px(:,backtrack);
%             znew(backtrack)   = z(backtrack) + alpha(backtrack) .* pz(backtrack);
%             f(backtrack)      = objective(xnew(:,backtrack));
%             c(backtrack)      = constraint(xnew(:,backtrack),[]);
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
%             fprintf('%d %d | %d %d | %d %d | %d %d\n',size(psinew),size(eta),size(alpha),size(dpsi));
            stopBacktracking = ~(c > 0) & psinew < psi + tau*(eta.*alpha).*dpsi;
            stopBacktracking = stopBacktracking & backtrack;
            % For those which have to stop, save the point
            x(:,stopBacktracking)     = xnew(:,stopBacktracking);
            z(stopBacktracking)     = znew(stopBacktracking);
            
            % Candidate points which do not meet our criteria, will have a
            % smaller alpha next backtracking search
            alpha(~stopBacktracking) = alpha(~stopBacktracking) * beta;
            % Mark the candidates which were backtracked and that have
            % a too small alpha
            newToSave = alpha<alphamin & backtrack;

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
            
% plot(positions(saved),innerIt(saved),'go');
% plot(positions(backtrack),innerIt(backtrack),'rx');
% plot(positions(stopBacktracking),innerIt(stopBacktracking),'bs');
% plot(positions(newToSave),innerIt(newToSave),'dk');
% innerIt = innerIt+1;
% if innerIt(1) == 500
%     break
% end
            
        end
%         fprintf('\n');
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

% ------------------------------------------------------------------
% Update the quasi-Newton approximation using Broyden-Fletcher-
% Goldfarb-Shanno (BFGS) formula.
function B = bfgsupdate (B, s, y)  
    if y'*s < 0
        error('dot(y,s) > 0 is not satisfied');
    end
    x = B*s;
    B = B - x*x'/(x'*s) + y*y'/(y'*s);
end
