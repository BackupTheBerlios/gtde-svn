function solution = branchAndBound(problem)

    % Check problem structure, mandatory and default optional fields.
    problem = checkOptions(problem);
    
    % Cost function and lipschitz constant
    f = problem.costFunction;
    L = problem.constantLipschitz;
    
    % Define the entry point to the problem
    P0.x = problem.startingPoint;
    P0.s = problem.startingSize;
    P0.l = lowerBound(f,P0,L);
    
    % Define the algorithm variables
    % Lists
    List = {P0};
    NList = cell(0);
    % Level
    lv = 1;
    % Current size
    s = P0.s;
    % Number of branches
    nb = 0;
    % Upper bound
    u = problem.upperBound;
    
    % Main loop, bound and branch
    while s > problem.tolerance && nb < problem.maxBranches,
        % Printf
        fprintf('Lv = %d;  |List| = %d; u = %1.2f .\n',lv,numel(List),u);
        % Process all the elements in list
        for pp=1:numel(List)
            P = List{pp};
            % Check bounds
            if P.l <= u,
                % Update problem upper bound
                u = upperBound(f,P,L);
                % Increase nb and branch
                nb = nb + 1;
                auxList = branch(P,problem,lv);
                % Save to NList
                NList = cat(1,NList,auxList);
            end
        end
        % The next level's list becomes the current list and we empty new
        % level's list
        List = NList;
        NList = cell(0);
        % Update s and lv
        s = List{1}.s;
        lv = lv + 1;
    end
    
    % Retrieve optimal solution
    for ii = 1:numel(List),
        if f(List{ii}.x) < u
            u = f(List{ii}.x);
            solution = List{ii}.x;
        end
    end
    for ii = 1:numel(NList),
        if f(NList{ii}.x) < u
            u = f(NList{ii}.x);
            solution = NList{ii}.x;
        end
    end
    

end

% Branch
function list = branch(P,problem,lv)
    
    % Define
    list = cell(0);
    % Variables
    dim = numel(P.x);
    s = P.s/2;

    % Constrained problem
    for d = 1:dim,
        v = zeros(size(P.x));
        v(d) = 1;
        % Add first point
        Q.x = P.x + s*v;
        Q.s = s;
        Q.l = lowerBound(problem.costFunction,Q,problem.constantLipschitz);
        if feasibleRegion(Q,problem),
            list = cat(1,list,Q);
        end
        % Add second point
        Q.x = P.x - s*v;
        Q.s = s;
        Q.l = lowerBound(problem.costFunction,Q,problem.constantLipschitz);
        list = cat(1,list,Q);
        if feasibleRegion(Q,problem),
            list = cat(1,list,Q);
        end
    end        

end

% Feasible region
function isIt = feasibleRegion(Q,problem)

    % Default is yes
    isIt = true;
    % Shall we check the constraint?
    if isfield(problem,'constrain') && lv > problem.levelConvex
        
    end

end

% Upper bound function
function u = upperBound(f,P,L)
    u = f(P.x) + L*P.s*sqrt(2);
end

% Lower bound function
function u = lowerBound(f,P,L)
    u = f(P.x) - L*P.s*sqrt(2);
end

% Check problem structure, mandatory and default optional fields.
function p = checkOptions(p)

    % Check the mandatory fields of problem
    mandatoryFields = {...
        'costFunction','constantLipschitz',...
        'startingPoint','startingSize'};
    for f = 1:numel(mandatoryFields),
        if ~isfield(p,mandatoryFields{f})
            error(['Field ' mandatoryFields{f} ' is mandatory.']);
        end
    end
    
    % Default options
    % Field names
    defaultValues = {...
        {'upperBound',      inf    },...
        {'tolerance',       p.startingSize/10  },...
        {'maxBranches',         200          },...
        {'levelConvex',         3       },...
        };
    
    % Loop on it
    for no = 1:numel(defaultValues),
        if ~isfield(p,defaultValues{no}{1})
            p.(defaultValues{no}{1}) = defaultValues{no}{2};
        end
    end
    
end