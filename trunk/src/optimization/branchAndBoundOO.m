function solution = branchAndBoundOO(problem)

    % Check problem structure, mandatory and default optional fields.
    problem = checkOptions(problem);
    
    % Cost function and lipschitz constant
    f = problem.costFunction;
    L = problem.constantLipschitz;
    
    % Define the entry point to the problem
    X0 = problem.startingPoint;
    
    % Define the algorithm variables
    % Lists
    List = X0;
    NList = [];
    % Level
    lv = 1;
    % Current size
    s = problem.startingSize;
    if numel(s) == 1,
        s = s*ones(size(X0));
    end
    % Number of branches
    nb = 0;
    % Upper bound
    u = problem.upperBound;
    
    % Main loop, bound and branch
    while max(s) > problem.tolerance && nb < problem.maxBranches,
        % Compute bounds
        lB = lowerBound(f,List,s,L);
        uB = upperBound(f,List,s,L);
        % Remove those who do not satisfy the constraint
        u = min(uB);
        
%         figure;
%         plot3(List(1,lB > u),List(2,lB > u),List(3,lB > u),'rx')
%         hold on;
%         plot3(List(1,lB <= u),List(2,lB <= u),List(3,lB <= u),'go')
        
        List(:,lB > u) = [];
        
        % Printf
%         fprintf('Lv = %d;  |List| = %d; u = %1.2g .\n',lv,size(List,2),u);
        
        % Increase nb and branch
        nb = nb + size(List,2);
        
        % Process all the elements in list
        List = branch(List,problem,s);

        % Update s and lv
        s = s/2;
        lv = lv + 1;
    end
    
%     % Retrieve optimal solution
%     fval = f(List);
%     
%     % Impose the constraint
%     if isfield(problem,'constraint')
%         ind = problem.constraint(List) <= 0;
%         % Remove those that do not satisfy the constraint
%         fval = fval(ind);
%         List = List(:,ind);
%     end
% 
%     [~,ind] = min(fval);
%     solution = List(:,ind);

    % Give solution
    solution = List;
    

end

% Branch
function olist = branch(ilist,problem,s)
    
    % Define
    olist = [];
    
    % Variables
    dim = size(ilist,1);
    s = s/2;
    
    % MyBasis
    B = (2*dec2bin(0:2^dim-1,dim)-97)';

    % Constrained problem
    for bb = 1:size(B,2),
        v = s.*B(:,bb);
        % Add first point
        x = ilist + repmat(v,1,size(ilist,2));
        % Cat
        olist = cat(2,olist,x);
    end        

end

% Feasible region
function isIt = feasibleRegion(Q,problem)

    % Default is yes
    isIt = true;
    % Shall we check the constraint?
    if isfield(problem,'constraint') && lv > problem.levelConvex
        
    end

end

% Upper bound function
function u = upperBound(f,x,s,L)
    u = f(x) + L*norm(s)*sqrt(2);
end

% Lower bound function
function u = lowerBound(f,x,s,L)
    u = f(x) - L*norm(s)*sqrt(2);
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
        {'tolerance',       min(p.startingSize)/10  },...
        {'maxBranches',         1000          },...
        {'levelConvex',         3       },...
        };
    
    % Loop on it
    for no = 1:numel(defaultValues),
        if ~isfield(p,defaultValues{no}{1})
            p.(defaultValues{no}{1}) = defaultValues{no}{2};
        end
    end
    
end