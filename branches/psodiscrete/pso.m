function [xOpt,fval,exitflag,output,population,scores] = ...
    pso(fitnessfcn,nvars,Aineq,bineq,Aeq,beq,LB,UB,nonlcon,options)
% Find the minimum of a function using Particle Swarm Optimization.
% 
% This is an implementation of Particle Swarm Optimization algorithm using
% the same syntax as the Genetic Algorithm Toolbox, with some additional
% options specific to PSO. Allows code-reusability when trying different
% population-based optimization algorithms. Certain GA-specific parameters
% such as cross-over and mutation functions will not be applicable to the
% PSO algorithm. Demo function included, with a small library of test
% functions. Requires Optimization Toolbox.
% 
% In development, new features will be added regularly until this is made
% redundant by an official MATLAB PSO toolbox.
% 
% S. Chen. Version 20100522.
% Available from http://www.mathworks.com/matlabcentral/fileexchange/25986
% Distributed under BSD license.
% 
% Syntax:
% psodemo
% pso
% x = pso(fitnessfcn,nvars)
% x = pso(fitnessfcn,nvars,Aineq,bineq)
% x = pso(fitnessfcn,nvars,Aineq,bineq,Aeq,beq)
% x = pso(fitnessfcn,nvars,Aineq,bineq,Aeq,beq,LB,UB)
% x = pso(fitnessfcn,nvars,Aineq,bineq,Aeq,beq,LB,UB,nonlcon)
% x = pso(fitnessfcn,nvars,Aineq,bineq,Aeq,beq,LB,UB,nonlcon,options)
% x = pso(problem)
% [x, fval] = pso(...)
% [x, fval,exitflag] = pso(...)
% [x, fval,exitflag,output] = pso(...)
% [x, fval,exitflag,output,population] = pso(...)
% [x, fval,exitflag,output,population,scores] = pso(...)
% 
% Description:
% psodemo
% Runs a demonstration of the PSO algorithm using test function specified
% by user.
% 
% pso
% Runs a default demonstration using Rosenbrock's banana function.
% 
% x = pso(fitnessfcn,nvars)
% Runs the particle swarm algorithm with no constraints and default
% options. fitnessfcn is a function handle, nvars is the number of
% parameters to be optimized, i.e. the dimensionality of the problem.
% 
% x = pso(fitnessfcn,nvars,Aineq,bineq)
% Linear constraints, such that Aineq*x <= bineq. Soft boundaries only.
% Aineq is a matrix of size nconstraints x nvars, while b is a column
% vector of length nvars.
% 
% x = pso(fitnessfcn,nvars,Aineq,bineq,Aeq,beq)
% Linear equality constraints, such that Aeq*x = beq. Soft boundaries only.
% If no inequality constraints exist, set Aineq and bineq to [].
% 
% x = pso(fitnessfcn,nvars,Aineq,bineq,Aeq,beq,LB,UB)
% Lower and upper bounds definted as LB and UB respectively. Use empty
% arrays [] for A, b, Aeq, or beq if no linear constraints exist.
% 
% x = pso(fitnessfcn,nvars,Aineq,bineq,Aeq,beq,LB,UB,nonlcon)
% Non-linear constraints. Nonlinear inequality constraints in the form c(x)
% <= 0 have now been implemented using 'soft' or 'penalize' boundaries. See
% the Optimization Toolbox documentation for the proper syntax for defining
% nonlinear constraints.
% 
% x = pso(fitnessfcn,nvars,Aineq,bineq,Aeq,beq,LB,UB,nonlcon,options)
% Default algorithm parameters replaced with those defined in the options
% structure:
% Use >> options = psooptimset('Param1,'value1','Param2,'value2',...) to
% generate the options structure. Type >> psooptimset with no input or
% output arguments to display a list of available options and their
% default values.
% 
% NOTE: the swarm will perform better if the PopInitRange option is defined
% so that it closely bounds the expected domain of the feasible region.
% This is automatically done for lower and upper bound constraints, but not
% for linear and nonlinear constraints.
% 
% NOTE 2: If options.HybridFcn is to be defined, and if it is necessary to
% pass a non-default options structure to the hybrid function, then the
% options structure may be included in a cell array along with the pointer
% to the hybrid function. Example:
% >> % Let's say that we want to use fmincon to refine the result from PSO:
% >> hybridoptions = optimset(@fmincon) ;
% >> options.HybridFcn = {@fmincon, hybridoptions} ;
%
% NOTE 3:
% Perez and Behdinan (2007) demonstrated that the particle swarm is only
% stable if the following conditions are satisfied:
% Given that C0 = particle inertia
%            C1 = options.SocialAttraction
%            C2 = options.CognitiveAttraction
%    1) 0 < (C1 + C2) < 4
%    2) (C1 + C2)/2 - 1 < C0 < 1
% If conditions 1 and 2 are satisfied, the system will be guaranteed to
% converge to a stable equilibrium point. However, whether or not this
% point is actually the global minimum cannot be guaranteed, and its
% acceptability as a solution should be verified by the user.
% 
% x = pso(problem)
% Parameters imported from problem structure. Should work, but no error
% checking yet.
% 
% [x, fval] = pso(...)
% Returns the fitness value of the best solution.
% 
% [x, fval,exitflag] = pso(...)
% Returns information on outcome of pso run. Should match exitflag values
% for GA where applicable, for code reuseability.
% 
% [x, fval,exitflag,output] = pso(...)
% The structure output contains more detailed information about the PSO
% run. Should match output fields of GA, where applicable.
% 
% [x, fval,exitflag,output,population] = pso(...)
% A matrix population of size PopulationSize x nvars, with the locations of
% particles across the rows.
% 
% [x, fval,exitflag,output,population,scores] = pso(...)
% Final scores of the particles in population.
%
% Bibliography
% J Kennedy, RC Eberhart, YH Shi. Swarm Intelligence. Academic Press, 2001.
%
% SM Mikki, AA Kishk. Particle Swarm Optimization: A Physics-Based
% Approach. Morgan & Claypool, 2008.
%
% RE Perez and K Behdinan. "Particle swarm approach for structural
% design optimization." Computers and Structures, Vol. 85:1579-88, 2007.
% 
% See also:
% PSODEMO, PSOOPTIMSET, PSOBINARY.

if ~nargin % Rosenbrock's banana function by default, as a demonstration
    fitnessfcn = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2 ;
    nvars = 2 ;
    LB = [-1.5,-2] ;
    UB = [2,2] ;
    options.PopInitRange = [[-2;4],[-1;2]] ;
    options.PlotFcns = {@psoplotbestf,@psoplotswarmsurf} ;
    options.Generations = 200 ;
    options.DemoMode = 'on' ;
    options.KnownMin = [1 1] ;
elseif isstruct(fitnessfcn)
    nvars = fitnessfcn.nvars ;
    Aineq = fitnessfcn.Aineq ;
    bineq = fitnessfcn.bineq ;
    Aeq = fitnessfcn.Aeq ;
    beq = fitnessfcn.beq ;
    LB = fitnessfcn.LB ;
    UB = fitnessfcn.UB ;
    nonlcon = fitnessfcn.nonlcon ;
    if ischar(nonlcon) && ~isempty(nonlcon)
        nonlcon = str2func(nonlcon) ;
    end
    options = fitnessfcn.options ;
    fitnessfcn = fitnessfcn.fitnessfcn ;
elseif nargin < 2
    msg = 'PSO requires at least two input arguments' ;
    error('%s, or a problem structure. Type >> help pso for details',...
        msg)
end % if ~nargin

if ~exist('options','var') % Set default options
    options = struct ;
end % if ~exist

options = psooptimset(options) ;

options.Verbosity = 1 ; % For options.Display == 'final' (default)
if strcmpi(options.Display,'off')
    options.Verbosity = 0 ;
elseif strncmpi(options.Display,'iter',4)
    options.Verbosity = 2 ;
elseif strncmpi(options.Display,'diag',4)
    options.Verbosity = 3 ;
end

if ~exist('Aineq','var'), Aineq = [] ; end
if ~exist('bineq','var'), bineq = [] ; end
if ~exist('Aeq','var'), Aeq = [] ; end
if ~exist('beq','var'), beq = [] ; end
if ~exist('LB','var'), LB = [] ; end
if ~exist('UB','var'), UB = [] ; end
if ~exist('nonlcon','var'), nonlcon = [] ; end

% Check for swarm stability
% -------------------------------------------------------------------------
if options.SocialAttraction + options.CognitiveAttraction >= 4
    msg = 'Warning: Swarm is unstable and may not converge ' ;
    msg = [msg 'since the sum of the cognitive and social attraction'] ;
    msg = [msg ' parameters is greater than or equal to 4.'] ;
    msg = [msg ' Suggest adjusting options.CognitiveAttraction and/or'] ;
    sprintf('%s options.SocialAttraction.',msg)
end
% -------------------------------------------------------------------------

% Check for constraints and bit string population type
% -------------------------------------------------------------------------
if strncmpi(options.PopulationType,'bitstring',2) && ...
        (~isempty([Aineq,bineq]) || ~isempty([Aeq,beq]) || ...
        ~isempty(nonlcon) || ~isempty([LB,UB]))
    Aineq = [] ; bineq = [] ; Aeq = [] ; beq = [] ; nonlcon = [] ;
    LB = [] ; UB = [] ;
    msg = sprintf('Warning: Constraints will be ignored') ;
    msg = sprintf('%s for options.PopulationType ''bitstring''',msg) ;
    disp(msg)
end
% -------------------------------------------------------------------------

% Change this when nonlcon gets fully implemented:
% -------------------------------------------------------------------------
if ~isempty(nonlcon) && strcmpi(options.ConstrBoundary,'reflect')
    msg = 'Non-linear constraints don''t have ''reflect'' boundaries' ;
    msg = [msg, ' implemented yet.'] ;
    warning('pso:main:nonlcon',...
        '%s Changing options.ConstrBoundary to ''penalize''.',...
        msg)
    options.ConstrBoundary = 'penalize' ;
end
% -------------------------------------------------------------------------

% Is options.PopInitRange reconcilable with LB and UB constraints?
% -------------------------------------------------------------------------
% Resize PopInitRange in case it was given as one range for all dimensions
if size(options.PopInitRange,2) == 1 && nvars > 1
    options.PopInitRange = repmat(options.PopInitRange,1,nvars) ;
end

% Check initial population with respect to bound constraints
% Is this really desirable? Maybe there are some situations where the user
% specifically does not want an uniform inital population covering all of
% LB and UB?
if ~isempty(LB) || ~isempty(UB)
    options.LinearConstr.type = 'boundconstraints' ;
    options.PopInitRange = ...
        psocheckpopulationinitrange(options.PopInitRange,LB,UB) ;
end
% -------------------------------------------------------------------------

% Check validity of VelocityLimit
% -------------------------------------------------------------------------
if all(~isfinite(options.VelocityLimit))
    options.VelocityLimit = [] ;
elseif isscalar(options.VelocityLimit)
    options.VelocityLimit = repmat(options.VelocityLimit,1,nvars) ;
elseif ~isempty(length(options.VelocityLimit)) && ...
        ~isequal(length(options.VelocityLimit),nvars)
    msg = 'options.VelocityLimit must be either a positive scalar' ;
    error('%s, or a vector of size 1xnvars.',msg)
end % if isscalar
options.VelocityLimit = abs(options.VelocityLimit) ;

% Generate swarm initial state (this line must not be moved)
if strncmpi(options.PopulationType,'double',2)
    state = psocreationuniform(options,nvars) ;
elseif strncmpi(options.PopulationType,'bi',2)
    state = psocreationbinary(options,nvars) ;
end
% -------------------------------------------------------------------------

% Check initial population with respect to linear and nonlinear constraints
% -------------------------------------------------------------------------
if ~isempty(Aeq) || ~isempty(Aineq) || ~isempty(nonlcon)
    options.LinearConstr.type = 'linearconstraints' ;
    if ~isempty(nonlcon)
        options.LinearConstr.type = 'nonlinearconstraints' ;
    end
    if strcmpi(options.ConstrBoundary,'reflect')
        options.ConstrBoundary = 'soft' ;
        msg = sprintf('Constraint boundary behavior ''reflect''') ;
        msg = sprintf('%s is not yet supported for linear constraints.',...
            msg) ;
        msg = sprintf('%s Switching boundary behavior to ''soft''.',msg) ;
        warning('pso:mainfcn:constraintbounds',...
            '%s',msg)
    end
    [state,options] = psocheckinitialpopulation(state,...
        Aineq,bineq,Aeq,beq,...
        LB,UB,...
        nonlcon,...
        options) ;
end
% -------------------------------------------------------------------------

% Check constraint type
% -------------------------------------------------------------------------
if isa(options.ConstrBoundary,'function_handle')
    boundcheckfcn = options.ConstrBoundary ;
elseif strcmpi(options.ConstrBoundary,'soft')
    boundcheckfcn = @psoboundssoft ;
elseif strcmpi(options.ConstrBoundary,'penalize')
    boundcheckfcn = @psoboundspenalize ;
    state.Penalty = zeros(options.PopulationSize,1) ;
elseif strcmpi(options.ConstrBoundary,'reflect')
    boundcheckfcn = @psoboundsreflect ;
elseif strcmpi(options.ConstrBoundary,'absorb')
    boundcheckfcn = @psoboundsabsorb ;
end
% -------------------------------------------------------------------------

n = options.PopulationSize ;

if ~isempty(options.PlotFcns)
    close(findobj('Tag','Swarm Plots','Type','figure'))
    state.hfigure = figure('NumberTitle','off',...
        'Name','PSO Progress',...
        'Tag','Swarm Plots') ;
end % if ~isempty

if options.Verbosity > 0, fprintf('\nSwarming...'), end
exitflag = 0 ; % Default exitflag, for max iterations reached.
flag = 'init' ;

state.fitnessfcn = fitnessfcn ;
state.LastImprovement = 1 ;
state.ParticleInertia = 0.9 ; % Initial inertia
state.OutOfBounds = false(options.PopulationSize,1) ;
% alpha = 0 ;


% Iterate swarm
for k = 1:options.Generations
    state.Score = inf*ones(n,1) ; % Reset fitness vector
    state.Generation = k ;
    state.OutOfBounds = false(options.PopulationSize,1) ;
    
    % Check bounds before proceeding
    % ---------------------------------------------------------------------
    if ~isequal(boundcheckfcn,@psoboundspenalize) && ...
            ~all([isempty([Aineq,bineq]), isempty([Aeq,beq]), ...
            isempty([LB;UB]), isempty(nonlcon)])
        state = boundcheckfcn(state,Aineq,bineq,Aeq,beq,LB,UB,nonlcon,...
            options) ;
    end % if ~isempty
    % ---------------------------------------------------------------------
    
    % Evaluate fitness
    % ---------------------------------------------------------------------
    if strcmpi(options.Vectorized,'off')
        for i = setdiff(1:n,find(state.OutOfBounds))
            state.Score(i) = fitnessfcn(state.Population(i,:)) ;
        end % for i
    else % Vectorized fitness function
        state.Score(setdiff(1:n,find(state.OutOfBounds))) = ...
            fitnessfcn(state.Population(setdiff(1:n,...
            find(state.OutOfBounds)),:)) ;
    end % if strcmpi
    % ---------------------------------------------------------------------
    
    % Check bounds before proceeding ('penalize' needs special treatment)
    % ---------------------------------------------------------------------
    if isequal(boundcheckfcn,@psoboundspenalize) && ...
            ~all([isempty([Aineq,bineq]), isempty([Aeq,beq]), ...
            isempty([LB;UB]), isempty(nonlcon)])
        state = boundcheckfcn(state,Aineq,bineq,Aeq,beq,LB,UB,nonlcon,...
            options) ;
    end % if ~isempty
    % ---------------------------------------------------------------------
    
    % Update the local bests
    % ---------------------------------------------------------------------
    betterindex = state.Score < state.fLocalBests ;
    state.fLocalBests(betterindex) = state.Score(betterindex) ;
    state.xLocalBests(betterindex,:) = ...
        state.Population(betterindex,:) ;
    % ---------------------------------------------------------------------
    
    % Update the global best and its fitness, then check for termination
    % ---------------------------------------------------------------------
    [minfitness, minfitnessindex] = min(state.Score) ;
    
%     alpha = alpha + (1/k) * ...
%         ((1/n)*sum((state.Velocities*state.Velocities')^2) ./ ...
%         ((1/n)*sum(state.Velocities*state.Velocities')).^2) ;
%     tempchk = alpha <= 1.6 ;
    if minfitness < state.fGlobalBest
        state.fGlobalBest(k) = minfitness ;
        state.xGlobalBest = state.Population(minfitnessindex,:) ;
        state.LastImprovement = k ;
        imprvchk = k > options.StallGenLimit && ...
            (state.fGlobalBest(k - options.StallGenLimit) - ...
                state.fGlobalBest(k)) / (k - options.StallGenLimit) < ...
                options.TolFun ;
        if imprvchk
            exitflag = 1 ;
            flag = 'done' ;
        elseif state.fGlobalBest(k) < options.FitnessLimit
            exitflag = 2 ;
            flag = 'done' ;
        end % if k
    elseif k > 1 % No improvement from last iteration
        state.fGlobalBest(k) = state.fGlobalBest(k-1) ;
    end % if minfitness
    
    stallchk = k - state.LastImprovement >= options.StallGenLimit ;
    if stallchk
        % No improvement for StallGenLimit generations
        exitflag = 3 ;
        flag = 'done' ;
    end
    % ---------------------------------------------------------------------
    
    % Update flags, state and plots before updating positions
    % ---------------------------------------------------------------------
    if k == 2, flag = 'iter' ; end
    if k == options.Generations
        flag = 'done' ;
        exitflag = 0 ;
    end
    
    if ~isempty(options.PlotFcns) && ~mod(k,options.PlotInterval)
        % Exit gracefully if user has closed the figure
        if isempty(findobj('Tag','Swarm Plots','Type','figure'))
            exitflag = -1 ;
            break
        end % if isempty
        % Find a good size for subplot array
        rows = floor(sqrt(length(options.PlotFcns))) ;
        cols = ceil(length(options.PlotFcns) / rows) ;
        % Cycle through plotting functions
        if strcmpi(flag,'init')
            haxes = zeros(length(options.PlotFcns),1) ;
        end % if strcmpi
        for i = 1:length(options.PlotFcns)
            if strcmpi(flag,'init')
                haxes(i) = subplot(rows,cols,i,...
                    'Parent',state.hfigure) ;
                set(gca,'NextPlot','replacechildren')
            else
                subplot(haxes(i))
            end % if strcmpi
            state = options.PlotFcns{i}(options,state,flag) ;
        end % for i
        drawnow
    end % if ~isempty
    
    if ~isempty(options.OutputFcns) && ~mod(k,options.PlotInterval)
        for i = 1:length(options.Output)
            state = options.OutputFcns{i}(options,state,flag) ;
        end % for i
    end % if ~isempty
    
    if strcmpi(flag,'done')
        break
    end % if strcmpi
    % ---------------------------------------------------------------------
    
    % Update the particle velocities and positions
    state = psoiterate(state,options) ;
end % for k

% Assign output variables and generate output
% -------------------------------------------------------------------------
xOpt = state.xGlobalBest ;
fval = state.fGlobalBest(k) ; % Best fitness value
% Final population: (hopefully very close to each other)
population = state.Population ;
scores = state.Score ; % Final scores (NOT local bests)
output.generations = k ; % Number of iterations performed
clear state

output.message = psogenerateoutputmessage(options,output,exitflag) ;
if options.Verbosity > 0, fprintf('\n\n%s\n',output.message) ; end
% -------------------------------------------------------------------------

% Check for hybrid function, run if necessary
% -------------------------------------------------------------------------
if ~isempty(options.HybridFcn) && exitflag ~= -1
    [xOpt,fval] = psorunhybridfcn(fitnessfcn,xOpt,Aineq,bineq,...
        Aeq,beq,LB,UB,nonlcon,options) ;
end
% -------------------------------------------------------------------------

% Wrap up
% -------------------------------------------------------------------------
if options.Verbosity > 0
    if exitflag == -1
        fprintf('\nBest point found: %s\n\n',mat2str(xOpt,5))
    else
        fprintf('\nFinal best point: %s\n\n',mat2str(xOpt,5))
    end
end % if options.Verbosity

if ~nargout, clear all, end
% -------------------------------------------------------------------------