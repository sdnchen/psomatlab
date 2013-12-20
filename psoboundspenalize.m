function state = ...
    psoboundspenalize(state,Aineq,bineq,Aeq,beq,LB,UB,nonlcon,options)
% Penalty-based constraint enforcement method.

% Unpack variables from state structure (necessary for parfor)
x = state.Population ;
v = state.Velocities ;
n = size(state.Population,1) ;
nvars = size(x,1) ;
OutOfBounds = false(n,1) ;
tol = options.TolCon ;
% state.InBounds = false(n,1) ;

% Check for existence of constraints. Don't bother running rest of the
% constraint-checking code if problem is unconstrained.
% This could be done much more efficiently, by deciding early on if the
% problem is actually constrained or not. If there are actually no
% constraints defined, then don't even bother going into the constraint
% handling function.
[c,ceq] = nonlconwrapper(nonlcon,Aineq,bineq,Aeq,beq,LB,UB,...
    tol,x(1,:)) ;
if isempty([c,ceq])
    state.ConstrViolations = zeros(n,1) ;
    return
end

nbrConstraints = size([c ceq],2) ;
ConstrViolations = zeros(nvars,nbrConstraints) ;
if ~strcmpi(options.UseParallel,'always') % Default option
    for i = 1:n
        [c,ceq] = nonlconwrapper(nonlcon,Aineq,bineq,Aeq,beq,LB,UB,...
            tol,x(i,:)) ;

        % Tolerances already dealt with in nonlconwrapper
        if sum([c,ceq]) ~= 0
            OutOfBounds(i) = true ;
             % Sticky boundaries: kills the inertia of the particle if
             % it is in a non-feasible region of the design space.
            v(i,:) = 0 ;
        end % if sum
        ConstrViolations(i,:) = [c,ceq] ;
    end % for i
else % Parallel computing option
    parfor i = 1:n
        [c,ceq] = nonlconwrapper(nonlcon,Aineq,bineq,Aeq,beq,LB,UB,...
            tol,x(i,:)) ;
        % Tolerances already dealt with in nonlconwrapper
        if sum([c,ceq]) ~= 0
            OutOfBounds(i) = true ;
            % Sticky boundaries: kills the inertia of the particle if
            % it is in a non-feasible region of the design space.
            v(i,:) = 0 ;
        end % if sum
        ConstrViolations(i,:) = [c,ceq] ;
    end % parfor i
end % if strcmpi

% state.InBounds(setdiff((1:n)',find(state.OutOfBounds))) = true ;
% Repack variables into state structure
state.Population = x ;
state.Velocities = v ;
state.OutOfBounds = OutOfBounds ;
state.ConstrViolations = ConstrViolations ;

function [c,ceq] = nonlconwrapper(nonlcon,Aineq,bineq,Aeq,beq,LB,UB,...
    TolCon,x)
% Wrapper function for combining evaluation of all constraints

c = [] ; ceq = [] ;
if ~isempty(nonlcon)
    [c,ceq] = nonlcon(x) ;
    c = reshape(c,1,[]) ; ceq = reshape(ceq,1,[]) ; % Robustness
end

if ~isempty(Aineq)
    c = [c, (Aineq*x' - bineq)'] ;
end

if ~isempty(LB)
    c = [c, LB - x] ;
end

if ~isempty(UB)
    c = [c, x - UB] ;
end

if ~isempty(Aeq)
    ceq = [ceq, abs(Aeq*x' - beq)'] ;
end

% Tolerances
if ~isempty(c), c(c < 0) = 0 ; end
if ~isempty(ceq), ceq(abs(ceq) < TolCon) = 0 ; end