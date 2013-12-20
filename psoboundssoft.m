function state = ...
    psoboundssoft(state,Aineq,bineq,Aeq,beq,LB,UB,nonlcon,options)
% Internal toolbox function.
%
% Decide whether a particle violates any constraints, and assign it a score
% of +infinity if it does.

x = state.Population ;
% v = state.Velocities ;
n = size(x,1) ;
OutOfBounds = false(n,1) ;
tol = options.TolCon ;
score = state.Score ;

if ~strcmpi(options.UseParallel,'always')
    for i = 1:n
        lowindex = [] ; highindex = [] ;
        if ~isempty(LB), lowindex = x(i,:) < LB ; end
        if ~isempty(UB), highindex = x(i,:) > UB ; end

        outofbounds = any([lowindex,highindex]) ;
        if ~outofbounds && ~isempty(Aineq) % Check linear inequalities
            outofbounds = any(Aineq*x(i,:)' - bineq > tol) ;
        end % if ~isempty
        if ~outofbounds && ~isempty(Aeq) % Check linear equalities
            outofbounds = any(abs(Aeq*x(i,:)' - beq) > tol) ;
        end % if ~isempty
        if ~outofbounds && ~isempty(nonlcon) % Nonlinear constraint check
            [c,ceq] = nonlcon(x(i,:)) ;
            outofbounds = any(c > tol) ;
            outofbounds = outofbounds || any(abs(ceq) > tol) ;
        end

        if outofbounds
            score(i) = inf ;
        end % if outofbounds

        OutOfBounds(i) = outofbounds ;
    end % for i
else % Parallel computing option
    parfor i = 1:n
        lowindex = [] ; highindex = [] ;
        if ~isempty(LB), lowindex = x(i,:) < LB ; end
        if ~isempty(UB), highindex = x(i,:) > UB ; end

        outofbounds = any([lowindex,highindex]) ;
        if ~outofbounds && ~isempty(Aineq) % Check linear inequalities
            outofbounds = any(Aineq*x(i,:)' - bineq > tol) ;
        end % if ~isempty
        if ~outofbounds && ~isempty(Aeq) % Check linear equalities
            outofbounds = any(abs(Aeq*x(i,:)' - beq) > tol) ;
        end % if ~isempty
        if ~outofbounds && ~isempty(nonlcon) % Nonlinear constraint check
            [c,ceq] = nonlcon(x(i,:)) ;
            outofbounds = any(c > tol) ;
            outofbounds = outofbounds || any(abs(ceq) > tol) ;
        end

        if outofbounds
            score(i) = inf ;
        end % if outofbounds

        OutOfBounds(i) = outofbounds ;
    end % for i
end

% Repack variables into state structure
state.Population = x ;
% state.Velocities = v ;
state.OutOfBounds = OutOfBounds ;
state.Score = score ;
