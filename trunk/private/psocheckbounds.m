function state = ...
    psocheckbounds(options,state,Aineq,bineq,Aeq,beq,LB,UB,nonlcon)
% Check the the swarm population against boundary and linear constraints.

x = state.Population ;
v = state.Velocities ;

state.OutOfBounds = zeros(size(state.Population,1),1) ;

for i = 1:size(state.Population,1)
    lowindex = [] ; highindex = [] ;
    if ~isempty(LB), lowindex = x(i,:) < LB ; end
    if ~isempty(UB), highindex = x(i,:) > UB ; end
    % Three constraint types
    if strcmpi(options.ConstrBoundary,'soft')
        outofbounds = any([lowindex,highindex]) ;
        if ~outofbounds && ~isempty(Aineq)
            outofbounds = any(Aineq*x(i,:)' - bineq > options.TolCon) ;
        end % if ~isempty
        if ~outofbounds && ~isempty(Aeq)
            outofbounds = any(abs(Aeq*x(i,:)' - beq) > options.TolCon) ;
        end % if ~isempty
        
        if ~outofbounds && ~isempty(nonlcon) % Nonlinear constraint check
            c = nonlcon(x(i,:)) ;
            outofbounds = any(c > options.TolCon) ;
%                 any(abs(ceq) > options.TolCon);
        end
        
        if outofbounds
            state.Score(i) = inf ;
        end % if outofbounds
        
        state.OutOfBounds(i) = outofbounds ;
    elseif strcmpi(options.ConstrBoundary,'reflect')
        x(i,lowindex) = LB(lowindex) ;
        x(i,highindex) = UB(highindex) ;
        v(i,lowindex) = -v(i,lowindex) ;
        v(i,highindex) = -v(i,highindex);
    elseif strcmpi(options.ConstrBoundary,'absorb')
        % Check against bounded constraints
        x(i,lowindex) = LB(lowindex) ;
        x(i,highindex) = UB(highindex) ;
        v(i,lowindex) = 0 ;
        v(i,highindex) = 0 ;
        
        % Linear and nonlinear constraints
        if ~isempty(Aineq) || ~isempty(Aeq) || ~isempty(nonlcon)
            % "Sticky" linear inequality constraints
            if ~isempty(Aineq)
                if max(Aineq*x(i,:)' - bineq) > options.TolCon
                    v(i,:) = 0 ;
                end % if Aineq
            end % if ~isempty
            
            % Won't do set velocities to zero for particles outside of
            % equality constraints, or else particles will rarely ever
            % move. This could change if "slippery" bounds are implemented
            % for linear constraints.
            
            % Finally update all particle positions
            if isempty(nonlcon)
                x(i,:) = linprog([],Aineq,bineq,Aeq,beq,LB,UB,...
                    x(i,:),state.LinprogOptions) ;
            else % Check nonlinear constraints
                [c,ceq] = nonlcon(state.Population(i,:)) ;
                if any(c > options.TolCon) || ...
                        any(abs(ceq) > options.TolCon)
                    v(i,:) = 0 ; % Sticky
                    x(i,:) = fmincon(@void,state.Population(i,:),...
                        Aineq,bineq,Aeq,beq,LB,UB,...
                        nonlcon,state.LinprogOptions) ;
                end % if any
            end % if isempty
        end % if ~isempty
    end % if strcmpi
end % for i

state.Population = x ;
state.Velocities = v ;