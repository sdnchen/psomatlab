function state = ...
    psoboundspenalize(state,Aineq,bineq,Aeq,beq,LB,UB,nonlcon,options)
% This is like "soft" boundaries, except that some kind of penalty value
% must be calculated from the degree of each constraint violation.

x = state.Population ;
% v = state.Velocities ;

state.OutOfBounds = false(size(state.Population,1),1) ;
state.Penalty = zeros(size(state.Population,1),1) ;

nLB = size(LB,2) ;
nUB = nLB + size(UB,2) ;
nineq = nUB + size(Aineq,1) ;
neq = nineq + size(Aeq,1) ;

if ~isempty(nonlcon)
    [ctest,ceqtest] = nonlcon(zeros(1,options.PopulationSize)) ;
    ctest = ctest(:) ; ceqtest = ceqtest(:) ;
    nnonl = neq + size([ctest;ceqtest],1) ;
else
    nnonl = neq ;
end

nconstr = nnonl ;
f = abs(mean(state.Score)) ;
g = zeros(options.PopulationSize,nconstr) ;

for i = 1:options.PopulationSize
    if ~isempty(LB)
        g(i,1:nLB) = max([LB - x(i,:) ;
            zeros(1,size(state.Population,2))]) ;
    end
    
    if ~isempty(UB)
        g(i,nLB+1:nUB) = max([x(i,:) - UB ;
            zeros(1,size(state.Population,2))]) ;
    end
    
    if ~isempty(Aineq) % Check linear inequalities
        g(i,nUB+1:nineq) = max(Aineq*x(i,:)' - bineq,...
            zeros(size(bineq))) ;
    end % if ~isempty
    
    if ~isempty(Aeq) % Check linear equalities
        g(i,nineq+1:neq) = max(abs(Aeq*x(i,:)' - beq),...
            zeros(size(beq))) ;
    end % if ~isempty
    
    if ~isempty(nonlcon) % Nonlinear constraint check
        [c,ceq] = nonlcon(x(i,:)) ; 
        g(i,neq+1:nnonl) = ...
            [max(c(:)' - options.TolCon,zeros(size(c(:)'))) , ...
             max(abs(ceq(:)') > options.TolCon,zeros(size(ceq(:)')))] ;
    end
    
    if any(g(i,:)), state.Velocities(i,:) = 0 ; end
end % for i

state.Penalty = calculatepenalties(f,g,options) ;
state.Score = state.Score + state.Penalty ;
state.Population = x ;
% state.Velocities = v ;

function penalty = calculatepenalties(f,g,options)

ssg = sum(mean(g,1).^2,2) ;
penalty = zeros(options.PopulationSize,1) ;
if ssg > options.TolCon ;
    k = zeros(1,size(g,2)) ;
    for i = 1:size(g,1)
        for j = 1:size(g,2)
            k(i,j) = f*mean(g(:,j),1)/ssg ;
        end
        penalty(i) = k(i,:)*g(i,:)' ;
    end
end