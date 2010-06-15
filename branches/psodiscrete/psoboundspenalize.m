function state = ...
    psoboundspenalize(state,Aineq,bineq,Aeq,beq,LB,UB,nonlcon,options)
% This is like "soft" boundaries, except that some kind of penalty value
% must be calculated from the degree of each constraint violation.

x = state.Population ;
% v = state.Velocities ;

state.OutOfBounds = false(size(state.Population,1),1) ;
state.Penalty = zeros(size(state.Population,1),1) ;

if ~isempty(nonlcon)
    [ctest,ceqtest] = nonlcon(zeros(1,options.PopulationSize)) ;
    ctest = ctest(:) ; ceqtest = ceqtest(:) ;
    nnonl = size([ctest;ceqtest],1) ;
else
    nnonl = 0 ;
end

nLB = size(LB,2) ;
nUB = size(UB,2) ;
nineq = size(Aineq,1) ; neq = size(Aeq,1) ;

nconstr = nLB + nUB + nineq + neq + nnonl ;
f = abs(mean(state.Score)) ;
g = zeros(options.PopulationSize,nconstr) ;

for i = 1:options.PopulationSize
    if ~isempty(LB)
        g(i,1:nLB) = max([LB - x(i,:) - options.TolCon ;
            zeros(1,size(state.Population,2))]) ;
    end
    
    if ~isempty(UB)
        g(i,nLB+1:nLB+nUB) = max([x(i,:) - UB  - options.TolCon ;
            zeros(1,size(state.Population,2))]) ;
    end
    
    if ~isempty(Aineq) % Check linear inequalities
        g(i,nUB+1:nUB+nineq) = max(Aineq*x(i,:)' - bineq - options.TolCon,...
            zeros(size(bineq))) ;
    end % if ~isempty
    
    if ~isempty(Aeq) % Check linear equalities
        g(i,nineq+1:nineq+neq) = max(abs(Aeq*x(i,:)' - beq) - ...
            options.TolCon,zeros(size(beq))) ;
    end % if ~isempty
    
    if ~isempty(nonlcon) % Nonlinear constraint check
        [c,ceq] = nonlcon(x(i,:)) ; 
        g(i,neq+1:neq+nnonl) = ...
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
    penalty = zeros(size(g,1),1) ;
    for i = 1:size(g,1)
        for j = 1:size(g,2)
            k(i,j) = f*mean(g(:,j),1)/ssg ;
        end
        penalty(i) = k(i,:)*g(i,:)' ;
    end
end