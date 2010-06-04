function state = ...
    psocheckbounds(state,Aineq,bineq,Aeq,beq,LB,UB,nonlcon,options)
% Check the the swarm population against constraints.

if isa(options.ConstrBoundary,'function_handle')
    boundcheckfcn = options.ConstrBoundary ;
elseif strcmpi(options.ConstrBoundary,'soft') || ...
        strcmpi(options.ConstrBoundary,'penalize')
    boundcheckfcn = @psoboundspenalize ;
elseif strcmpi(options.ConstrBoundary,'reflect')
    boundcheckfcn = @psoboundsreflect ;
elseif strcmpi(options.ConstrBoundary,'absorb')
    boundcheckfcn = @psoboundsabsorb ;
end

state = boundcheckfcn(state,Aineq,bineq,Aeq,beq,LB,UB,nonlcon,options) ;