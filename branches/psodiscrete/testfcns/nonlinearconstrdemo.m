function f = nonlinearconstrdemo(x)
% Nonlinear constraints demo with Rosenbrock's function.
% Takes row inputs. If input is not row, will attempt to correct it.
% Syntax:
% y = rosenbrocksfcn(x)
% options = rosenbrocksfcn('init')

if strcmp(x,'init')
    f.Aineq = [] ;
    f.bineq = [] ;
    f.Aeq = [] ;
    f.beq = [] ;
    f.LB = [] ; f.UB = [] ;
    f.nonlcon = 'quadrifolium' ; % Could also use 'heart' or 'unitdisk'
    f.options.PopInitRange = [-2, -2; 2, 2] ;
    f.options.KnownMin = [1,1] ;
    f.options.PopulationSize = 100 ;
    f.options.ConstrBoundary = 'soft' ;
else
    x = reshape(x,1,[]) ;
    if size(x,2) >= 2
%         x1 = x(1:end-1) ; x2 = x(end) ;
        f = 0 ;
        for i = 1:size(x,2)-1
            f = f + (1-x(i))^2 + 100*(x(i+1) - x(i)^2)^2 ;
        end
    else
        error('Rosenbrock''s function requires at least two inputs')
    end
end