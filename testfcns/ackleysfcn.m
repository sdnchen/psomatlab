function f = ackleysfcn(x)
% Ackley's Function.

if strcmp(x,'init')
    f.PopInitRange = [-32.768; 32.768] ;
    f.KnownMin = [0,0] ; % For plotting only
    f.Aineq = [] ;
    f.bineq = [] ;
    f.Aeq = [] ;
    f.beq = [] ;
    f.LB = [] ; f.UB = [] ;
    f.LB = [-10 -10] ;
    f.UB = [-1 -1] ;
    f.nonlcon = [] ;
    f.options.DemoMode = 'pretty' ;
else
    x = reshape(x,1,[]) ;
    a = 20 ;
    b = 0.2 ;
    n = size(x,2) ;
    c = 2*pi ;
    f = -a*exp(-b*sqrt(1/n*(x*x'))) - exp(1/n*(cos(c*x)*cos(c*x'))) + ...
        a + exp(1) ;
end