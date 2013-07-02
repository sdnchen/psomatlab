function f = griewangksfcn(x)
% Griewangk's function.
% 
% Formula: f = 1/4000*sum(x.*x,2) - prod(cos(x./idx),2) + 1
% 
% Global minimum at (0,0).

if strcmp(x,'init')
    f.PopInitRange = [-600; 600] ;
%     f.PopInitRange = [-2; 2] ;
    f.KnownMin = [0,0] ;
    f.Vectorized = 'on' ;
else
%     x = reshape(x,1,[]) ;
    idx = repmat(1:size(x,2),size(x,1),1) ;
    f = 1/4000*sum(x.*x,2) - prod(cos(x./idx),2) + 1 ;
end