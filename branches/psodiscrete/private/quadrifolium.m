function [c,ceq] = quadrifolium(x)
% For demonstrating nonlinear constraints.
%
% Reference:
% Weisstein, Eric W. "Quadrifolium." From MathWorld--A Wolfram Web
% Resource. http://mathworld.wolfram.com/Quadrifolium.html

a = 2 ;
c = (x(1)^2 + x(2)^2)^3 - 4*a*x(1)^2*x(2)^2 ;
ceq = [] ;