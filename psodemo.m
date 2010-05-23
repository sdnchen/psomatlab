function psodemo(DemoMode)
% Runs the PSO on a few demonstration functions, which should be located
% in the <./testfcns> directory.
%
% For less intensive 3D graphics, run with input argument 'fast', viz.:
% >> psodemo('fast')
%
% Note that PSODEMO was included as a quick way to showcase the behaviour
% of the particle swarm on some commonly used optimizer test functions. To
% use the particle swarm algorithm your own problems, use PSO, not PSODEMO.
%
% S. Chen, Dec 2009.
% Available as part of "Another Particle Swarm Toolbox" at:
% http://www.mathworks.com/matlabcentral/fileexchange/25986
% Distributed under BSD license.
%
% See also: PSO

workingdir = pwd ;
testdir = ls('testf*') ;
if ~isempty(testdir), cd(testdir), end

[testfcn,testdir] = uigetfile('*.m','Load demo function for PSO') ;
if ~testfcn
    cd(workingdir)
    return
elseif isempty(regexp(testfcn,'\.m(?!.)','once'))
    error('Test function must be m-file')
else
    cd(testdir)
end

fitnessfcn = str2func(testfcn(1:regexp(testfcn,'\.m(?!.)')-1)) ;
cd(workingdir)

options = fitnessfcn('init') ;

if any(isfield(options,{'options','Aineq','Aeq','LB'}))
    % Then the test function gave us a (partial) problem structure.
    problem = options ;
else
    % Aineq = [1 1] ; bineq = [1.2] ; % Test case for linear constraint
    problem.options = options ;
    problem.Aineq = [] ; problem.bineq = [] ;
    problem.Aeq = [] ; problem.beq = [] ;
    problem.LB = [] ; problem.UB = [] ;
    problem.nonlcon = [] ;
end

problem.fitnessfcn = fitnessfcn ;
problem.nvars = 2 ;

if ~nargin
    problem.options.DemoMode = 'pretty' ;
else
    problem.options.DemoMode = DemoMode ;
end
problem.options.PlotFcns = {@psoplotbestf,@psoplotswarmsurf} ;
% problem.options.VelocityLimit = 0.2 ;
problem.options.HybridFcn = @fmincon ;
% problem.options.Display = 'off' ;

pso(problem)