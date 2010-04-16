function psobinary(fitnessfcn,nvars,options)

if ~exist('options','var') % Set default options
    options = struct ;
end % if ~exist

Vmax = 4 ;

cd('../psopt')
options = psooptimset(options) ;
cd('../psodiscrete')

n = options.PopulationSize ;
itr = options.Generations ;

 % Initialize Population
x = randi([0 1],n,nvars) ;
v = zeros(n,nvars) ;
if ~isempty(options.InitialVelocities)
    v = options.InitialVelocities ;
end

idx_ngbrs = 1:n ; % Index of neighborhood
p = x ; % Matrix with local best genomes along its rows.
f = inf*ones(size(p,1),1) ;
localbestf = inf*ones(size(p,1),1) ;
state.ParticleInertia = 0.9 ;
% bestf = inf ;

for k = 1:itr
    % Cycle through population
    R1 = rand(n,nvars) ; % Every iteration, regenerate random number array
    for i = 1:n
        % Evaluate fitness function
        f(i) = fitnessfcn(x(i,:)) ;
        
        % Check new fitness values against local bests
        if f < localbestf(i)
            localbestf(i) = f ;
            p(i,:) = x(i,:) ; % p(i,:) is the local best so far
        end % if
        
        % Get current global best genome
        g = i ; % Index of the global best genome
        for j = idx_ngbrs
            if f(j) < p(g)
                g = j ; % g is index of best performer in the neighborhood
            end % if f(j)
        end % for j
        
        v(i,:) = state.ParticleInertia*v(i,:) + ...
            options.CognitiveAttraction*(p(i,:) - x(i,:)) + ...
            options.SocialAttraction*(p(g,:) - x(i,:)) ;
        for d = 1:nvars
            if abs(v(i,d)) > Vmax, v(i,d) = sign(v(i,d))*Vmax ; end
            if R1(i,d) < sigmoid(v(i,d))
                x(i,d) = 1 ;
            else
                x(i,d) = 0 ;
            end % if rho
        end % for d
        
        state.ParticleInertia = 0.9 - 0.2*(state.Generation-1)/ ...
            (options.Generations-1) ;
    end % for i
    
    % Check convergence here
end

function s = sigmoid(v)
% Sigmoid function

s = 1/(1+exp(v)) ;