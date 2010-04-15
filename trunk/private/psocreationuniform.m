function state = psocreationuniform(options,nvars)
% Generates uniformly distributed swarm based on options.PopInitRange.

n = options.PopulationSize ;
itr = options.Generations ;

nbrtocreate = n ;
state.Population = zeros(n,nvars) ;
if ~isempty(options.InitialPopulation)
    nbrtocreate = nbrtocreate - size(options.InitialPopulation,1) ;
    state.Population(1:n-nbrtocreate,:) = options.InitialPopulation ;
    disp('Found initial population')
end

nbrtocreate = n ;
% Initial particle velocities
state.Velocities = zeros(n,nvars) ;
if ~isempty(options.InitialVelocities)
    nbrtocreate = nbrtocreate - size(options.InitialVelocities,1) ;
    state.Velocities(1:n-nbrtocreate,:) = options.InitialVelocities ;
    disp('Found initial velocities')
end

% Initialize particle positions
state.Population(n-nbrtocreate+1:n,:) = ...
    repmat(options.PopInitRange(1,:),nbrtocreate,1) + ...
    repmat((options.PopInitRange(2,:) - options.PopInitRange(1,:)),...
    nbrtocreate,1).*rand(nbrtocreate,nvars) ;

% Initialize the global and local fitness to the worst possible
state.fGlobalBest = ones(itr,1)*inf; % Global best fitness score
state.fLocalBests = ones(n,1)*inf ; % Individual best fitness score

% Initialize global and local best positions
state.xGlobalBest = ones(1,nvars)*inf ;
state.xLocalBests = ones(n,nvars)*inf ;