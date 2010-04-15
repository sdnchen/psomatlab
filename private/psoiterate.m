function state = psoiterate(state,options)
% Updates swarm positions and velocities.

% Weightings for inertia, local, and global influence.
C0 = state.ParticleInertia ;
C1 = options.CognitiveAttraction ; % Local (self best point)
C2 = options.SocialAttraction ; % Global (overall best point)
n = size(state.Population,1) ;
nvars = size(state.Population,2) ;

% Random number seed
R1 = rand(n,nvars) ;
R2 = rand(n,nvars) ;

R1(isinf(state.fLocalBests),:) = 0 ;

% Calculate matrix of velocities state.Velocities for entire population
state.Velocities = C0.*state.Velocities + ...
    C1.*R1.*(state.xLocalBests - state.Population) + ...
    C2.*R2.*(repmat(state.xGlobalBest,n,1) - state.Population) ;

if ~isempty(options.VelocityLimit) && any(isfinite(options.VelocityLimit))
    state.Velocities(state.Velocities > options.VelocityLimit) = ...
        options.VelocityLimit ;
end

state.Population = state.Population + state.Velocities ;

% Update behavioral parameters: reduced inertial term
state.ParticleInertia = 0.9 - 0.2*(state.Generation-1)/ ...
    (options.Generations-1) ;