#Introduction
Particle swarm optimization (PSO) is a derivative-free global optimum solver. It is inspired by the surprisingly organized behaviour of large groups of simple animals, such as flocks of birds, schools of fish, or swarms of locusts. The individual creatures, or "particles", in this algorithm are primitive, knowing only four simple things: 1 & 2) their own current location in the search space and fitness value, 3) their previous personal best location, and  4) the overall best location found by all the particles in the "swarm". There are no gradients or Hessians to calculate. Each particle continually adjusts its speed and trajectory in the search space based on this information, moving closer towards the global optimum with each iteration. As seen in nature, this computational swarm displays a remarkable level of coherence and coordination despite the simplicity of its individual particles.

#Ease of Use
If you are already using the Genetic Algorithm (GA) included with MATLAB's Global Optimization Toolbox, then this PSO toolbox will save you a great deal of time. It can be called from the MATLAB command line using the same syntax as the GA, with some additional options specific to PSO. This will allow a high degree of code re-usability between the PSO toolbox and the GA toolbox. Certain GA-specific parameters such as cross-over and mutation functions will obviously not be applicable to the PSO algorithm. However, many of the commonly used options for the Genetic Algorithm Toolbox may be used interchangeably with PSO since they are both iterative population-based solvers. See >> help pso (from the ./psopt directory) for more details.

#Features
  * NEW: support for distributed computing using MATLAB's parallel computing toolbox.
  * Full support for bounded, linear, and nonlinear constraints.
  * Modular and customizable.
  * Binary optimization. See PSOBINARY function for details.
  * Vectorized fitness functions.
  * Solver parameters controlled using 'options' structure similar to existing MATLAB optimization solvers.
  * User-defined custom plots may be written using same template as GA plotting functions.
  * Another optimization solver may be called as a "hybrid function" to refine PSO results.

A demo function is included, with a small library of test functions. To run the demo, from the psopt directory, call >> psodemo with no inputs or outputs.

New features and bug fixes will continue to be released until this is made redundant by the release of an official MATLAB PSO toolbox. Bug reports and feature requests are welcome.

Special thanks to the following people for contributing code and bug fixes:
  * Ben Xin Kang of the University of Hong Kong
  * Christian Hansen of the University of Hannover
  * Erik Schreurs from the MATLAB Central community
  * J. Oliver of Brigham Young University
  * Michael Johnston of the IRIS toolbox
  * Ziqiang (Kevin) Chen

#Bibliography
  * J Kennedy, RC Eberhart, YH Shi. Swarm Intelligence. Academic Press, 2001.
  * Particle Swarm Optimization. http://en.wikipedia.org/wiki/Particle_swarm_optimization
  * RE Perez, K Behdinan. Particle swarm approach for structural design optimization. Computers and Structures 85 (2007) 1579–1588.
  * SM Mikki, AA Kishk. Particle Swarm Optimization: A Physics-Based Approach. Morgan & Claypool, 2008.

#Addendum A
Nonlinear inequality constraints in the form c(x) ≤ 0 and nonlinear equality constraints of the form ceq(x) = 0 have now been fully implemented. The 'penalize' constraint boundary enforcement method is now default. It has been redesigned and tested extensively, and should work with all types of constraints.

See the following document for the proper syntax for defining nonlinear constraint functions: http://www.mathworks.com/help/optim/ug/writing-constraints.html#brhkghv-16.
To see a demonstration of nonlinear inequality constraints using a quadrifolium overlaid on Rosenbrock's function, run PSODEMO and choose 'nonlinearconstrdemo' as the test function.

#Addendum B
See the following guide in the GA toolbox documentation to get started on using the parallel computing toolbox.
http://www.mathworks.com/help/gads/genetic-algorithm-options.html#f17234

#Addendum C
If you are a beginner hoping to learn to use this toolbox for work or school, here are some essential readings:
  * MATLAB's Optimization Toolbox: http://www.mathworks.com/help/optim/index.html
  * MATLAB's Global Optimization Toolbox: http://www.mathworks.com/help/gads/index.html
  * MATLAB's Genetic Algorithm: http://www.mathworks.com/help/gads/genetic-algorithm.html

#Addendum D
There is now a particle swarm optimizer included with the Global Optimization Toolbox. It does not seem to handle constraints at this time. If you have a recent version of the Global Optimization Toolbox installed, you will need to set the path appropriately in your code to use this toolbox.
