function [X, W, opt, time, iter, dGap] = QUIC(mode, S, L, tol, msg, maxIter, X0, W0)
% [X W opt time iter dGap] = QUIC(mode, ...)
% [X W opt time iter dGap] = QUIC("default", S, L, tol, msg, ...
%                                 maxIter, X0, W0)
% [X W opt time iter dGap] = QUIC("path", S, L, path, tol, msg, ...
%                                 maxIter, X0, W0)
% [X W opt time iter dGap] = QUIC("trace", S, L, tol, msg, ...
%                                 maxIter, X0, W0)
% Arguments:
% mode      selects computation
% S         empirical covariance matrix
% L         matrix of the regularization parameters
% path      a vector of values used for scaling L in path mode
% tol       convergence threshold
% msg       verbosity of print statistics
% maxIter   maximum number of Newton iterations to execute
% X0        initial value for the inverse covariance matrix
% W0        initial value for the covariance matrix
%
% Returned values:
% X         inverse covariance matrix, in path mode a vector of matrices
%           the solution of
%             min_X f(X), where f(X) = -logdet(X) + trace(XS) + |Lambda.*X|
% W         inverse of X, in path mode a vector of matrices
% opt       value of f(X) at optimum found; an array in "path" and "trace"
%           modes
% time      cpu time; an array in "path" and "trace" modes
% iter      the number of Newton iterations executed in "default" or
%           "trace" mode, an array of values in "path" mode
% dGap      duality gap
%      
% The "trace" mode allows the computation of approximations to the optimal
% value and the cputime used to acquire them along the way.  Useful for
% plotting graphs showing the speed of convergence.
    
disp('Please compile the executable with make QUIC.[mex|mexglx|mexa64]!');

                                         
