function opts = CMR_opts(opts)

% Options for CMR
%
% If opts is empty upon input, opts will be returned containing the default
% options for CMR.m.  
%
% Alternatively, if opts is passed with some fields already defined, those
% fields will be checked for errors, and the remaining fields will be added
% and initialized to their default values.
%
% Table of Options.  ** indicates default value.
%
% FIELD             DESCRIPTION
% .lambda           regularization parameter for estimator W
%                   ** 1 **
% .GroupType        Type of group sparsity regularization on W
%                   ** 'Nonoverlapping' **
% .group            Set of groups when GroupType = 'Overlapping'
%                   Must be provided.
% .sigma            Initial Lagrangian penalty parameter
%                   ** 1 **
% .sigmaMax         Upper bound of Lagrangian penalty parameter
%                   ** 100 **
% .sigmaMin         Lower bound of Lagrangian penalty parameter
%                   ** 0.01 **
% .tol1             Tolerance on primal and dual residual
%                   ** 1E-4 **
% .tol2             Tolerance on primal and dual gap
%                   ** 1E-4 **
% .tau              Scaling factor of the Larangian multipliers updating step
%                   ** 1.618 **
% .mxitr            Maximum number of iterations.
%                   ** 2000 **
% .delta            Ratio for deciding if updating .sigma, If
%                   primal_residual / dual_residual > delta, increase .sigma; if
%                   dual_residual / primal_residual > delta, decrease .sigma. 
%                   ** 10 **
% .alpha            Scaling factor for increasing or decreasing .sigma
%                   ** 2 **
% .verbose          logical variable deciding if printing intermediate results
%                   ** 0 **
%
% 	
% Written by Xiaowei Zhang, June/2016
%
if isfield(opts,'lambda')
    if opts.lambda < 0
        error('opts.lambda must be nonnegative.');
    end
else
    opts.lambda = 1;
end

if ~isfield(opts, 'GroupType')
    opts.GroupType = 'Nonoverlapping';
end

if strcmpi(opts.GroupType, 'Overlapping')  && (~isfield(opts,'group'))
    error('The group structure must not be empty when overlapping group sparse regularizer is used.')
end

if isfield(opts,'sigma') % initial penalty parameter
    if opts.sigma <= 0
        error('opts.sigma must be positive.');
    end
else
    opts.sigma = 1;
end

if isfield(opts,'sigmaMax')
    if opts.sigmaMax < opts.sigma
        error('opts.sigmaMax must be no less than opts.sigma.');
    end
else
    opts.sigmaMax = 1E+2;
end

if isfield(opts,'sigmaMin')
    if opts.sigma < opts.sigmaMin
        error('opts.sigmaMin must be no more than opts.sigma.');
    end
else
    opts.sigmaMin = 1E-2;
end

if isfield(opts,'tol1')
    if (opts.tol1 <= 0) || (opts.tol1 >= 1)
        error('opts.tol1 is tolerance on primal and dual residual. Should be in (0,1).');
    end
else
    opts.tol1 = 1E-4;
end

if isfield(opts,'tol2')
    if (opts.tol2 <= 0) || (opts.tol2 >= 1)
        error('opts.tol2 is tolerance on primal and dual residual. Should be in (0,1).');
    end
else
    opts.tol2 = 1E-4;
end
 
if isfield(opts,'tau')
    if (opts.tau <= 0) || (opts.tau >= (1+sqrt(5))/2)
        error('If used, opts.tau must be in (0,(1+\sqrt{5})/2).');
    end
else
    opts.tau = 1.618; 
end

if isfield(opts,'maxitr')
    if (opts.maxitr <= 0) || (opts.maxitr ~= floor(opts.maxitr)) 
        error('opts.maxitr should be a positive integer.');
    end
else
    opts.maxitr = 2000;
end

if isfield(opts,'delta')
    if (opts.delta <= 1)
        error('opts.delta should be greater than 1.');
    end
else
    opts.delta = 10;
end

if isfield(opts,'alpha')
    if (opts.alpha <= 1)
        error('opts.alpha should be greater than 1.');
    end
else
    opts.alpha = 2;
end

if ~isfield(opts,'verbose'); 
    opts.verbose = 0; 
end

return
