function [W, Output] = OMR(X, Y, opts)
%
% Solves optimization problem
%  min 0.5*\|Y-XW\|_{F}^2 + lambda*R(W)
%
%
% ADMM method for the following optimization problem of Ordinary Multivariate Regression with gross
% errors
%    W = argmin_{W} 0.5*\|Y - X * W\|_{F}^2 + lambda * R(W)
%
% USAGE:
% [W, Output] = OMR(X, Y, opts)
%
% INPUTS:
%
%       X :   n x d observation matrix, where n is sample size and d is the variable dimension.
%       Y :   n x p response matrix, p is the number of responses. 
%     opts:   opts describes the available options.  If opts is passed as empty,
%             OMR_opts will be called to obtain the default values.
%
% OUTPUTS:
%     W :  d x p the coefficient matrix
% Output:  all other outputs for the optimization problem
%
% 	
% Written by Xiaowei Zhang, June/2016
%==========================================================================

if nargin < 3
    opts = [];
    if nargin < 2;
        error('OMR requires at least two input arguments: X and Y.');
    end
end

% Check compatibility
if size(X,1) ~= size(Y,1)
    error('Dimensions of X and Y are not compatble.');
else
    [n,d] = size(X); p = size(Y,2);
end

% get or check opts
opts = OMR_opts(opts);

t0 = tic; % record time

% Construction of matrix C
if strcmpi(opts.GroupType, 'Overlapping') 
    NumGroup = length(opts.group); 
    gSize = zeros(1+NumGroup,1);
    C = [];
    for i=1:NumGroup
        gSize(i+1) = length(opts.group{i});
        Ctemp = zeros(gSize(i+1),d);
        Ctemp(1:gSize(i+1), opts.group{i}) = eye(gSize(i+1));
        C = [C; Ctemp];
    end
    gSize = cumsum(gSize);
    dnew = gSize(end);
    Ddiag = sum(C); 
elseif strcmpi(opts.GroupType, 'Nonoverlapping')
    C = eye(d); % fro nonoverlapping group sparse or low-rank regularizer
    Ddiag = ones(1,d);
    dnew = d; NumGroup = d; gSize = 0:d;
else
    error('Group sets must be either overlapping or nonoverlapping.')
end

% set parameters
lambda    = opts.lambda;
sigma     = opts.sigma;   
sigmaMax  = opts.sigmaMax;
sigmaMin  = opts.sigmaMin;
tau       = opts.tau; 
tol1      = opts.tol1; 
tol2      = opts.tol2;
maxitr    = opts.maxitr;            
alpha     = opts.alpha;          
delta     = opts.delta;          
verbose   = opts.verbose;

% Compute the inverse of sigma * D + X^T * X
if d < n
    R = chol(X' * X + sigma * diag(Ddiag));
    R = R \ (R \ eye(d))';
else
    Xhat = bsxfun(@rdivide, X, sigma * Ddiag);
    R = chol(eye(n) + sigma * Xhat * diag(Ddiag) * Xhat', 'lower');
    R = R \ Xhat;
    R = diag(1 ./ Ddiag) ./ sigma - R' * R;
end

% Initialization 
Lambda = zeros(dnew,p);
W      = zeros(d,p); 

if verbose
    fprintf('Iteration    primal fun.val. \t dual fun.val. \t primal resd. \t dual resd. \n');
end

% Main loop
Output.flag = 0;
for Numit = 1:maxitr        
    % update Theta
    DTheta = C * W + Lambda / sigma;
    DTheta_norm = zeros(dnew, 1);
    for i = 1:NumGroup
        ind = gSize(i)+1:gSize(i+1);
        temp = norm(DTheta(ind,:),'fro');
        DTheta_norm(ind,:) = max(0, 1 - lambda / (sigma * temp)) * ones(length(ind),1);
    end
    Theta = bsxfun(@times, DTheta_norm, DTheta);
        
    % update W
    W = R * (X' * Y + C' * (sigma * Theta - Lambda));
    
    % update Lambda
    EqResd = Theta - C * W;
    Lambda = Lambda - tau * sigma * EqResd;
    
    if mod(Numit,50) == 0 % compute residuals every 20 iterations
        % primal residual
        etap = norm(EqResd,'fro') / (1 + norm(Theta, 'fro'));
        
        % dual residual        
        etad1 = 0;
        Theta_norm = zeros(NumGroup, 1);
        for i = 1:NumGroup
            ind = gSize(i)+1:gSize(i+1);
            Theta_norm(i) = norm(Theta(ind,:),'fro');
            if Theta_norm(i) == 0
                etad1 = etad1 + max(norm(Lambda(ind,:),'fro') - lambda, 0) ^ 2;
            else
                etad1 = etad1 + norm(Lambda(ind,:) - lambda * Theta(ind,:) / Theta_norm(i),'fro') ^ 2; 
            end
        end
        etad1 = sqrt(etad1);
        
        Z = X' * Y - C' * Lambda;        
        etad2 = norm(X' * X * W - Z, 'fro') / (1 + norm(Z, 'fro'));
        
        etad = max(etad1, etad2);
        
        % duality gap
        Objvalp = 0.5 * norm(Y - X * W, 'fro')^2 + lambda * sum( Theta_norm );
        [U, S, V] = svd(X,0);
        Sdiag = diag(S);
        rankX = rank(X);
        Sdiag(rankX + 1:end) = [];
        S = diag(1 ./ Sdiag);
        Objvald = 0.5 * (norm(Y, 'fro')^2 - norm(U(:,1:rankX) * S * V(:, 1:rankX)' * Z, 'fro')^2);
        Output.dualgap = abs(Objvalp - Objvald) / max([1, abs(Objvalp), abs(Objvald)]);
        
        if verbose
            fprintf(' %i \t\t %8.6f \t %8.6f \t %8.6e \t %8.6e\n', Numit, Objvalp,...
                Objvald, etap, etad);
        end
        % update sigma
        mark = 0;
        if etap > delta * etad
            sigma = min(alpha * sigma, sigmaMax); % increase sigma
            mark = 1;
        elseif etad > delta * etap
            sigma = max(sigma / alpha, sigmaMin); % decrease sigma
            mark = 1;
        end
        
        % if sigma changes, update the inverse of sigma * D + X^T * X
        if mark == 1
            if d < n
                R = chol(X' * X + sigma * diag(Ddiag));
                R = R \ (R \ eye(d))';
            else
                Xhat = bsxfun(@rdivide, X, sigma * Ddiag);
                R = chol(eye(n) + sigma * Xhat * diag(Ddiag) * Xhat', 'lower');
                R = R \ Xhat;
                R = diag(1 ./ Ddiag) ./ sigma - R' * R;
            end
        end
        
        % stopping criterion
        if max(etap, etad) <= tol1 && Output.dualgap <= tol2 
            Output.flag = 1;
            break
        end
    end
end

% Compute W from equality W = D^(-1) * C' * Theta;
W = bsxfun(@ldivide, Ddiag', C' * Theta);

Output.Theta = Theta;
Output.etap = etap;
Output.etad = etad;
Output.Objvalp = Objvalp;
Output.Objvald = Objvald;
Output.time = toc(t0);
Output.Numit = Numit;
Output.sigma = sigma;

% Print results
if (verbose)
    fprintf(1,'\nFinished the main algorithm!\n');
    if Output.flag ==1
        fprintf(1,'The algorithm coverges with tolerance!\n');
    else
        fprintf(1,'The algorithm exits after %2.0d iterations!\n',maxitr);
    end
    fprintf(1,' Primal residual                    = %6.5e\n',Output.etap);
    fprintf(1,' Dual residual                      = %6.5e\n',Output.etad);
    fprintf(1,' Duality gap                        = %6.5e\n',Output.dualgap);
    fprintf(1,' Lgrangian penalty parameter        = %6.5e\n',Output.sigma);
    fprintf(1,' Number of iterations               = %2.0d\n',Output.Numit);
    fprintf(1,' CPU time                           = %3.2e\n',Output.time);
    fprintf(1,'\n');
end

end % end of OMR
