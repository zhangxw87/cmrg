function [W, G, Output] = CMR_Gross(X, Y, opts)
%
% ADMM method for the following optimization problem of CMR with gross
% errors
%    (W, G) = argmin_{W,G} \|Y - X * W - G\|_{2,1} + lambda * R(W) + rho * \|G\|_1
%
% USAGE:
% [W, G, Output] = CMR_Gross(X, Y, opts)
%
% INPUTS:
%
%       X :   n x d observation matrix, where n is sample size and d is the variable dimension.
%       Y :   n x p response matrix, p is the number of responses. 
%     opts:   opts describes the available options.  If opts is passed as empty,
%             CMR_Gross_opts will be called to obtain the default values.
%
% OUTPUTS:
%     W :  d x p the coefficient matrix
%     G :  n x p the gross noise
% Output:  all other outputs for the optimization problem
%
% 	
% Written by Xiaowei Zhang, June/2016
%==========================================================================

if nargin < 3
    opts = [];
    if nargin < 2;
        error('CMR_Gross requires at least two input arguments: X and Y.');
    end
end

% Check compatibility
if size(X,1) ~= size(Y,1)
    error('Dimensions of X and Y are not compatble.');
else
    [n,d] = size(X); p = size(Y,2);
end

% get or check opts
opts = CMR_Gross_opts(opts);

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
 rho      = opts.rho;
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

% Compute the inverse of D + X^TX
if d < n
    R = chol(X' * X + diag(Ddiag));
    R = R \ (R \ eye(d))';
else
    Xhat = bsxfun(@rdivide, X, Ddiag);
    R = chol(eye(n) + Xhat * diag(Ddiag) * Xhat', 'lower');
    R = R \ Xhat;
    R = diag(1 ./ Ddiag) - R' * R;
end

% Initialization 
Lambda1 = zeros(n,p); 
Lambda2 = zeros(dnew,p); 
G       = zeros(n,p); 
Z       = zeros(n,p); 
Theta   = zeros(dnew,p);
W       = R * (X' * (Y - Z - G) + C' * Theta); 

if verbose
    fprintf('Iteration    primal fun.val. \t dual fun.val. \t primal resd. \t dual resd. \n');
end

% Main loop
Output.flag = 0;
for Numit = 1:maxitr    
    % update Z
    DeltaZ    = Y - X * W -G + Lambda1 / sigma;
    DZdiag    = sqrt(sum(DeltaZ.^2));
    ind       = (DZdiag > 1/sigma);
    Z(:,~ind) = 0;
    Z(:,ind)  = bsxfun( @times, DeltaZ(:,ind), 1 - 1 ./ (sigma * DZdiag(ind)) );
    
    % update Theta
    DTheta = C * W + Lambda2 / sigma;
    DTheta_norm = zeros(dnew, 1);
    for i = 1:NumGroup
        ind = gSize(i)+1:gSize(i+1);
        temp = norm(DTheta(ind,:),'fro');
        DTheta_norm(ind,:) = max(0, 1 - lambda / (sigma * temp)) * ones(length(ind),1);
    end
    Theta = bsxfun(@times, DTheta_norm, DTheta);
    
    % update W
    W = R * (X'*(Y-Z-G) + C'*Theta);
    
    % update G
    DeltaG = Y - X * W -Z + Lambda1 / sigma;
    G      = sign(DeltaG) .* max(abs(DeltaG) - rho / sigma, 0);
    
    % update W
    W = R * (X'*(Y-Z-G) + C'*Theta);
    
    % update Lambda1 and Lambda2
    EqResd1 = Z + X * W + G - Y;
    EqResd2 = Theta - C * W;
    Lambda1 = Lambda1 - tau * sigma * EqResd1;
    Lambda2 = Lambda2 - tau * sigma * EqResd2;
    
    if mod(Numit,50) == 0 % compute residuals every 50 iterations
        % primal residual
        etap = max( norm(EqResd1,'fro') / (1 + norm(Y,'fro')), ...
            norm(EqResd2,'fro') / (1 + norm(Theta, 'fro')) );
        
        % dual residual
        Z_norm = sqrt(sum(Z.^2));
        ind = (Z_norm == 0);
        etad1 = sum( (max(sqrt(sum((Lambda1(:,ind)).^2)) - 1,0)).^2 );
        if any(~ind)
            etad1 = etad1 + norm(Lambda1(:,~ind) - Z(:,~ind) ./ repmat(Z_norm(~ind), n, 1),'fro')^2;
        end
        etad1 = sqrt(etad1);
        
        etad2 = 0;
        Theta_norm = zeros(NumGroup, 1);
        for i = 1:NumGroup
            ind = gSize(i)+1:gSize(i+1);
            Theta_norm(i) = norm(Theta(ind,:),'fro');
            if Theta_norm(i) == 0
                etad2 = etad2 + max(norm(Lambda2(ind,:),'fro') - lambda, 0) ^ 2;
            else
                etad2 = etad2 + norm(Lambda2(ind,:) - lambda * Theta(ind,:) / Theta_norm(i),'fro') ^ 2; 
            end
        end
        etad2 = sqrt(etad2);
        
        ind = (G ~= 0);
        etad3 = sum( (Lambda1(ind) - rho * sign(G(ind))) .^ 2 );
        etad3 = etad3 + sum( max(abs(Lambda1(~ind)) - rho,0) .^ 2 );
        etad3 = sqrt(etad3);
        etad = max([etad1, etad2, etad3]);
        
        % duality gap
        Objvalp = sum(Z_norm) + lambda * sum( Theta_norm ) + rho * sum( abs(G(:)) );
        Objvald = trace(Lambda1' * Y);
        Output.dualgap = abs(Objvalp - Objvald) / max([1, abs(Objvalp), abs(Objvald)]);
        
        if verbose
            fprintf(' %i \t\t %8.6f \t %8.6f \t %8.6e \t %8.6e\n', Numit, Objvalp,...
                Objvald, etap, etad);
        end
        % update sigma
        if etap > delta * etad
            sigma = min(alpha * sigma, sigmaMax); % increase sigma
        elseif etad > delta * etap
            sigma = max(sigma / alpha, sigmaMin); % decrease sigma
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

end % end of CMR_Gross

