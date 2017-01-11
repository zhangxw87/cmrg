function [minloc_row, minloc_col] = CMR_Gross_CV(X,Y,lambdaSet,rhoSet,opts)
%
% Cross validation for finding the optimal regularization parameter for
% CMR_Gross
% 
% 	
% Written by Xiaowei Zhang, June/2016
%
if nargin < 5
    opts = [];
end

n = size(X,1);
Numlbda = length(lambdaSet); Numrho = length(rhoSet);
Error = zeros(Numlbda,Numrho);
num_cv = 5;
% lambda0 = sqrt(log(d)) + sqrt(p);

% start of cross-validation
index = crossvalind('Kfold',n,num_cv);
for i = 1:num_cv
    % get the training and test data for cross-validation
    ind = (index ~= i);
    Xtr = X(ind,:); Ytr = Y(ind,:); % training data
    Xte = X(~ind,:); Yte = Y(~ind,:); % test data
    
    for j=1:Numlbda
        for k=1:Numrho
            opts.lambda = lambdaSet(j);
            opts.rho = rhoSet(k);
            
            [W, ~, ~] = CMR_Gross(Xtr, Ytr, opts);
            
            % error on validation data
            Error(j,k) = Error(j,k) + norm(Yte - Xte*W,'fro') / norm(Yte,'fro');
        end        
    end
end
[~,minloc] = min(Error(:));
[minloc_row, minloc_col] = ind2sub(size(Error), minloc);
