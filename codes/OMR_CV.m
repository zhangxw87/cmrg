function minloc = OMR_CV(X, Y, lambdaSet, opts)
%
% Cross validation for finding the optimal regularization parameter for
% OMR
%
% 	
% Written by Xiaowei Zhang, June/2016
%
if nargin < 4
    opts = [];
end

n = size(X,1);
Numlbda = length(lambdaSet);
Error = zeros(Numlbda,1);
num_cv = 5;

% start of cross-validation
index = crossvalind('Kfold',n,num_cv);
for i=1:num_cv
    % get the training and test data for cross-validation
    ind = (index~=i);
    Xtr = X(ind,:); Ytr = Y(ind,:); % training data
    Xte = X(~ind,:); Yte = Y(~ind,:); % test data
    
    for j=1:Numlbda
        opts.lambda = lambdaSet(j);
        [W, ~] = OMR(Xtr, Ytr, opts);
        
        % error on validation data
        Error(j) = Error(j) + norm(Yte - Xte*W,'fro') / norm(Yte,'fro');
    end
end
[~,minloc] = min(Error);
