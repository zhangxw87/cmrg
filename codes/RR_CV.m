function minloc = RR_CV(X, Y, lambdaSet)
% ridge regression with cross-validation
% 	
% Written by Xiaowei Zhang, June/2016
%
[n, d] = size(X); p = size(Y,2);
Numlbda = length(lambdaSet);
Error = zeros(Numlbda,1);
num_cv = 5;
lambda0 = sqrt(log(d)) + sqrt(p);

% start of cross-validation
index = crossvalind('Kfold',n,num_cv);
for i=1:num_cv
    % get the training and test data for cross-validation
    ind = (index~=i);
    Xtr = X(ind,:); Ytr = Y(ind,:); % training data
    Xte = X(~ind,:); Yte = Y(~ind,:); % test data
    
    Xgram = Xtr'*Xtr;
    
    for j=1:Numlbda
        lambda = lambda0 * lambdaSet(j);
        W = (Xgram + lambda*eye(d)) \ Xtr'*Ytr;
        
        % error on validation data
        Error(j) = Error(j) + norm(Yte-Xte*W,'fro') / norm(Yte,'fro');
    end
end
[~,minloc] = min(Error);
