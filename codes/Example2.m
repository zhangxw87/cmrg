function Example2(randVar)
% Example 2 
% covariance matrix of stochastic noise is a scalar matrix 
% gross error rate is 20% size is 5

%% generate data
d = 1000; nTrain = 400; nValid = 400; nTest = 10^4; p = 13; 
sigmaMax = sqrt(2);
opts.NoiseType = 'SameNoise';
opts.gross_rate = 0.2;
opts.grossScale = 5;

lambdaSet = (sqrt(log(d)) + sqrt(p)) * 2.^(linspace(-5,5,21));
rhoSet = 2.^(linspace(-5,5,21));

[Xtr, Ytr, Xvd, Yvd, Xte, Yte, W0, G0, group] = SimuGenerator(d, nTrain, nValid, nTest, p, sigmaMax, randVar, opts);

%% OMR
OMR_opts.GroupType = 'Overlapping';
OMR_opts.group = group;

% select parameter lambda on validation data
OMR_CV_opts = OMR_opts;
OMR_CV_opts.maxitr = 500;
Error = zeros(size(lambdaSet)); 
for i=1:length(lambdaSet)
    OMR_CV_opts.lambda = lambdaSet(i);
    [W, ~] = OMR(Xtr, Ytr, OMR_CV_opts);
    
    % error on validation data
    Error(i) = norm(Yvd-Xvd*W,'fro') / norm(Yvd, 'fro'); 
end
[~,minloc] = min(Error);

OMR_opts.lambda = lambdaSet(minloc);
fprintf('lambda of OMR: %f \n', OMR_opts.lambda);
% OMR_opts.verbose = 1;
[W, ~] = OMR(Xtr, Ytr, OMR_opts);

OMR_ErrPre = norm(Yte - Xte * W,'fro') / norm(Yte,'fro');
OMR_ErrEstW = norm(W - W0,'fro') / norm(W0,'fro');
OMR_RecRateW = sum( sign(W(:)) == sign(W0(:)) ) / numel(W0);

clear Error W; 

%% CMR

CMR_opts.GroupType = 'Overlapping';
CMR_opts.group = group;

% select parameter lambda on validation data
CMR_CV_opts = CMR_opts;
CMR_CV_opts.maxitr = 500;
Error = zeros(size(lambdaSet)); 
for i=1:length(lambdaSet)
    CMR_CV_opts.lambda = lambdaSet(i);
    [W, ~] = CMR(Xtr, Ytr, CMR_CV_opts);
    
    % error on validation data
    Error(i) = norm(Yvd - Xvd * W,'fro') / norm(Yvd, 'fro'); 
end
[~, minloc] = min(Error);

CMR_opts.lambda = lambdaSet(minloc);
fprintf('lambda of CMR: %f \n', CMR_opts.lambda);
% CMR_opts.verbose = 1;
[W, ~] = CMR(Xtr, Ytr, CMR_opts);

CMR_ErrPre = norm(Yte - Xte * W,'fro') / norm(Yte,'fro');
CMR_ErrEstW = norm(W - W0,'fro') / norm(W0,'fro');
CMR_RecRateW = sum( sign(W(:))==sign(W0(:)) ) / numel(W0);

clear Error W; 

%% OMR_Gross
OMR_Gross_opts.GroupType = 'Overlapping';
OMR_Gross_opts.group = group;

% select parameter lambda on validation data
OMR_Gross_CV_opts = OMR_Gross_opts;
OMR_Gross_CV_opts.maxitr = 500;
Error = zeros(length(lambdaSet),length(rhoSet));
for i=1:length(lambdaSet)
    for j=1:length(rhoSet)
        OMR_Gross_CV_opts.lambda = lambdaSet(i);
        OMR_Gross_CV_opts.rho = rhoSet(j);
        [W, ~, ~] = OMR_Gross(Xtr, Ytr, OMR_Gross_CV_opts);
        
        % error on validation data
        Error(i,j) = norm(Yvd - Xvd * W,'fro') / norm(Yvd, 'fro'); 
    end
end
[~,minloc] = min(Error(:));
[minloc_row, minloc_col] = ind2sub(size(Error), minloc);

OMR_Gross_opts.lambda = lambdaSet(minloc_row);
OMR_Gross_opts.rho = rhoSet(minloc_col);
fprintf('lambda of OMR_Gross: %f, rho of OMR_Gross: %f \n', OMR_Gross_opts.lambda, OMR_Gross_opts.rho);
% OMR_Gross_opts.verbose = 1;
[W, G, ~] = OMR_Gross(Xtr, Ytr, OMR_Gross_opts);

OMR_Gross_ErrPre = norm(Yte - Xte * W,'fro') / norm(Yte,'fro');
OMR_Gross_ErrEstW = norm(W - W0,'fro') / norm(W0,'fro');
OMR_Gross_RecRateW = sum( sign(W(:)) == sign(W0(:)) ) / numel(W0);
OMR_Gross_ErrEstG = norm(G - G0,'fro') / norm(G0,'fro');
OMR_Gross_RecRateG = sum( sign(G(:)) == sign(G0(:)) ) / numel(G0);
% clear W G Error;
 
%% CMR_Gross

CMR_Gross_opts.GroupType = 'Overlapping';
CMR_Gross_opts.group = group;

% select parameter lambda on validation data
CMR_Gross_CV_opts = CMR_Gross_opts;
CMR_Gross_CV_opts.maxitr = 500;
Error = zeros(length(lambdaSet),length(rhoSet));
for i=1:length(lambdaSet)
    for j=1:length(rhoSet)
       CMR_Gross_CV_opts.lambda = lambdaSet(i);
       CMR_Gross_CV_opts.rho = rhoSet(j);
        [W, ~, ~] = CMR_Gross(Xtr, Ytr, CMR_Gross_CV_opts);
        
        % error on validation data
        Error(i,j) = norm(Yvd - Xvd * W,'fro') / norm(Yvd, 'fro'); 
   end
end
[~,minloc] = min(Error(:));
[minloc_row, minloc_col] = ind2sub(size(Error), minloc);

CMR_Gross_opts.lambda = lambdaSet(minloc_row);
CMR_Gross_opts.rho = rhoSet(minloc_col);
fprintf('lambda of CMR_Gross: %f, rho of CMR_Gross: %f \n', CMR_Gross_opts.lambda, CMR_Gross_opts.rho);
% CMR_Gross_opts.verbose = 1;
[W, G, ~] = CMR_Gross(Xtr, Ytr, CMR_Gross_opts);

CMR_Gross_ErrPre = norm(Yte - Xte * W,'fro') / norm(Yte,'fro');
CMR_Gross_ErrEstW = norm(W - W0,'fro') / norm(W0,'fro');
CMR_Gross_RecRateW = sum( sign(W(:)) == sign(W0(:)) ) / numel(W0);
CMR_Gross_ErrEstG = norm(G - G0,'fro') / norm(G0,'fro');
CMR_Gross_RecRateG = sum( sign(G(:)) == sign(G0(:)) ) / numel(G0);

% print results
fid = fopen(['../Results/Example2_Results_' num2str(randVar) '.txt'],'a+');
fprintf(fid, ['%8.6f  %8.6f  %8.6f  %8.6f  '...
              '%8.6f  %8.6f  %8.6f  %8.6f  '...
              '%8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  '...
              '%8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f\n'],...
    OMR_ErrPre, OMR_ErrEstW, OMR_RecRateW, OMR_opts.lambda,...
    CMR_ErrPre, CMR_ErrEstW, CMR_RecRateW, CMR_opts.lambda,...
    OMR_Gross_ErrPre, OMR_Gross_ErrEstW, OMR_Gross_RecRateW, OMR_Gross_opts.lambda, OMR_Gross_opts.rho, OMR_Gross_ErrEstG, OMR_Gross_RecRateG,...
    CMR_Gross_ErrPre, CMR_Gross_ErrEstW, CMR_Gross_RecRateW, CMR_Gross_opts.lambda, OMR_Gross_opts.rho, CMR_Gross_ErrEstG, CMR_Gross_RecRateG);
fprintf(fid, '%8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f\n', OMR_Gross_ErrPre, OMR_Gross_ErrEstW, OMR_Gross_RecRateW, OMR_Gross_opts.lambda, OMR_Gross_opts.rho, OMR_Gross_ErrEstG, OMR_Gross_RecRateG);
fclose(fid);
