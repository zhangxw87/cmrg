% Experiment on synthetic hand pose depth images with 20% gross error
clear; clc;
gross_rate = 0.2;

%% generate data
load('../Datasets/syn_hand_cvprw_features/preds_4kd.mat');
Xtr = preds(:, splits == 1); Ytr = hand_xyz(:, splits == 1);
Xte = preds(:, splits ~= 1); Yte = hand_xyz(:, splits ~= 1);
Xtr = Xtr'; Ytr = Ytr'; Xte = Xte'; Yte = Yte';
nValid = 20000; % number of validation data
Xvd = Xtr(1:nValid,:); Yvd = Ytr(1:nValid, :);
Xtr(1:nValid,:) = []; Ytr(1:nValid, :) = [];

clear preds hand_xyz splits;

% add gross noise to Ytr
Ytr_wg = Ytr; 
NumGross = numel(Ytr) * gross_rate;
index = randperm(numel(Ytr));
Ytr_wg(index(1:NumGross)) = 0;
G0 = Ytr_wg - Ytr;

lambdaSet = 10.^(linspace(-5,5,21));
rhoSet = 10.^(linspace(-5,5,21));

%% OMR
% select parameter lambda on validation data
OMR_CV_opts.maxitr = 200;
Error = zeros(size(lambdaSet)); 
for i=1:length(lambdaSet)
    OMR_CV_opts.lambda = lambdaSet(i);
    [W, ~] = OMR(Xtr, Ytr_wg, OMR_CV_opts);
    
    % error on validation data
    Error(i) = sum( AveJointAccuracy(Xvd, W, Yvd) ); 
end
[~,minloc] = min(Error);

OMR_opts.lambda = lambdaSet(minloc);
fprintf('lambda of OMR: %f \n', OMR_opts.lambda);
OMR_opts.verbose = 1;
[W_OMR, ~] = OMR(Xtr, Ytr_wg, OMR_opts);

OMR_Output.AveJointAccuracy = AveJointAccuracy(Xte, W_OMR, Yte);

clear Error; 
save Syn_hand4kd_Results1.mat W_OMR OMR_Output OMR_opts gross_rate;

%% CMR
% select parameter lambda on validation data
CMR_CV_opts.maxitr = 200;
Error = zeros(size(lambdaSet)); 
for i=1:length(lambdaSet)
    CMR_CV_opts.lambda = lambdaSet(i);
    [W, ~] = CMR(Xtr, Ytr_wg, CMR_CV_opts);
    
    % error on validation data
    Error(i) = sum( AveJointAccuracy(Xvd, W, Yvd) );
end
[~, minloc] = min(Error);

CMR_opts.lambda = lambdaSet(minloc);
fprintf('lambda of CMR: %f \n', CMR_opts.lambda);
CMR_opts.verbose = 1;
[W_CMR, ~] = CMR(Xtr, Ytr_wg, CMR_opts);

CMR_Output.AveJointAccuracy = AveJointAccuracy(Xte, W_CMR, Yte);

clear Error W; 
save Syn_hand4kd_Results1.mat -append W_CMR CMR_Output CMR_opts;

%% OMR_Gross
% select parameter lambda on validation data
OMR_Gross_CV_opts.maxitr = 200;
Error = zeros(length(lambdaSet),length(rhoSet));
for i=1:length(lambdaSet)
    for j=1:length(rhoSet)
        OMR_Gross_CV_opts.lambda = lambdaSet(i);
        OMR_Gross_CV_opts.rho = rhoSet(j);
        [W, ~, ~] = OMR_Gross(Xtr, Ytr_wg, OMR_Gross_CV_opts);
        
        % error on validation data
        Error(i,j) = sum( AveJointAccuracy(Xvd, W, Yvd) );
    end
end
[~,minloc] = min(Error(:));
[minloc_row, minloc_col] = ind2sub(size(Error), minloc);

OMR_Gross_opts.lambda = lambdaSet(minloc_row);
OMR_Gross_opts.rho = rhoSet(minloc_col);
fprintf('lambda of OMR_Gross: %f, rho of OMR_Gross: %f \n', OMR_Gross_opts.lambda, OMR_Gross_opts.rho);

OMR_Gross_opts.verbose = 1;
[W_OMRG, G_OMRG, ~] = OMR_Gross(Xtr, Ytr_wg, OMR_Gross_opts);

OMR_Gross_Output.AveJointAccuracy = AveJointAccuracy(Xte, W_OMRG, Yte);
OMR_Gross_Output.ErrEstG = norm(G_OMRG - G0,'fro') / (1 + norm(G0,'fro'));
OMR_Gross_Output.RecRateG = sum( sign(G_OMRG(:)) == sign(G0(:)) ) / numel(G0);

clear Error;
save Syn_hand4kd_Results1.mat -append W_OMRG G_OMRG OMR_Gross_Output OMR_Gross_opts;
 
%% CMR_Gross
% select parameter lambda on validation data
CMR_Gross_CV_opts.maxitr = 200;
Error = zeros(length(lambdaSet),length(rhoSet));
for i=1:length(lambdaSet)
    for j=1:length(rhoSet)
        CMR_Gross_CV_opts.lambda = lambdaSet(i);
        CMR_Gross_CV_opts.rho = rhoSet(j);
        [W, ~, ~] = CMR_Gross(Xtr, Ytr_wg, CMR_Gross_CV_opts);
        
        % error on validation data
        Error(i,j) = sum( AveJointAccuracy(Xvd, W, Yvd) );
    end
end
[~,minloc] = min(Error(:));
[minloc_row, minloc_col] = ind2sub(size(Error), minloc);

CMR_Gross_opts.lambda = lambdaSet(minloc_row);
CMR_Gross_opts.rho = rhoSet(minloc_col);
fprintf('lambda of CMR_Gross: %f, rho of CMR_Gross: %f \n', CMR_Gross_opts.lambda, CMR_Gross_opts.rho);

CMR_Gross_opts.verbose = 1;
[W_CMRG, G_CMRG, ~] = CMR_Gross(Xtr, Ytr_wg, CMR_Gross_opts);

CMR_Gross_Output.AveJointAccuracy = AveJointAccuracy(Xte, W_CMRG, Yte);
CMR_Gross_Output.ErrEstG = norm(G_CMRG - G0,'fro') / (1 + norm(G0,'fro'));
CMR_Gross_Output.RecRateG = sum( sign(G_CMRG(:)) == sign(G0(:)) ) / numel(G0);

save Syn_hand4kd_Results1.mat -append W_CMRG G_CMRG CMR_Gross_Output CMR_Gross_opts;
exit
