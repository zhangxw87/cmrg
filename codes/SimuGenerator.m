function [Xtr, Ytr, Xvd, Yvd, Xte, Yte, W0, G0, group] = SimuGenerator(d, nTrain, nValid, nTest, p, sigmaMax, randVar, opts)
% generating simulated data
% 	
% Written by Xiaowei Zhang, June/2016
%
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100 * randVar * clock)));

% design matrix
SIGMA = diag(.5 * ones(1,d)) + .5 * ones(d,1) * ones(d,1)';
X = mvnrnd(zeros(1,d),SIGMA, nTrain + nValid + nTest);
Xtr = X(1:nTrain,:);
Xvd = X((nTrain + 1):(nTrain + nValid),:);
Xte = X((nTrain + nValid + 1):end,:);

% construct regression coefficient matrix groups
W0 = zeros(d,p);
NumOLGroup = 19; % can change this value
J = 5 * NumOLGroup + 5;
if J > d
    error('The number of overlapping groups is two large.');
else
    NumGroup = NumOLGroup + d - J;
    group = cell(NumGroup,1);
    for j = 1:NumGroup
        if j <= NumOLGroup
            group{j} = (5 * (j-1) + 1):(5 * (j + 1));
        else
            group{j} = J - NumOLGroup + j;
        end
    end
end
W0(1:J,:) = ((-1) .^ (1:J) .* exp(-(0:J-1) / J))' * ones(1,p);


% noise matrix
if strcmpi(opts.NoiseType, 'SameNoise')
    Z = sigmaMax * randn(nTrain,p);
elseif strcmpi(opts.NoiseType, 'DiffNoise')
    d = sigmaMax .* 2.^(-(0:(p-1))/4);
    Z = randn(nTrain,p) * diag(d);    
else
    error('opts.NoiseType must be either ''SameNoise'' or ''DiffNoise''.');
end


% gross noise
G0 = zeros(nTrain,p);
Numel_gross = floor(opts.gross_rate * nTrain * p);
ind = randperm(nTrain * p, Numel_gross);
G0(ind) = opts.grossScale * sigmaMax * sign(randn(Numel_gross,1));

% response matrix
Ytr = Xtr * W0 + Z + G0;
Yvd = Xvd * W0;
Yte = Xte * W0;
