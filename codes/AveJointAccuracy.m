function MeanAcc = AveJointAccuracy(Xte, W, Yte)
% Computes average joint mean error
% 	
% Written by Xiaowei Zhang, June/2016
%
NumJoint = size(Yte, 2) / 3;
MeanAcc = zeros(1, NumJoint);
Ypred = Xte * W;
for i = 1:NumJoint
    temp = Ypred(:, (3 * (i-1) + 1):(3 * i)) - Yte(:, (3 * (i-1) + 1):(3 * i));
    MeanAcc(i) = mean( sqrt( sum(temp .^ 2, 2) ) );
end
