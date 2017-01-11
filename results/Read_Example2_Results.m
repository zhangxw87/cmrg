% Read Example2 Results
clear; clc;
load Example2_Results.mat

for i = 1:40
    str = ['Example2_NewResults_' int2str(i) '.txt'];
    fid = fopen(str);
    if fid >= 0
        temp = fscanf(fid, '%f');
        fclose(fid);
        ResMat = [ResMat; temp(1:22)'];
    end   
end

Mean = mean(ResMat); Std = std(ResMat, 0, 1);

OMR.ErrPre = Mean(:,1); OMR.ErrEstW = Mean(:,2); OMR.RecRateW = Mean(:,3); OMR.lambda = Mean(:,4);
OMR.ErrPre_std = Std(:,1); OMR.ErrEstW_std = Std(:,2); OMR.RecRateW_std = Std(:,3); OMR.lambda_std = Std(:,4);

CMR.ErrPre = Mean(:,5); CMR.ErrEstW = Mean(:,6); CMR.RecRateW = Mean(:,7); CMR.lambda = Mean(:,8);
CMR.ErrPre_std = Std(:,5); CMR.ErrEstW_std = Std(:,6); CMR.RecRateW_std = Std(:,7); CMR.lambda_std = Std(:,8);

OMR_Gross.ErrPre = Mean(:,9); OMR_Gross.ErrEstW = Mean(:,10); OMR_Gross.RecRateW = Mean(:,11); 
OMR_Gross.ErrPre_std = Std(:,9); OMR_Gross.ErrEstW_std = Std(:,10); OMR_Gross.RecRateW_std = Std(:,11); 
OMR_Gross.lambda = Mean(:,12); OMR_Gross.rho = Mean(:,13); OMR_Gross.ErrEstG = Mean(:,14); OMR_Gross.RecRateG = Mean(:,15);
OMR_Gross.lambda_std = Std(:,12); OMR_Gross.rho_std = Std(:,13); OMR_Gross.ErrEstG_std = Std(:,14); OMR_Gross.RecRateG_std = Std(:,15);

CMR_Gross.ErrPre = Mean(:,16); CMR_Gross.ErrEstW = Mean(:,17); CMR_Gross.RecRateW = Mean(:,18); 
CMR_Gross.ErrPre_std = Std(:,16); CMR_Gross.ErrEstW_std = Std(:,17); CMR_Gross.RecRateW_std = Std(:,18); 
CMR_Gross.lambda = Mean(:,19); CMR_Gross.rho = Mean(:,20); CMR_Gross.ErrEstG = Mean(:,21); CMR_Gross.RecRateG = Mean(:,2);
CMR_Gross.lambda_std = Std(:,19); CMR_Gross.rho_std = Std(:,20); CMR_Gross.ErrEstG_std = Std(:,21); CMR_Gross.RecRateG_std = Std(:,22);