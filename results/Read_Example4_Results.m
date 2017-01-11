% Read Example3 Results
clear; clc;
load Example4_Results.mat

for i = 1:90
    str1 = ['Example41_Results_' int2str(i) '.txt'];
    fid1 = fopen(str1);
    if fid1 >= 0
        temp = fscanf(fid1, '%f');
        fclose(fid1);
        temp = reshape(temp, 27, 3)';
        ResMat1 = [ResMat1; temp(1,:)];
        ResMat2 = [ResMat2; temp(2,:)];
        ResMat3 = [ResMat3; temp(3,:)];
    end
    
    str2 = ['Example42_Results_' int2str(i) '.txt'];
    fid2 = fopen(str2);
    if fid2 >= 0
        temp = fscanf(fid2, '%f');
        fclose(fid2);
        temp = reshape(temp, 27, 3)';
        ResMat4 = [ResMat4; temp(1,:)];
        ResMat5 = [ResMat5; temp(2,:)];
        ResMat6 = [ResMat6; temp(3,:)];
    end
end

Mean(1,:) = mean(ResMat1); Std(1,:) = std(ResMat1, 0, 1);
Mean(2,:) = mean(ResMat2); Std(2,:) = std(ResMat2, 0, 1);
Mean(3,:) = mean(ResMat3); Std(3,:) = std(ResMat3, 0, 1);
Mean(4,:) = mean(ResMat4); Std(4,:) = std(ResMat4, 0, 1);
Mean(5,:) = mean(ResMat5); Std(5,:) = std(ResMat5, 0, 1);
Mean(6,:) = mean(ResMat6); Std(6,:) = std(ResMat6, 0, 1);

OMR.ErrPre = Mean(:,2); OMR.ErrPreAdj = Mean(:,3); OMR.ErrEstW = Mean(:,4); 
OMR.ErrPre_std = Std(:,2); OMR.ErrPreAdj_std = Std(:,3); OMR.ErrEstW_std = Std(:,4); 
OMR.RecRateW = Mean(:,5); OMR.lambda = Mean(:,6);
OMR.RecRateW_std = Std(:,5); OMR.lambda_std = Std(:,6);

CMR.ErrPre = Mean(:,7); CMR.ErrPreAdj = Mean(:,8); CMR.ErrEstW = Mean(:,9); 
CMR.ErrPre_std = Std(:,7); CMR.ErrPreAdj_std = Std(:,8); CMR.ErrEstW_std = Std(:,9); 
CMR.RecRateW = Mean(:,10); CMR.lambda = Mean(:,11);
CMR.RecRateW_std = Std(:,10); CMR.lambda_std = Std(:,11);

OMR_Gross.ErrPre = Mean(:,12); OMR_Gross.ErrPreAdj = Mean(:,13); OMR_Gross.ErrEstW = Mean(:,14); 
OMR_Gross.ErrPre_std = Std(:,12); OMR_Gross.ErrPreAdj_std = Std(:,13); OMR_Gross.ErrEstW_std = Std(:,14); 
OMR_Gross.RecRateW = Mean(:,15); OMR_Gross.lambda = Mean(:,16); OMR_Gross.rho = Mean(:,17); 
OMR_Gross.RecRateW_std = Std(:,15); OMR_Gross.lambda_std = Std(:,16); OMR_Gross.rho_std = Std(:,17); 
OMR_Gross.ErrEstG = Mean(:,18); OMR_Gross.RecRateG = Mean(:,19);
OMR_Gross.ErrEstG_std = Std(:,18); OMR_Gross.RecRateG_std = Std(:,19);

CMR_Gross.ErrPre = Mean(:,20); CMR_Gross.ErrPreAdj = Mean(:,21); CMR_Gross.ErrEstW = Mean(:,22); 
CMR_Gross.ErrPre_std = Std(:,20); CMR_Gross.ErrPreAdj_std = Std(:,21); CMR_Gross.ErrEstW_std = Std(:,22); 
CMR_Gross.RecRateW = Mean(:,23); CMR_Gross.lambda = Mean(:,24); CMR_Gross.rho = Mean(:,25); 
CMR_Gross.RecRateW_std = Std(:,23); CMR_Gross.lambda_std = Std(:,24); CMR_Gross.rho_std = Std(:,25);
CMR_Gross.ErrEstG = Mean(:,26); CMR_Gross.RecRateG = Mean(:,27);
CMR_Gross.ErrEstG_std = Std(:,26); CMR_Gross.RecRateG_std = Std(:,27);
