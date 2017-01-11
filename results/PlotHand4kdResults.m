% plot Average Joint Error
clear;clc
load('Hand4kd_FinalResults2.mat');

Method = {'OMR', 'CMR', 'OMRG', 'CMRG', 'DHand', 'AlexNet'};
JointError = [OMR_Output.AveJointAccuracy; CMR_Output.AveJointAccuracy; ...
    OMR_Gross_Output.AveJointAccuracy; CMR_Gross_Output.AveJointAccuracy;
    DHand_Output.AveJointAccuracy; AlexNet_Output.AveJointAccuracy]';
h = boxplot(JointError, Method);
set(h,{'linew'},{2})
set(findobj(gca,'Type','text'),'FontSize',15)
ylabel('Average joint error e (mm)')