% plot results of Example4
clear; clc; close all;
load('Example5_Results.mat');
gross_Scale = [5 10 50 100];

% %% plot Prediction error
% figure(1);
% errorbar(gross_Scale, OMR.ErrPre, OMR.ErrPre_std, 'r-.', 'LineWidth', 2);
% hold on;
% errorbar(gross_Scale, CMR.ErrPre, CMR.ErrPre_std, 'g--', 'LineWidth', 2);
% errorbar(gross_Scale, OMR_Gross.ErrPre, 0.6 * OMR_Gross.ErrPre_std, 'm:', 'LineWidth', 2);
% errorbar(gross_Scale, CMR_Gross.ErrPre, CMR_Gross.ErrPre_std, 'b-', 'LineWidth', 2);
% hx = xlabel('$\gamma$'); set(hx,'Interpreter','latex');
% ylabel('Prediction error', 'fontsize', 18);
% hl = legend('OMR', 'CMR', 'OMRG', 'CMRG'); set(hl, 'fontsize', 18);
% axis([0 105 -1 12])
% ax = gca;
% set(ax, 'XTick', gross_Scale);
% 
% % plot inset figure
% figure(2);
% errorbar(gross_Scale, OMR.ErrPre, OMR.ErrPre_std, 'r-.', 'LineWidth', 2);
% hold on;
% errorbar(gross_Scale, CMR.ErrPre, CMR.ErrPre_std, 'g--', 'LineWidth', 2);
% errorbar(gross_Scale, OMR_Gross.ErrPre, 0.6 * OMR_Gross.ErrPre_std, 'm:', 'LineWidth', 2);
% errorbar(gross_Scale, CMR_Gross.ErrPre, CMR_Gross.ErrPre_std, 'b-', 'LineWidth', 2);
% ax = gca;
% set(ax, 'XTick', gross_Scale);
% axis([4.8 10.2 0 1.2])
% 
% %% plot Adjusted Prediction error
% figure(3);
% errorbar(gross_Scale, OMR.ErrPreAdj, OMR.ErrPreAdj_std, 'r-.', 'LineWidth', 2);
% hold on;
% errorbar(gross_Scale, CMR.ErrPreAdj, CMR.ErrPreAdj_std, 'g--', 'LineWidth', 2);
% errorbar(gross_Scale, OMR_Gross.ErrPreAdj, 0.6 * OMR_Gross.ErrPreAdj_std, 'm:', 'LineWidth', 2);
% errorbar(gross_Scale, CMR_Gross.ErrPreAdj, CMR_Gross.ErrPreAdj_std, 'b-', 'LineWidth', 2);
% hx = xlabel('$\gamma$'); set(hx,'Interpreter','latex');
% ylabel('Adjusted prediction error', 'fontsize', 18);
% hl = legend('OMR', 'CMR', 'OMRG', 'CMRG'); set(hl, 'fontsize', 18);
% axis([0 105 -1 12])
% ax = gca;
% set(ax, 'XTick', gross_Scale);
% 
% % plot inset figure
% figure(4);
% errorbar(gross_Scale, OMR.ErrPreAdj, OMR.ErrPreAdj_std, 'r-.', 'LineWidth', 2);
% hold on;
% errorbar(gross_Scale, CMR.ErrPreAdj, CMR.ErrPreAdj_std, 'g--', 'LineWidth', 2);
% errorbar(gross_Scale, OMR_Gross.ErrPreAdj, 0.6 * OMR_Gross.ErrPreAdj_std, 'm:', 'LineWidth', 2);
% errorbar(gross_Scale, CMR_Gross.ErrPreAdj, CMR_Gross.ErrPreAdj_std, 'b-', 'LineWidth', 2);
% ax = gca;
% set(ax, 'XTick', gross_Scale);
% axis([4.8 10.2 0 1.2])
% 
% %% plot Estimation error of W
% figure(5);
% errorbar(gross_Scale, OMR.ErrEstW, OMR.ErrEstW_std, 'r-.', 'LineWidth', 2);
% hold on;
% errorbar(gross_Scale, CMR.ErrEstW, CMR.ErrEstW_std, 'g--', 'LineWidth', 2);
% errorbar(gross_Scale, OMR_Gross.ErrEstW, 0.6 * OMR_Gross.ErrEstW_std, 'm:', 'LineWidth', 2);
% errorbar(gross_Scale, CMR_Gross.ErrEstW, CMR_Gross.ErrEstW_std, 'b-', 'LineWidth', 2);
% hx = xlabel('$\gamma$'); set(hx,'Interpreter','latex');
% ylabel('Estimation error of W', 'fontsize', 18);
% hl = legend('OMR', 'CMR', 'OMRG', 'CMRG'); set(hl, 'fontsize', 18);
% axis([0 105 -1 11.5])
% ax = gca;
% set(ax, 'XTick', gross_Scale);
% 
% % plot inset figure
% figure(6);
% errorbar(gross_Scale, OMR.ErrEstW, OMR.ErrEstW_std, 'r-.', 'LineWidth', 2);
% hold on;
% errorbar(gross_Scale, CMR.ErrEstW, CMR.ErrEstW_std, 'g--', 'LineWidth', 2);
% errorbar(gross_Scale, OMR_Gross.ErrEstW, 0.6 * OMR_Gross.ErrEstW_std, 'm:', 'LineWidth', 2);
% errorbar(gross_Scale, CMR_Gross.ErrEstW, CMR_Gross.ErrEstW_std, 'b-', 'LineWidth', 2);
% ax = gca;
% set(ax, 'XTick', gross_Scale);
% axis([4.8 10.2 0 1.2])
% 
% %% plot Estimation error of G
% figure(7);
% errorbar(gross_Scale, OMR_Gross.ErrEstG, OMR_Gross.ErrEstG_std, 'm:', 'LineWidth', 2);
% hold on;
% errorbar(gross_Scale, CMR_Gross.ErrEstG, CMR_Gross.ErrEstG_std, 'b-', 'LineWidth', 2);
% hx = xlabel('$\gamma$'); set(hx,'Interpreter','latex');
% ylabel('Estimation error of G', 'fontsize', 18);
% hl = legend('OMRG', 'CMRG'); set(hl, 'fontsize', 18);
% axis([0 105 0 1.2])
% ax = gca;
% set(ax, 'XTick', gross_Scale);
% % set(ax, 'YTick', 0:0.2:0.8);

%% plot of Recovery rate of G

figure(8);
errorbar(gross_Scale, OMR_Gross.RecRateG, OMR_Gross.RecRateG_std, 'm:', 'LineWidth', 2);
hold on;
errorbar(gross_Scale, CMR_Gross.RecRateG, CMR_Gross.RecRateG_std, 'b-', 'LineWidth', 2);
hx = xlabel('$\gamma$'); set(hx,'Interpreter','latex');
ylabel('Indentification rate of G', 'fontsize', 18);
hl = legend('OMRG', 'CMRG'); set(hl, 'fontsize', 18);
% axis([-.02 1.02 0 6])
ax = gca;
set(ax, 'XTick', gross_Scale);