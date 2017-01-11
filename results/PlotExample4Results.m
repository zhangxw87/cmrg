% plot results of Example4
clear; clc; close all;
load('Example4_Results.mat');
gross_rate = 0:0.2:1;

% %% plot Prediction error
% figure(1);
% errorbar(gross_rate, OMR.ErrPre, OMR.ErrPre_std, 'r-.', 'LineWidth', 2);
% hold on;
% errorbar(gross_rate, CMR.ErrPre, CMR.ErrPre_std, 'g--', 'LineWidth', 2);
% errorbar(gross_rate, OMR_Gross.ErrPre, OMR_Gross.ErrPre_std, 'm:', 'LineWidth', 2);
% errorbar(gross_rate, CMR_Gross.ErrPre, CMR_Gross.ErrPre_std, 'b-', 'LineWidth', 2);
% hx = xlabel('$\gamma$'); set(hx,'Interpreter','latex');
% ylabel('Prediction error', 'fontsize', 18);
% hl = legend('OMR', 'CMR', 'OMRG', 'CMRG'); set(hl, 'fontsize', 18);
% axis([-.02 1.02 0 0.9])
% ax = gca;
% set(ax, 'XTick', gross_rate);
% set(ax, 'YTick', 0:0.2:0.8);
% 
% % plot inset figure
% figure(2);
% errorbar(gross_rate, OMR.ErrPre, OMR.ErrPre_std, 'r-.', 'LineWidth', 2);
% hold on;
% errorbar(gross_rate, CMR.ErrPre, CMR.ErrPre_std, 'g--', 'LineWidth', 2);
% errorbar(gross_rate, OMR_Gross.ErrPre, OMR_Gross.ErrPre_std, 'm:', 'LineWidth', 2);
% errorbar(gross_rate, CMR_Gross.ErrPre, CMR_Gross.ErrPre_std, 'b-', 'LineWidth', 2);
% ax = gca;
% set(ax, 'XTick', gross_rate);
% set(ax, 'YTick', 0:0.2:0.8);
% axis([-0.02 0.22 0.1 0.45])
% 
% %% plot Adjusted Prediction error
% figure(3);
% errorbar(gross_rate, OMR.ErrPreAdj, OMR.ErrPreAdj_std, 'r-.', 'LineWidth', 2);
% hold on;
% errorbar(gross_rate, CMR.ErrPreAdj, CMR.ErrPreAdj_std, 'g--', 'LineWidth', 2);
% errorbar(gross_rate, OMR_Gross.ErrPreAdj, OMR_Gross.ErrPreAdj_std, 'm:', 'LineWidth', 2);
% errorbar(gross_rate, CMR_Gross.ErrPreAdj, CMR_Gross.ErrPreAdj_std, 'b-', 'LineWidth', 2);
% hx = xlabel('$\gamma$'); set(hx,'Interpreter','latex');
% ylabel('Adjusted prediction error', 'fontsize', 18);
% hl = legend('OMR', 'CMR', 'OMRG', 'CMRG'); set(hl, 'fontsize', 18);
% axis([-.02 1.02 0 0.9])
% ax = gca;
% set(ax, 'XTick', gross_rate);
% set(ax, 'YTick', 0:0.2:0.8);
% 
% % plot inset figure
% figure(4);
% errorbar(gross_rate, OMR.ErrPreAdj, OMR.ErrPreAdj_std, 'r-.', 'LineWidth', 2);
% hold on;
% errorbar(gross_rate, CMR.ErrPreAdj, CMR.ErrPreAdj_std, 'g--', 'LineWidth', 2);
% errorbar(gross_rate, OMR_Gross.ErrPreAdj, OMR_Gross.ErrPreAdj_std, 'm:', 'LineWidth', 2);
% errorbar(gross_rate, CMR_Gross.ErrPreAdj, CMR_Gross.ErrPreAdj_std, 'b-', 'LineWidth', 2);
% ax = gca;
% set(ax, 'XTick', gross_rate);
% set(ax, 'YTick', 0:0.2:0.8);
% axis([-0.02 0.22 0.05 0.45])
% 
% %% plot Estimation error of W
% figure(5);
% errorbar(gross_rate, OMR.ErrEstW, OMR.ErrEstW_std, 'r-.', 'LineWidth', 2);
% hold on;
% errorbar(gross_rate, CMR.ErrEstW, CMR.ErrEstW_std, 'g--', 'LineWidth', 2);
% errorbar(gross_rate, OMR_Gross.ErrEstW, OMR_Gross.ErrEstW_std, 'm:', 'LineWidth', 2);
% errorbar(gross_rate, CMR_Gross.ErrEstW, CMR_Gross.ErrEstW_std, 'b-', 'LineWidth', 2);
% hx = xlabel('$\gamma$'); set(hx,'Interpreter','latex');
% ylabel('Estimation error of W', 'fontsize', 18);
% hl = legend('OMR', 'CMR', 'OMRG', 'CMRG'); set(hl, 'fontsize', 18);
% axis([-.02 1.02 0 0.9])
% ax = gca;
% set(ax, 'XTick', gross_rate);
% set(ax, 'YTick', 0:0.2:0.8);
% 
% % plot inset figure
% figure(6);
% errorbar(gross_rate, OMR.ErrEstW, OMR.ErrEstW_std, 'r-.', 'LineWidth', 2);
% hold on;
% errorbar(gross_rate, CMR.ErrEstW, CMR.ErrEstW_std, 'g--', 'LineWidth', 2);
% errorbar(gross_rate, OMR_Gross.ErrEstW, OMR_Gross.ErrEstW_std, 'm:', 'LineWidth', 2);
% errorbar(gross_rate, CMR_Gross.ErrEstW, CMR_Gross.ErrEstW_std, 'b-', 'LineWidth', 2);
% ax = gca;
% set(ax, 'XTick', gross_rate);
% set(ax, 'YTick', 0:0.2:0.8);
% axis([-0.02 0.22 0.1 0.45])

%% plot Estimation error of G
% figure(7);
% OMR_Gross.ErrEstG_std(1) = 0.1 * OMR_Gross.ErrEstG_std(1);
% CMR_Gross.ErrEstG_std(1) = 0.1 * CMR_Gross.ErrEstG_std(1);
% 
% errorbar(gross_rate, OMR_Gross.ErrEstG, OMR_Gross.ErrEstG_std, 'm:', 'LineWidth', 2);
% hold on;
% errorbar(gross_rate, CMR_Gross.ErrEstG, CMR_Gross.ErrEstG_std, 'b-', 'LineWidth', 2);
% hx = xlabel('$\gamma$'); set(hx,'Interpreter','latex');
% ylabel('Estimation error of G', 'fontsize', 18);
% hl = legend('OMRG', 'CMRG'); set(hl, 'fontsize', 18);
% axis([-.02 1.02 0 6])
% ax = gca;
% set(ax, 'XTick', gross_rate);

%% plot of Recovery rate of G

figure(8);
errorbar(gross_rate, OMR_Gross.RecRateG, OMR_Gross.RecRateG_std, 'm:', 'LineWidth', 2);
hold on;
errorbar(gross_rate, CMR_Gross.RecRateG, CMR_Gross.RecRateG_std, 'b-', 'LineWidth', 2);
hx = xlabel('$\gamma$'); set(hx,'Interpreter','latex');
ylabel('Indentification rate of G', 'fontsize', 18);
hl = legend('OMRG', 'CMRG'); set(hl, 'fontsize', 18);
% axis([-.02 1.02 0 6])
ax = gca;
set(ax, 'XTick', gross_rate);