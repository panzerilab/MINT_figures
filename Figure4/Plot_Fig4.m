
clc, clear, close all; 
load('Results/FIT_results_1410.mat')
rng(0)

% -------------------------------------------------------------------------
% Mean Values  
% -------------------------------------------------------------------------
TE_mean = mean(TE_v,3);

FIT1_mean = mean(FIT1,3);
FIT2_mean = mean(FIT2,3);

FIT_shift_12_1_mean = mean(FIT12_1_shift,3);
FIT_shift_12_2_mean = mean(FIT12_2_shift,3);
FIT_shift_21_2_mean = mean(FIT21_2_shift,3);

% -------------------------------------------------------------------------
% Null Samples 
% -------------------------------------------------------------------------
null_samples = 500;

TE_simple =  btsp_shuffle_helper(TE_sh.simple,null_samples);

FIT1_simple = btsp_shuffle_helper(FIT1_sh.simple,null_samples);
FIT1_cond = btsp_shuffle_helper(FIT1_sh.conditioned,null_samples);

FIT2_simple = btsp_shuffle_helper(FIT2_sh.simple,null_samples);
FIT2_cond = btsp_shuffle_helper(FIT2_sh.conditioned,null_samples);

FIT_shift_12_1_simple = btsp_shuffle_helper(FIT12_1_sh.simple,null_samples);
FIT_shift_12_1_cond = btsp_shuffle_helper(FIT12_1_sh.cond,null_samples);

FIT_shift_12_2_simple = btsp_shuffle_helper(FIT12_2_sh.simple,null_samples);
FIT_shift_12_2_cond = btsp_shuffle_helper(FIT12_2_sh.cond,null_samples);

FIT_shift_21_2_simple = btsp_shuffle_helper(FIT21_2_sh.simple,null_samples);
FIT_shift_21_2_cond = btsp_shuffle_helper(FIT21_2_sh.cond,null_samples);

% -------------------------------------------------------------------------
% Take the element-wise maximum across the two null hypotheses 
% -------------------------------------------------------------------------
FIT1_sh_max = max(FIT1_simple, FIT1_cond);
FIT2_sh_max = max(FIT2_simple, FIT2_cond);

FIT12_1_sh_max = max(FIT_shift_12_1_simple, FIT_shift_12_1_cond);
FIT12_2_sh_max = max(FIT_shift_12_2_simple, FIT_shift_12_2_cond);
FIT21_2_sh_max = max(FIT_shift_21_2_simple, FIT_shift_21_2_cond);

% -------------------------------------------------------------------------
% Identify Cluster for time-delay map 
% -------------------------------------------------------------------------
alpha = 0.01;
cluster_mask_FIT12_1 = (clusterStatistics(FIT_shift_12_1_mean, FIT12_1_sh_max, alpha, (1-alpha),0))';
cluster_mask_FIT12_2 = (clusterStatistics(FIT_shift_12_2_mean, FIT12_2_sh_max, alpha, (1-alpha),0))';
cluster_mask_FIT21_2 = (clusterStatistics(FIT_shift_21_2_mean, FIT21_2_sh_max, alpha, (1-alpha),0))';

% -------------------------------------------------------------------------
% Determine significance thresholds 
% -------------------------------------------------------------------------
sig_thresh_TE = prctile(TE_simple,99.9,3);

sig_thresh_FIT1 = prctile(FIT1_sh_max,99,3);
sig_thresh_FIT2 = prctile(FIT2_sh_max,99,3);

sig_thresh_FIT12_1 = prctile(FIT12_1_sh_max,99,3); 
sig_thresh_FIT12_2 = prctile(FIT12_2_sh_max,99,3);
sig_thresh_FIT21_2 = prctile(FIT21_2_sh_max,99,3);



% -------------------------------------------------------------------------
%% Maximal Information Population
% -------------------------------------------------------------------------
infoS1_max = zeros(4,1);
infoS2_max = zeros(4,1);

InfoS1_mean = mean(infoS1, 3);
InfoS2_mean = mean(infoS2, 3);

[~, indexS1] = max(InfoS1_mean, [], 2);
[~, indexS2] = max(InfoS2_mean, [], 2);

for neuron = 1:nPopulations
    t = indexS1(neuron);
    infoS1_max(neuron,:) = InfoS1_mean(neuron, t);
end 

for neuron = 1:nPopulations
     t = indexS2(neuron);
    infoS2_max(neuron,:) = InfoS2_mean(neuron, t);
end 


% -------------------------------------------------------------------------
%% Plot Maximal Information Population 
% -------------------------------------------------------------------------
for neuron = 1:nPopulations
    figure;
    meanS1 = infoS1_max(neuron,:);
    meanS2 = infoS2_max(neuron,:);
    colorS1 = [123 148 146]/255;
    colorS2 = [190 161 127]/255;
    bar(1, meanS1, 'FaceColor', colorS1); 
    hold on;
    bar(2, meanS2, 'FaceColor', colorS2);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'', ''});
    yLimits = [0, 0.65]; 
    ylim(yLimits);
    set(gca, 'YColor', 'none');
    set(gca, 'YTickLabel', []);
    set(gca, 'Box', 'off');
    saveas(gcf, ['Plots/Information_FITsim_Neuron' num2str(neuron) '.svg']);
    close(gcf);
 end


% -------------------------------------------------------------------------
%% Plot Information Time Course 
% -------------------------------------------------------------------------
figure('Position',[270,197,798,277])
subplot(1,2,1)
plot(squeeze(mean(infoS1,3))','LineWidth',2)
xlabel('t [ms]')
ylabel('Information [bits]')
title('Information S^1')

subplot(1,2,2)
plot(squeeze(mean(infoS2,3))','LineWidth',2)
legend('X_1','X_2','X_3','X_4','AutoUpdate','off')
xlabel('t [ms]')
ylabel('[bits]')
title('Information S^2')

saveas(gcf, 'Plots/Information_time_course.svg')

% -------------------------------------------------------------------------
%% Plot TE and FIT Map
% -------------------------------------------------------------------------
figure('Position',[25,230,1179,296])
subplot(1,3,1)
imagesc(TE_mean)
for i = 1:nPopulations
    for j = 1:nPopulations
        if(TE_mean(i,j)>sig_thresh_TE(i,j))
            text(j,i,'*')
        end
    end
end
xlabel('receiver')
ylabel('sender')
yticks(1:nPopulations)
yticklabels_array = cellstr(num2str((1:nPopulations)', '%d'));
title('TE')
colorbar()

subplot(1,3,2)
imagesc(FIT1_mean)
for i = 1:nPopulations
    for j = 1:nPopulations
        if(FIT1_mean(i,j)>sig_thresh_FIT1(i,j))
            text(j,i,'*')
        end
    end
end
xlabel('receiver')
ylabel('sender')
title('FIT_{S1}')
yticks(1:nPopulations)
yticklabels_array = cellstr(num2str((1:nPopulations)', '%d'));
colorbar()

subplot(1,3,3)
imagesc(FIT2_mean)
for i = 1:nPopulations
    for j = 1:nPopulations
        if(FIT2_mean(i,j)>sig_thresh_FIT2(i,j))
            text(j,i,'*')
        end
    end
end
xlabel('receiver')
ylabel('sender')
yticks(1:nPopulations)
yticklabels_array = cellstr(num2str((1:nPopulations)', '%d'));
title('FIT_{S2}')
colorbar()

saveas(gcf, 'Plots/TE_and_FIT_FITsim.svg')

% -------------------------------------------------------------------------
%% Plot time-delay Map
% -------------------------------------------------------------------------

figure('Position',[25,600,1100,900]);

fontsize=14.5;
fontsize2=15.5;

% -------------------------------------------------------------------------
% Prepare Data 
% -------------------------------------------------------------------------
infoS1_mean_1 = squeeze(mean(infoS1(1,:,:), 3));
infoS1_mean_2 = squeeze(mean(infoS1(2,:,:), 3));

infoS1_mean_1 = infoS1_mean_1(:,1:20);
infoS1_mean_2 = infoS1_mean_2(:,1:20);

SEM_1 = squeeze(std(infoS1(1,1:20,:), 0, 3) / sqrt(size(infoS1, 3)));
SEM_2 = squeeze(std(infoS1(2,1:20,:), 0, 3) / sqrt(size(infoS1, 3)));

FIT12_1_Shift_mean = mean(FIT12_1_shift,3);
y_limits_1 = [min(min(infoS1_mean_1 - SEM_1), min(infoS1_mean_2 - SEM_2)) max(max(infoS1_mean_1 + SEM_1), max(infoS1_mean_2 + SEM_2))+0.1];

infoS2_mean_1 = squeeze(mean(infoS2(1,:,:), 3));
infoS2_mean_2 = squeeze(mean(infoS2(2,:,:), 3));

infoS2_mean_1 = infoS2_mean_1(:,1:20);
infoS2_mean_2 = infoS2_mean_2(:,1:20);

SEM_1_2 = squeeze(std(infoS2(1,1:20,:), 0, 3) / sqrt(size(infoS2, 3)));
SEM_2_2 = squeeze(std(infoS2(2,1:20,:), 0, 3) / sqrt(size(infoS2, 3)));

FIT12_2_Shift_mean = mean(FIT12_2_shift,3);
FIT21_2_Shift_mean = mean(FIT21_2_shift,3);
y_limits_2 = [min(min(infoS2_mean_1 - SEM_1_2), min(infoS2_mean_2 - SEM_2_2)) max(max(infoS2_mean_1 + SEM_1_2), max(infoS2_mean_2 + SEM_2_2))+0.1];
y_limits = [min(y_limits_1(1), y_limits_2(1)), max(y_limits_1(2), y_limits_2(2))];

colorS1 = [123 148 146]/255;
colorS2 = [190 161 127]/255;

minVal = min([min(FIT12_1_Shift_mean(:)), min(FIT12_2_Shift_mean(:)), min(FIT21_2_Shift_mean(:))]);
maxVal = max([max(FIT12_1_Shift_mean(:)), max(FIT12_2_Shift_mean(:)), max(FIT21_2_Shift_mean(:))]);


% -------------------------------------------------------------------------
% Plot information time course
% -------------------------------------------------------------------------
t = tiledlayout(11, 6, 'TileSpacing', 'tight', 'Padding', 'tight');
ax2 = [];
ax1 = []; 

nexttile([3, 3])
plot(squeeze(mean(infoS1,3))','LineWidth',2)
xlabel('t')
ylabel('[bits]', 'FontSize', fontsize, 'FontName', 'Arial')
set(gca, 'YLim', [0,0.7], 'FontSize', fontsize, 'FontName', 'Arial', 'TickDir', 'out', 'Box', 'off')
title('Information S1', 'FontSize',fontsize2)
ax1 = [ax1, gca];

nexttile([3, 3])
plot(squeeze(mean(infoS2,3))','LineWidth',2)
legend('X_1','X_2','X_3','X_4','AutoUpdate','off')
xlabel('t')
set(gca, 'YLim', [0,0.7],'FontSize', fontsize,'TickDir', 'out', 'Box', 'off', 'YTick', [])
title('Information S2', 'FontSize',fontsize2)
ax1 = [ax1, gca];

% -------------------------------------------------------------------------
% Plot 1 - FIT1, MI(X1;S1)
% -------------------------------------------------------------------------
nexttile([2 2])
plot(infoS1_mean_1, 'LineWidth', 1.5, 'Color', colorS1)
hold on
fill([1:length(infoS1_mean_1), length(infoS1_mean_1):-1:1], ...
    [infoS1_mean_1 - SEM_1, fliplr(infoS1_mean_1 + SEM_1)], colorS1, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold off
ylabel('[bits]', 'FontSize', fontsize, 'FontName', 'Arial')
set(gca, 'YLim', y_limits, 'FontSize', fontsize, 'FontName', 'Arial')
title('MI(X_1;S^1)', 'FontSize',fontsize2)
xticks([1, 5:5:15, 20])
set(gca, 'XTickLabel', [], 'TickDir', 'out', 'Box', 'off')
ax2 = [ax2, gca];

% -------------------------------------------------------------------------
% Plot 2 - FIT1, MI(X1;S2)
% -------------------------------------------------------------------------
nexttile([2 2])
plot(infoS2_mean_1, 'LineWidth', 1.5, 'Color', colorS2)
hold on
fill([1:length(infoS2_mean_1), length(infoS2_mean_1):-1:1], ...
    [infoS2_mean_1 - SEM_1_2, fliplr(infoS2_mean_1 + SEM_1_2)], colorS2, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold off
set(gca, 'YLim', y_limits,'YTick', [], 'TickDir', 'out', 'Box', 'off')
title('MI(X_1;S^2)', 'FontSize',fontsize2)
xticks([1, 5:5:15, 20])
set(gca, 'XTickLabel', [], 'TickDir', 'out', 'Box', 'off')
ax2 = [ax2, gca];

% -------------------------------------------------------------------------
% Plot 3 - FIT2, MI(X2;S2)
% -------------------------------------------------------------------------
nexttile([2 2])
plot(infoS2_mean_2, 'LineWidth', 1.5, 'Color', colorS2)
hold on
fill([1:length(infoS2_mean_2), length(infoS2_mean_2):-1:1], ...
    [infoS2_mean_2 - SEM_2_2, fliplr(infoS2_mean_2 + SEM_2_2)], colorS2, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold off
set(gca, 'YLim', y_limits,'YTick', [], 'TickDir', 'out', 'Box', 'off')
title('MI(X_2;S^2)', 'FontSize',fontsize2)
xticks([1, 5:5:15, 20])
set(gca, 'XTickLabel', [], 'TickDir', 'out', 'Box', 'off')
ax2 = [ax2, gca];

% -------------------------------------------------------------------------
% Heatmap - FIT1
% -------------------------------------------------------------------------
nexttile([4 2])
h = imagesc(1:(nTimepoints-10), delay_shifts, FIT12_1_Shift_mean);
set(gca, 'CLim', [minVal, maxVal]);
ylabel(' Δt', 'FontSize', fontsize, 'FontName', 'Arial')
colormap(jet);
yticks([1, 5:5:15, 19])
xticks([1, 5:5:15, 20])
set(gca, 'XTickLabel', [], 'FontSize', fontsize, 'FontName', 'Arial')
title('FIT_{S1}(X_1;X_2)', 'FontSize',fontsize2)
hold on 
alpha_data = ~isnan(FIT12_1_Shift_mean) & cluster_mask_FIT12_1';
set(h, 'AlphaData', alpha_data);
ax2 = [ax2, gca];
hold off

% -------------------------------------------------------------------------
% Heatmap - FIT2
% -------------------------------------------------------------------------
nexttile([4 2])
h = imagesc(1:(nTimepoints-10), delay_shifts, FIT12_2_Shift_mean);
set(gca, 'CLim', [minVal, maxVal]);
colormap(jet);
yticks([1, 5:5:15, 19])
xticks([1, 5:5:15, 20])
set(gca, 'XTickLabel', [], 'YTickLabel', [], 'FontSize', fontsize, 'FontName', 'Arial')
hold on 
alpha_data = ~isnan(FIT12_2_Shift_mean) & cluster_mask_FIT12_2';
title('FIT_{S2}(X_1;X_2)', 'FontSize',fontsize2)
set(h, 'AlphaData', alpha_data);
hold off
ax2 = [ax2, gca];

% -------------------------------------------------------------------------
% Heatmap - FIT3
% -------------------------------------------------------------------------
nexttile([4 2])
h = imagesc(1:(nTimepoints-10), delay_shifts, FIT21_2_Shift_mean);
set(gca, 'CLim', [minVal, maxVal]);
colorbar
set(get(colorbar, 'Title'), 'String', '[bits]', 'FontSize', fontsize, 'FontName', 'Arial');
colormap(jet);
yticks([1, 5:5:15, 19])
xticks([1, 5:5:15, 20])
set(gca, 'XTickLabel', [], 'YTickLabel', [], 'FontSize', fontsize, 'FontName', 'Arial')
alpha_data = ~isnan(FIT21_2_Shift_mean) & cluster_mask_FIT21_2';
title('FIT_{S2}(X_2;X_1)', 'FontSize',fontsize2)
set(h, 'AlphaData', alpha_data);
hold off
ax2 = [ax2, gca];


% -------------------------------------------------------------------------
% Plot4 - FIT1, MI(X2;S1)
% -------------------------------------------------------------------------
nexttile([2 2])
plot(infoS1_mean_2, 'LineWidth', 1.5, 'Color', colorS1)
hold on
fill([1:length(infoS1_mean_2), length(infoS1_mean_2):-1:1], ...
    [infoS1_mean_2 - SEM_2, fliplr(infoS1_mean_2 + SEM_2)], colorS1, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold off
xlabel('t', 'FontSize', fontsize, 'FontName', 'Arial')
ylabel('[bits]', 'FontSize', fontsize, 'FontName', 'Arial')
set(gca, 'YLim', y_limits, 'FontSize', fontsize, 'FontName', 'Arial')
title('MI(X_2;S^1)', 'FontSize',fontsize2)
xticks([1, 5:5:15, 20])
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', fontsize)
ax2 = [ax2, gca];

% -------------------------------------------------------------------------
% Plot 5 - FIT2, MI(X2;S2)
% -------------------------------------------------------------------------
nexttile([2 2])
plot(infoS2_mean_2, 'LineWidth', 1.5,'Color', colorS2)
hold on
fill([1:length(infoS2_mean_2), length(infoS2_mean_2):-1:1], ...
    [infoS2_mean_2 - SEM_2_2, fliplr(infoS2_mean_2 + SEM_2_2)], colorS2, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold off
xlabel('t', 'FontSize', fontsize, 'FontName', 'Arial')
set(gca, 'YLim', y_limits,'YTick', [], 'TickDir', 'out', 'Box', 'off')
title('MI(X_2;S^2)', 'FontSize',fontsize2)
xticks([1, 5:5:15, 20])
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', fontsize)
ax2 = [ax2, gca];

% -------------------------------------------------------------------------
% Plot 5 - FIT3, MI(X1;S2)
% -------------------------------------------------------------------------
nexttile([2 2])
plot(infoS2_mean_1, 'LineWidth', 1.5,'Color', colorS2)
hold on
fill([1:length(infoS2_mean_1), length(infoS2_mean_1):-1:1], ...
    [infoS2_mean_1 - SEM_1_2, fliplr(infoS2_mean_1 + SEM_1_2)], colorS2, ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold off
xlabel('t', 'FontSize', fontsize, 'FontName', 'Arial')
set(gca, 'YLim', y_limits,'YTick', [], 'TickDir', 'out', 'Box', 'off')
title('MI(X_1;S^2)', 'FontSize',fontsize2)
xticks([1, 5:5:15, 20])
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', fontsize)
ax2 = [ax2, gca];
linkaxes(ax2, 'x');
xlim(ax2, [1 20])
xticks(ax2, [1, 5:5:15, 20])
% save
saveas(gcf, 'Plots/Combined_DelaySweep_FIT1_FIT2.svg');


% -------------------------------------------------------------------------
% Helper Function 
% -------------------------------------------------------------------------
function data_out = btsp_shuffle_helper(data,n_boot)
    data_in = permute(data,[3 1 2 4]);
    data_out = create_NullDistribution_groupLevel(data_in,n_boot);
end