function plot_fig3reduced(filename1, filename2, filename3, filename4, output_filename)
load(filename1)
% figureproperties()
xvals = -1:0.1:1;
line_shift = 0.0;

slope_s_real = -1;

theta = zeros(2,1);
estim_theta = zeros(2,1);

hf=figure();
hf.Units = 'centimeters';
hf.Position = [0 0 16 10];

subjects_list = [5,6,14,16,24,27,28,29,31,32,33,34];
tiledlayout(2,3)
nexttile([1,3])
% 
% load(['.\data\c_cleaned\info_analysis_with_population_sizebtsp30angle.mat'], ...
load(['.\data\c_cleaned\info_analysis_with_population_sizeupto192folds.mat'], ...
    "avgSI_popbtsp_partitions_pop_size_ind", ...
    "avgCI_popbtsp_partitions_pop_size_ind",...
    "avgII_popbtsp_partitions_pop_size_ind", ...
    "avgSI_popnaive_partitions_pop_size_ind", ...
    "avgCI_popnaive_partitions_pop_size_ind",...
    "avgII_popnaive_partitions_pop_size_ind", ...
    "Nsessions","sample_combinations_pop_size","SI_CI_info_ind")
load(['.\data\c_cleaned\info_analysis_with_population_sizebtsp30angle.mat'], ...
    "avgbetastimchoice_partitions_pop_size_ind");

% load(['.\data\c_cleaned\info_analysis_with_population_sizeupto192folds.mat'], ...
load(['.\data\c_cleaned\info_analysis_with_population_sizebtsp30direct.mat'], ...
    "avgSI_popnaive_partitions_pop_size_ind_direct", ...
    "avgCI_popnaive_partitions_pop_size_ind_direct",...
    "avgII_popnaive_partitions_pop_size_ind_direct", ...
    "avgSI_popbtsp_partitions_pop_size_ind_direct", ...
    "avgCI_popbtsp_partitions_pop_size_ind_direct",...
    "avgII_popbtsp_partitions_pop_size_ind_direct");
% load(['.\data\c_cleaned\info_analysis_with_population_sizeupto192folds.mat'], ...
load(['.\data\c_cleaned\Results_current_analysis_info_analysis_with_population_sizebtsp30direct.mat'], ...
    'SI_CI_info_ind_btsp');

pop_size_ALL = 1:19;

% Pool the data by taking the average across repetitions for each session and then find the mean and SEM across sessions
[avgSI_pop_partitions_pop_size_ind_across_subj_pairs] = pool_infodata_to_plot_as_info_vs_popsize(avgSI_popbtsp_partitions_pop_size_ind,Nsessions,sample_combinations_pop_size);
[avgCI_pop_partitions_pop_size_ind_across_subj_pairs] = pool_infodata_to_plot_as_info_vs_popsize(avgCI_popbtsp_partitions_pop_size_ind,Nsessions,sample_combinations_pop_size);
[avgII_pop_partitions_pop_size_ind_across_subj_pairs] = pool_infodata_to_plot_as_info_vs_popsize(avgII_popbtsp_partitions_pop_size_ind,Nsessions,sample_combinations_pop_size);
[avgang_pop_partitions_pop_size_ind_across_subj_pairs] = pool_infodata_to_plot_as_info_vs_popsize(avgbetastimchoice_partitions_pop_size_ind,Nsessions,sample_combinations_pop_size);

avgII_pop_partitions_pop_size_ind_across_subj_pairs(:,1,:) = 5*avgII_pop_partitions_pop_size_ind_across_subj_pairs(:,1,:);
avgSI_pop_partitions_pop_size_ind_across_subj_pairs(:,1,:) = 5*avgSI_pop_partitions_pop_size_ind_across_subj_pairs(:,1,:);
avgCI_pop_partitions_pop_size_ind_across_subj_pairs(:,1,:) = 5*avgCI_pop_partitions_pop_size_ind_across_subj_pairs(:,1,:);

avgII_pop_partitions_pop_size_ind_across_subj_pairs(:,19,:) = 5*avgII_pop_partitions_pop_size_ind_across_subj_pairs(:,19,:);
avgSI_pop_partitions_pop_size_ind_across_subj_pairs(:,19,:) = 5*avgSI_pop_partitions_pop_size_ind_across_subj_pairs(:,19,:);
avgCI_pop_partitions_pop_size_ind_across_subj_pairs(:,19,:) = 5*avgCI_pop_partitions_pop_size_ind_across_subj_pairs(:,19,:);

meanSI = mean(avgSI_pop_partitions_pop_size_ind_across_subj_pairs,1);
stdSI = std(avgSI_pop_partitions_pop_size_ind_across_subj_pairs)/sqrt(length(Nsessions)-1);
meanCI = mean(avgCI_pop_partitions_pop_size_ind_across_subj_pairs,1);
stdCI = std(avgCI_pop_partitions_pop_size_ind_across_subj_pairs)/sqrt(length(Nsessions)-1);
meanII = mean(avgII_pop_partitions_pop_size_ind_across_subj_pairs,1);
stdII = std(avgII_pop_partitions_pop_size_ind_across_subj_pairs)/sqrt(length(Nsessions)-1);

avgSI_pop_partitions_pop_size_ind_across_subj_pairs_direct = squeeze(avgSI_popbtsp_partitions_pop_size_ind_direct(subjects_list,:,1:20));
avgCI_pop_partitions_pop_size_ind_across_subj_pairs_direct = squeeze(avgCI_popbtsp_partitions_pop_size_ind_direct(subjects_list,:,1:20));
avgII_pop_partitions_pop_size_ind_across_subj_pairs_direct = squeeze(avgII_popbtsp_partitions_pop_size_ind_direct(subjects_list,:,1:20));
meanSI_direct = mean(avgSI_pop_partitions_pop_size_ind_across_subj_pairs_direct,'all');
stdSI_direct = std(avgSI_pop_partitions_pop_size_ind_across_subj_pairs_direct,[],'all')/sqrt(length(Nsessions)-1);
meanCI_direct = mean(avgCI_pop_partitions_pop_size_ind_across_subj_pairs_direct,'all');
stdCI_direct = std(avgCI_pop_partitions_pop_size_ind_across_subj_pairs_direct,[],'all')/sqrt(length(Nsessions)-1);
meanII_direct = mean(avgII_pop_partitions_pop_size_ind_across_subj_pairs_direct,'all');
stdII_direct = std(avgII_pop_partitions_pop_size_ind_across_subj_pairs_direct,[],'all')/sqrt(length(Nsessions)-1);

transparency_areas = 0.1;
purple = [178 87 189] ./ 255;
brown = [162 117 76] ./ 255;
behav_perf = yline(nanmean(SI_CI_info_ind_btsp),'k-.','LineWidth',2);hold on
errorbar(pop_size_ALL, meanSI(pop_size_ALL), stdSI(pop_size_ALL),'.',"Markersize",16,"Color",brown); hold on;
errorbar(pop_size_ALL, meanCI(pop_size_ALL), stdCI(pop_size_ALL),'.',"Markersize",16,"Color",purple); hold on;
errorbar(pop_size_ALL, meanII(pop_size_ALL), stdII(pop_size_ALL),'.',"Markersize",16,"Color",[107 136 98]./255); hold on;
errorbar(-1, meanSI_direct, stdSI_direct,'.',"Markersize",16,"Color",brown); hold on;
errorbar(-1, meanCI_direct, stdCI_direct,'.',"Markersize",16,"Color",purple); hold on;
errorbar(-1, meanII_direct, stdII_direct,'.',"Markersize",16,"Color",[107 136 98]./255); hold on;
xline(0)
box on

t = axis;
% text(.55*t(2), 1.2*t(4),'Information scaling of intersection information', 'FontWeight', 'Bold', 'FontSize', 8, HorizontalAlignment='center')
%title({'Information scaling of'; 'intersection information'}, 'FontWeight', 'Bold')
legend_names = {'I(S,C)',
    sprintf('I(S,%s)', char(348)),
    sprintf('I(C,%s)', char(264)),
    sprintf('II(S,[%s,%s],C)', char(348), char(264))};
leg = legend(legend_names,...
    'Location','southeast'); 
leg.ItemTokenSize = [15,9];
ylabel('Information [bits]');xlabel('Population size');
xlim([-2 19]);
xticks(-1:19);
xticklabels({'Dr' 0:19});
ylim([0 0.43])

yticks(linspace(0,.5, 6))
% Plot the information scaling curves as a function of the population size
nexttile;
image(imread('./figures/schema.png'));
axis image
xticks([])
yticks([])
nexttile;
load(['.\data\c_cleaned\info_analysis_with_population_sizebtsp30angletancross.mat'], "avgbetastim_partitions_pop_size_ind", "avgbetachoice_partitions_pop_size_ind", "avgbetastimchoice_partitions_pop_size_ind");
betaweights_stim = squeeze(avgbetastim_partitions_pop_size_ind(subjects_list,:,:,:));
betaweights_choice = squeeze(avgbetachoice_partitions_pop_size_ind(subjects_list,:,:,:));
angles= squeeze(avgbetastimchoice_partitions_pop_size_ind(subjects_list,:,:,:));
x1 = betaweights_stim(:,:,1);
x2 = betaweights_stim(:,:,2);
y1 = betaweights_choice(:,:,1);
y2 = betaweights_choice(:,:,2);
betaweights_angle= atan2d(y2, y1) - atan2d(x2, x1);
desired_angle = 3;
chosen_pair = find((desired_angle-.1)<betaweights_angle & betaweights_angle<(desired_angle+.1),1);
[subject, repetition] = ind2sub(size(betaweights_angle),chosen_pair);
% subject = 5;
% repetition =59;
disp('Weights stim')
disp(betaweights_stim(subject,repetition))
disp('Weights choice')
disp(betaweights_choice(subject,repetition))
disp('Angel stim choice')
disp(betaweights_angle(subject,repetition))
disp(angles(subject,repetition))

slope_s_realdata = -betaweights_stim(subject,repetition,1)/betaweights_stim(subject,repetition,2);
slope_c_realdata = -betaweights_choice(subject,repetition,1)/betaweights_choice(subject,repetition,2);

plot(xvals,slope_s_realdata*xvals,'Color',[133 159 188]./255);hold on;
plot(xvals,slope_c_realdata*xvals,'Color',[225 149 113]./255)

xlabel('Neuron 1')
ylabel('Neuron 2')
%title('Example boundaries')
xticks([xvals(1) xvals(end)])
xticklabels({'L'; 'H'})
yl = ylim;
yticks([yl(1) yl(end)])
yticklabels({'L'; 'H'})
% patch([xvals(1) xvals(end) xvals(end) xvals(1) ], ...
%     [slope_s_realdata*xvals(1) slope_s_realdata*xvals(end) yl(2) yl(2)], 0, ...
%     'EdgeColor', 'none','FaceColor',[133 159 188]./255,'FaceAlpha',.15)
% patch([xvals(1) xvals(end) xvals(end) xvals(1) ], ...
%     [slope_c_realdata*xvals(1) slope_c_realdata*xvals(end) yl(2) yl(2)], 0,...
%     'EdgeColor', 'none','FaceColor',[225 149 113]./255,'FaceAlpha',.15)

text_real = ['\theta_{estim} = ' num2str(betaweights_angle(subject,repetition),3) '^{o}'];
text(-0.3,1.3,text_real)
leg = legend({'dec S boundary','dec C boundary'}, 'Location', 'Southwest');
leg.ItemTokenSize = [15,9];


load(['.\data\c_cleaned\info_analysis_with_population_sizebtsp30angletan.mat'], "avgbetastimchoice_partitions_pop_size_ind");
nexttile;
angles = avgbetastimchoice_partitions_pop_size_ind(subjects_list,:,:);
angles = acosd(cosd(angles));
angles(angles>90) = 180 - angles(angles>90);
histogram(angles)
xticks(-180:90:180)
xlabel('\theta_{estim} [deg]')
%title({'Boundary'; 'angle distribution'})
exportgraphics(hf, output_filename, 'Resolution',300)
end





