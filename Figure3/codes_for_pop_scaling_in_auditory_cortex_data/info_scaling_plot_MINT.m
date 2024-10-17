function info_scaling_plot_MINT()
%%
% Plot the results
% load(['Results\current_analysis\info_analysis_with_population_sizeMatias.mat'],"avgSI_pop_partitions_pop_size_ind","avgCI_pop_partitions_pop_size_ind",...
%     "avgII_pop_partitions_pop_size_ind","Nsessions","sample_combinations_pop_size","SI_CI_info_ind");
load('Results\current_analysis\info_analysis_with_population_sizebtspqe.mat', ...
    "avgSI_popqe_partitions_pop_size_ind", ...
    "avgCI_popqe_partitions_pop_size_ind",...
    "avgII_popqe_partitions_pop_size_ind", ...
    "avgSI_popnaive_partitions_pop_size_ind", ...
    "avgCI_popnaive_partitions_pop_size_ind",...
    "avgII_popnaive_partitions_pop_size_ind", ...
    "sample_combinations_pop_size", ...
    "SI_CI_info_ind", ...
    "pop_size_ALL", ...
    "Nsessions");
% pop_size_ALL = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];
% Nsessions = 1:34;
pop_size_ALL = [1:19];
avgSI_pop_partitions_pop_size_ind = avgSI_popqe_partitions_pop_size_ind(:,1:19,:);
avgCI_pop_partitions_pop_size_ind = avgCI_popqe_partitions_pop_size_ind(:,1:19,:);
avgII_pop_partitions_pop_size_ind = avgII_popqe_partitions_pop_size_ind(:,1:19,:);


% Pool the data by taking the average across repetitions for each session and then find the mean and SEM across sessions
[avgSI_pop_partitions_pop_size_ind_across_subj_pairs] = pool_infodata_to_plot_as_info_vs_popsize(avgSI_pop_partitions_pop_size_ind,Nsessions,sample_combinations_pop_size);
[avgCI_pop_partitions_pop_size_ind_across_subj_pairs] = pool_infodata_to_plot_as_info_vs_popsize(avgCI_pop_partitions_pop_size_ind,Nsessions,sample_combinations_pop_size);
[avgII_pop_partitions_pop_size_ind_across_subj_pairs] = pool_infodata_to_plot_as_info_vs_popsize(avgII_pop_partitions_pop_size_ind,Nsessions,sample_combinations_pop_size); %[nonzeros(IIpersubject(:,1)) nonzeros(IIpersubject(:,2)) nonzeros(IIpersubject(:,3)) nonzeros(IIpersubject(:,4))] ;%

meanSI = mean(avgSI_pop_partitions_pop_size_ind_across_subj_pairs,1);
stdSI = std(avgSI_pop_partitions_pop_size_ind_across_subj_pairs)/sqrt(length(Nsessions)-1);
meanCI = mean(avgCI_pop_partitions_pop_size_ind_across_subj_pairs,1);
stdCI = std(avgCI_pop_partitions_pop_size_ind_across_subj_pairs)/sqrt(length(Nsessions)-1);
meanII = mean(avgII_pop_partitions_pop_size_ind_across_subj_pairs,1);
stdII = std(avgII_pop_partitions_pop_size_ind_across_subj_pairs)/sqrt(length(Nsessions)-1);


meanSI(1) = meanSI(1)*5;
stdSI(1) = stdSI(1)*5;
meanCI(1) = meanCI(1)*5;
stdCI(1) = stdCI(1)*5;
meanII(1) = meanII(1)*5;
stdII(1) = stdII(1)*5;

meanSI(19) = meanSI(19)*5;
stdSI(19) = stdSI(19)*5;
meanCI(19) = meanCI(19)*5;
stdCI(19) = stdCI(19)*5;
meanII(19) = meanII(19)*5;
stdII(19) = stdII(19)*5;

% meanSI(20) = meanSI(20)*100;
% stdSI(20) = stdSI(20)*100;
% meanCI(20) = meanCI(20)*100;
% stdCI(20) = stdCI(20)*100;
% meanII(20) = meanII(20)*100;
% stdII(20) = stdII(20)*100;
% Fit the information curves with a logistic growth model or with a
% polynomial curve
fit_function = 'Ince2013_model_2params';
% fit_function = 'Ince2013_model_1param';
% fit_function = 'logistic';
info_scaling = 'SI';
[yFit_meanSI, beta_SI, R2_SI, RMSE_SI, asymp_SI] = fit_info_curves(pop_size_ALL, meanSI, fit_function, info_scaling, SI_CI_info_ind);
info_scaling = 'CI';
[yFit_meanCI, beta_CI, R2_CI, RMSE_CI, asymp_CI] = fit_info_curves(pop_size_ALL, meanCI, fit_function, info_scaling, SI_CI_info_ind);
info_scaling = 'II';
[yFit_meanII, beta_II, R2_II, RMSE_II, asymp_II] = fit_info_curves(pop_size_ALL, meanII, fit_function, info_scaling, SI_CI_info_ind);

beta_SI
beta_CI
beta_II
% Plot the information scaling curves as a function of the population size
figure('Visible','on'); set(gcf,'Position',[00, 00, 300, 300]); hold on;
% title('Information scaling'); hold on;
title(fit_function); hold on;
transparency_areas = 0.1;
errorbar(pop_size_ALL, meanSI, stdSI,'b','LineStyle', 'none')
errorbar(pop_size_ALL, meanCI, stdCI,'g','LineStyle', 'none')
errorbar(pop_size_ALL, meanII, stdII,'r','LineStyle', 'none')
% fill([pop_size_ALL fliplr(pop_size_ALL)],[meanSI+stdSI fliplr(meanSI-stdSI)], 'b','FaceAlpha',transparency_areas,'EdgeColor','none');hold on;
% fill([pop_size_ALL fliplr(pop_size_ALL)],[meanCI+stdCI fliplr(meanCI-stdCI)], 'g','FaceAlpha',transparency_areas,'EdgeColor','none');hold on;
% fill([pop_size_ALL fliplr(pop_size_ALL)],[meanII+stdII fliplr(meanII-stdII)], 'r','FaceAlpha',transparency_areas,'EdgeColor','none');hold on;
scatter(pop_size_ALL,meanSI,20,'filled','b'); hold on;
scatter(pop_size_ALL,meanCI,20,'filled','g'); hold on;
scatter(pop_size_ALL,meanII,20,'filled','r'); hold on;
yline(asymp_SI,'b--');
yline(asymp_CI,'g--');
yline(asymp_II,'r--');
behav_perf = yline(nanmean(SI_CI_info_ind),'k-.','LineWidth',2);

% Plot the fitted curve
fitSI = plot(pop_size_ALL, yFit_meanSI, 'b', 'LineWidth', 1.5); hold on;
fitCI = plot(pop_size_ALL, yFit_meanCI, 'g', 'LineWidth', 1.5); hold on;
fitII = plot(pop_size_ALL, yFit_meanII, 'r', 'LineWidth', 1.5); hold on;
% h = legend([fitSI,fitCI,fitII,behav_perf],{'I(S,$\hat{S}$)','I(C,$\hat{C}$)','II(S,[$\hat{S},\hat{C}$],C)','I(S,C)'},...
%     'Location','southeast','Interpreter', 'latex'); 
h = legend([fitSI,fitCI,fitII,behav_perf],{'I(S,$\hat{S}$)','I(C,$\hat{C}$)','II(S,[$\hat{S},\hat{C}$],C)','I(S,C)'},...
    'Location','southeast','Interpreter', 'latex'); 
ylabel('bits');xlabel('population size'); 
xlim([1 19]); ylim([0 0.55])
filename = ['Figures/Population_based_SI_CI_II_fixedextrnaive',char(fit_function),'.pdf']; 
exportgraphics(gcf,filename,'ContentType','vector');    
end