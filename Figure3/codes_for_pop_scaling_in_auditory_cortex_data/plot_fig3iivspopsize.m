load('Results\current_analysis\info_analysis_with_population_sizebtspqe.mat', ...
    "avgSI_popqe_partitions_pop_size_ind", ...
    "avgCI_popqe_partitions_pop_size_ind",...
    "avgII_popqe_partitions_pop_size_ind", ...
    "avgSI_popnaive_partitions_pop_size_ind", ...
    "avgCI_popnaive_partitions_pop_size_ind",...
    "avgII_popnaive_partitions_pop_size_ind", ...
    "sample_combinations_pop_size", ...
    "SI_CI_info_ind", ...
    "SI_CI_info_indqe", ...
    "pop_size_ALL", ...
    "Nsessions");
%% Load
Nsessionsclean = Nsessions(Nsessions>0);
% Nsessionsclean = Nsessionsclean(Nsessionsclean~=14);
% Nsessionsclean = Nsessionsclean(Nsessionsclean~=24);
% Nsessionsclean =[5, 6];
sample_combinations_pop_size =[20 repmat(100, 1, 17) 20 1];
ci_naive = sum(avgCI_popnaive_partitions_pop_size_ind(Nsessionsclean,:,:),3);
si_naive = sum(avgSI_popnaive_partitions_pop_size_ind(Nsessionsclean,:,:),3);
ii_naive = sum(avgII_popnaive_partitions_pop_size_ind(Nsessionsclean,:,:),3);

sizes = 1:19;

ci_naive = ci_naive(:,sizes)./sample_combinations_pop_size(sizes);
si_naive = si_naive(:,sizes)./sample_combinations_pop_size(sizes);
ii_naive = ii_naive(:,sizes)./sample_combinations_pop_size(sizes);

ci_qe = sum(avgCI_popqe_partitions_pop_size_ind(Nsessionsclean,:,:),3);
si_qe = sum(avgSI_popqe_partitions_pop_size_ind(Nsessionsclean,:,:),3);
% ii_qe = sum(avgII_popqe_partitions_pop_size_ind(Nsessionsclean,:,:,3),3);
ii_qe = sum(avgII_popbtsp_partitions_pop_size_ind(Nsessionsclean,:,:),3);

ci_qe = ci_qe(:,sizes)./sample_combinations_pop_size(sizes);
si_qe = si_qe(:,sizes)./sample_combinations_pop_size(sizes);
ii_qe = ii_qe(:,sizes)./sample_combinations_pop_size(sizes);

%% Naive
% ci_naive(:,2:end-1) = ci_naive(:,2:end-1)/100;
% si_naive(:,2:end-1) = si_naive(:,2:end-1)/100;
% ii_naive(:,2:end-1) = ii_naive(:,2:end-1)/100;
% 
% ci_naive(:,1)=ci_naive(:,1)/20;
% si_naive(:,1)=si_naive(:,1)/20;
% ii_naive(:,1)=ii_naive(:,1)/20;
% 
% ci_naive(:,end)=ci_naive(:,end)/20;
% si_naive(:,end)=si_naive(:,end)/20;
% ii_naive(:,end)=ii_naive(:,end)/20;
% 
% 
% %% QE
% 
% 
% ci_qe(:,2:end-1) = ci_qe(:,2:end-1)/100;
% si_qe(:,2:end-1) = si_qe(:,2:end-1)/100;
% ii_qe(:,2:end-1) = ii_qe(:,2:end-1)/100;
% 
% ci_qe(:,1)=ci_qe(:,1)/20;
% si_qe(:,1)=si_qe(:,1)/20;
% ii_qe(:,1)=ii_qe(:,1)/20;
% 
% ci_qe(:,end)=ci_qe(:,end)/20;
% si_qe(:,end)=si_qe(:,end)/20;
% ii_qe(:,end)=ii_qe(:,end)/20;


%% Fit Ince
pop_size_ALL = pop_size_ALL(sizes);
% polynomial curve
% fit_function = 'Ince2013_model_2params';
fit_function = 'Ince2013_model_1param';
% fit_function = 'logistic';
[yFit_meanSInaive, beta_SInaive, R2_SInaive, RMSE_SInaive, asymp_SInaive] = fit_info_curves(pop_size_ALL, mean(si_naive,1), fit_function, 'SI', SI_CI_info_ind);
[yFit_meanCInaive, beta_CInaive, R2_CInaive, RMSE_CInaive, asymp_CInaive] = fit_info_curves(pop_size_ALL, mean(ci_naive,1), fit_function, 'CI', SI_CI_info_ind);
[yFit_meanIInaive, beta_IInaive, R2_IInaive, RMSE_IInaive, asymp_IInaive] = fit_info_curves(pop_size_ALL, mean(ii_naive,1), fit_function, 'II', SI_CI_info_ind);

[yFit_meanSIqe, beta_SIqe, R2_SIqe, RMSE_SIqe, asymp_SIqe] = fit_info_curves(pop_size_ALL, mean(si_qe,1), fit_function, 'SI', SI_CI_info_indqe);
[yFit_meanCIqe, beta_CIqe, R2_CIqe, RMSE_CIqe, asymp_CIqe] = fit_info_curves(pop_size_ALL, mean(ci_qe,1), fit_function, 'CI', SI_CI_info_indqe);
[yFit_meanIIqe, beta_IIqe, R2_IIqe, RMSE_IIqe, asymp_IIqe] = fit_info_curves(pop_size_ALL, mean(ii_qe,1), fit_function, 'II', SI_CI_info_indqe);


%% Plot
figure

ax1 = subplot(1,2,1);

scatter(pop_size_ALL, mean(si_naive,1),'filled','b');hold on ;
scatter(pop_size_ALL, mean(ci_naive,1),'filled','g');
scatter(pop_size_ALL, mean(ii_naive,1),'filled','r');

errorbar(pop_size_ALL, mean(si_naive,1),std(si_naive,[],1)/sqrt(length(Nsessionsclean)),'b','LineStyle', 'none')
errorbar(pop_size_ALL, mean(ci_naive,1),std(ci_naive,[],1)/sqrt(length(Nsessionsclean)),'g','LineStyle', 'none')
errorbar(pop_size_ALL, mean(ii_naive,1),std(ii_naive,[],1)/sqrt(length(Nsessionsclean)),'r','LineStyle', 'none')

title('Naive')


% fitSInaive = plot(pop_size_ALL, yFit_meanSInaive, 'b', 'LineWidth', 1.5); hold on;
% fitCInaive = plot(pop_size_ALL, yFit_meanCInaive, 'g', 'LineWidth', 1.5); hold on;
% fitIInaive = plot(pop_size_ALL, yFit_meanIInaive, 'r', 'LineWidth', 1.5); hold on;
behav_perf = yline(nanmean(SI_CI_info_ind),'k-.','LineWidth',2);
% yline(asymp_SInaive,'b--');
% yline(asymp_CInaive,'g--');
% yline(asymp_IInaive,'r--');
xlim([1,19])
% ylim([0,.55])

legend('CI', 'SI', 'II','Location','southeast')

ax2=subplot(1,2,2);

scatter(pop_size_ALL, mean(si_qe),'filled','b');hold on ;
scatter(pop_size_ALL, mean(ci_qe),'filled','g');
scatter(pop_size_ALL, mean(ii_qe),'filled','r');

errorbar(pop_size_ALL, mean(si_qe,1),std(si_qe,[],1)/sqrt(length(Nsessionsclean)),'b','LineStyle', 'none')
errorbar(pop_size_ALL, mean(ci_qe,1),std(ci_qe,[],1)/sqrt(length(Nsessionsclean)),'g','LineStyle', 'none')
errorbar(pop_size_ALL, mean(ii_qe,1),std(ii_qe,[],1)/sqrt(length(Nsessionsclean)),'r','LineStyle', 'none')

title('BTSP')
xlim([1,19])
% ylim([0,.55])
linkaxes([ax1 ax2], 'xy')


% fitSIqe = plot(pop_size_ALL, yFit_meanSIqe, 'b', 'LineWidth', 1.5); hold on;
% fitCIqe = plot(pop_size_ALL, yFit_meanCIqe, 'g', 'LineWidth', 1.5); hold on;
% fitIIqe = plot(pop_size_ALL, yFit_meanIIqe, 'r', 'LineWidth', 1.5); hold on;
behav_perfqe = yline(nanmean(SI_CI_info_indqe),'k-.','LineWidth',2);
% yline(asymp_SIqe,'b--');
% yline(asymp_CIqe,'g--');
% yline(asymp_IIqe,'r--');

legend('CI',  'SI',  'II','Location','southeast')