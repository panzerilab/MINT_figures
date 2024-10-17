load()
pop_size_ALL = 1:19;
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
