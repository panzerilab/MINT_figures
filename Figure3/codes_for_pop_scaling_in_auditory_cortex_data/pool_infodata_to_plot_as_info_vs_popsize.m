function [avgSI_pop_partitions_pop_size_ind_across_subj_pairs_output] = pool_infodata_to_plot_as_info_vs_popsize(avgSI_pop_partitions_pop_size_ind,Nsessions,sample_combinations_pop_size)
%%
    idx = 0;
    avgSI_pop_partitions_pop_size_ind_1subj = cell(1,length(Nsessions)-1);
    for sessionIdx = Nsessions(2:end)
        idx = 1 + idx;
        avgSI_pop_partitions_pop_size_ind_1subj{idx} = avgSI_pop_partitions_pop_size_ind(sessionIdx,:,:);
    end
    
    avgSI_pop_partitions_pop_size_ind_across_subj = cell2mat(avgSI_pop_partitions_pop_size_ind_1subj');
    
%     avgSI_pop_partitions_pop_size_ind_across_subj_pairs = [];
%     for set = 1:sample_combinations_pop_size
%         avgSI_pop_partitions_pop_size_ind_across_subj_pairs =  vertcat(avgSI_pop_partitions_pop_size_ind_across_subj(:,:,set), avgSI_pop_partitions_pop_size_ind_across_subj_pairs);
%     end
%     avgSI_pop_partitions_pop_size_ind_across_subj_pairs_output = avgSI_pop_partitions_pop_size_ind_across_subj_pairs;
    avgSI_pop_partitions_pop_size_ind_across_subj_pairs_output = mean(avgSI_pop_partitions_pop_size_ind_across_subj,3);

    
end