function [spike_train_sig_cells_peak_info_ind] = load_data_II_cells_spk_DFF_alltrials_toneD(spike_train_alltrials)
    %%
    GC_type  = {'H','M','C','F'}; time_scale_GC = 1; %%%time_scale_GC = 1: short time-scales; 
    load([pwd,'/Results/current_analysis/GCNetOutputs_10tp_earlyIIpeak']);      
    deconv_type = 1; type_data = 2; bin_size = 10; type_sig_cells = 'Sig_cells';    
    load([pwd,'/Results/current_analysis/metadata.mat']);
    load([pwd,'/Results/current_analysis/Neurons_class_GC_sliding_window_',char(type_sig_cells),'_',num2str(bin_size),'tp.mat'],'II_CLASS_PEAK_value_mice'); 
    load([pwd,'/Results/current_analysis/significant_neurons_sliding_window_',num2str(bin_size),'tp.mat'],'significant_neurons');    
    spike_train_sig_cells_peak_info_ind = cell(1,length(II_CLASS_PEAK_value_mice));
    for ind   = 1:length(II_CLASS_PEAK_value_mice)        
        if  length(significant_neurons{ind}) > 20  
            GC_neurons  = GCNets.(GC_type{1}){time_scale_GC,ind}.cellids;                        
            c_vec        = metadata(ind,type_data).c_vec;
            n_trials    = length(find(c_vec'~=0));
            sigNeuronTime_pks = II_CLASS_PEAK_value_mice{ind}(GC_neurons); 
            sigNeurons = GC_neurons;
            spike_train_sig_cells_peak_info = zeros(n_trials, length(sigNeurons));
            for n = 1:length(sigNeurons)
                spike_train_sig_cells_peak_info(:,n)  = ...
                    mean(spike_train_alltrials(ind,type_data,deconv_type).spike_train_alltrials(sigNeuronTime_pks(n):sigNeuronTime_pks(n)+bin_size,c_vec'~=0,sigNeurons(n)));                    
            end
            spike_train_sig_cells_peak_info_ind{ind} = spike_train_sig_cells_peak_info;        
        end
    end
end