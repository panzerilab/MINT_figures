function extract_trials_stim_choice
%% Extract trials stim and choice
    processed_original_data                     = dir(['Data/data_for_information_analysis/*.mat']);  %%% Consider all experiments and mice
    type_data_name                              = {'corrected','smoothed'};
    metadata                                    = struct;
    for ind                                     = 1:length(processed_original_data)  %%% Generalize with a parfor here (Load issue not save)
        ind
        for type_data = 1:length(type_data_name)
            load(['Data/data_for_information_analysis/', char(processed_original_data(ind).name)]);       %%% Load file with Fluorescence responses, stimuli and choices
            metadata(ind,type_data).CellLoc                             = cell_location_mice(brightCells == 1,:);
            metadata(ind,type_data).name_exp_mice                       = name_exp_mice;
            DFF_data_Beh                                                = {DFF_corrected_Beh, DFF_smoothed_Beh};
            DeltaF_trials_stim                                          = DFF_data_Beh{type_data}(:, :, brightCells == 1);     %%% Consider only responsive cells at a given trial
            metadata(ind,type_data).nan_trials                          = sum(sum(isnan(DeltaF_trials_stim)),3);

            metadata(ind,type_data).hits_vec                            = hit_trials(31,:);   %%% Hits = c1
            metadata(ind,type_data).FA_vec                              = FA_trials(31,:);    %%% FA = c1

            metadata(ind,type_data).early_vec                           = early_trials(31,:);    %%% FA = c1
            metadata(ind,type_data).early_hit_vec                       = early_hit_trials(31,:);    %%% FA = c1
            metadata(ind,type_data).early_FA_vec                        = early_FA_trials(31,:);    %%% FA = c1
            
            miss_vec                                                    = miss_trials(31,:);  %%% Miss = c2
            miss_vec(miss_vec == 1)                                     = 2;                                %%% Re-label miss = c2
            metadata(ind,type_data).miss_vec                            = miss_vec;
            
            CR_vec                                                      = CR_trials(31,:);    %%%
            CR_vec(CR_vec == 1)                                         = 2;                            %%% Re-label CR = c2
            metadata(ind,type_data).CR_vec                              = CR_vec;

            c_vec                                                       = metadata(ind,type_data).hits_vec + metadata(ind,type_data).FA_vec + metadata(ind,type_data).miss_vec + metadata(ind,type_data).CR_vec; %%% choice vector (without the early trials) for the information toolbox
            metadata(ind,type_data).c_vec                               = c_vec;
            metadata(ind,type_data).choice_vector                       = c_vec(c_vec'~=0);              %%% Select those neurons that have made a choice (H,M,FA,CR) (no early trials)
            
            s_vec                                                       = StimFrequency(31,:); 
            s_vec(s_vec<11000)                                          = 1;                            %%% Re-label low freqs  = s1
            s_vec(s_vec>11000)                                          = 2;                            %%% Re-label high freqs = s2
            metadata(ind,type_data).stimulus_vector                     = s_vec(c_vec'~=0);              %%% Select those stimuli that correspond to the choice (H,M,FA,CR) (no early trials)
            metadata(ind,type_data).considered_trials                   = sum(c_vec'~=0);                %%% Count number of trials with a choice (H,M,FA,CR)
            metadata(ind, type_data).FirstResponse                      = FirstResponse;
        end
    end
    save('Results/metadata.mat','metadata');
end