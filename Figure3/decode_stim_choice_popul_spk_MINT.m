function decode_stim_choice_popul_spk_MINT(spike_train_sig_cells_peak_info_ind, spike_train_alltrials)
%%
% Load the stimuli and choices across trials
load('Results/current_analysis/metadata.mat');
deconv_type = 1; type_data = 2; bin_size = 10;     
load([pwd,'/Results/current_analysis/significant_neurons_sliding_window_',num2str(bin_size),'tp.mat'],'significant_neurons');    

pop_size_ALL = [5:5:15];
sample_combinations_pop_size = 100;
number_of_sessions = 34;
% Find NaN trials in spike train data to remove from the entire analysis
nan_trials_in_spike_trains = cell(1,34);
for ind = 1:number_of_sessions
    if  length(significant_neurons{ind}) > 20  
        c_vec            = metadata(ind,type_data).c_vec;
        spike_train = spike_train_alltrials(ind,type_data,deconv_type).spike_train_alltrials(:,c_vec' ~= 0,:);
        nan_trials_in_spike_trains{ind} = find((sum(sum(isnan(spike_train),1),3)));
    end
end

%% Collect all trials
stimulus_collected = NaN(1,475);
choice_collected = NaN(1,475);
spikes_collected = NaN(1,475);

for ind = 1:number_of_sessions
    if  length(significant_neurons{ind}) > 20  
        stimulus_vector  = metadata(ind,type_data).stimulus_vector;
        choice_vector   = metadata(ind,type_data).choice_vector;
        stimulus_vector(nan_trials_in_spike_trains{ind}) = [];
        choice_vector(nan_trials_in_spike_trains{ind}) = [];
    end
end


xtrp = 10;
% Parameters for computing the Shannon mutual information
xMI_opts.verbose = false;
xMI_opts.method = 'dr';
xMI_opts.bias = 'naive';
xMI_opts.btsp = 0; 
xMI_opts.bin_methodX = 'none'; 
xMI_opts.bin_methodY = 'none';
xMI_opts.n_binsX = 2;
xMI_opts.n_binsY = 2;

MIqe_opts = MI_opts;
MIqe_opts.bias = 'qe';
MIqe_opts.xtrp = xtrp;

% Parameters for computing the Intersection information
II_opts.bin_methodS = 'none';                             
II_opts.bin_methodR = 'none';                    
II_opts.bin_methodC = 'none';                    
II_opts.n_binsS = 2;                        
II_opts.n_binsR = 2;
II_opts.n_binsC = 2;  
II_opts.bias = 'qe';
II_opts.xtrp = xtrp;

% Parameters for computing the SVM
opts = []; 
opts.optimize_params = false; 
opts.parallel_optim = true;
opts.hp_C  = 1;
opts.cv_type = "KFold"; 
nfolds = 2;
opts.K = nfolds; % opts.cv_type = 'LeaveOneOut';
opts.svm_family = 'RBF'; % or 'RBF'
opts.libsvm = false;

% Compute the information about stimulus and choice variables across sessions
SI_CI_info_ind = NaN(1,number_of_sessions);
for ind = 1:number_of_sessions
    if  length(significant_neurons{ind}) > 20  
        stimulus_vector  = metadata(ind,type_data).stimulus_vector;
        choice_vector   = metadata(ind,type_data).choice_vector;
        stimulus_vector(nan_trials_in_spike_trains{ind}) = [];
        choice_vector(nan_trials_in_spike_trains{ind}) = [];
        [SI_CI_info] = information(stimulus_vector, choice_vector, MI_opts, {'I'});
        SI_CI_info_ind(ind) = SI_CI_info{1};
    end
end

% Compute the mutual and intersection information at the population level
% based on the decoded stimuli and choices obatined through the SVM analysis
tic;
Nsessions = 0;


avgSI_popqe_partitions_pop_size_ind = zeros(number_of_sessions,length(pop_size_ALL), sample_combinations_pop_size); 
avgCI_popqe_partitions_pop_size_ind = zeros(number_of_sessions,length(pop_size_ALL), sample_combinations_pop_size); 
avgII_popqe_partitions_pop_size_ind = zeros(number_of_sessions,length(pop_size_ALL), sample_combinations_pop_size); 
avgSI_popnaive_partitions_pop_size_ind = zeros(number_of_sessions,length(pop_size_ALL), sample_combinations_pop_size); 
avgCI_popnaive_partitions_pop_size_ind = zeros(number_of_sessions,length(pop_size_ALL), sample_combinations_pop_size); 
avgII_popnaive_partitions_pop_size_ind = zeros(number_of_sessions,length(pop_size_ALL), sample_combinations_pop_size); 
for ind = 14:number_of_sessions
    if  length(significant_neurons{ind}) > 20  
        ind
        stimulus_vector  = metadata(ind,type_data).stimulus_vector;
        choice_vector   = metadata(ind,type_data).choice_vector;
        stimulus_vector(nan_trials_in_spike_trains{ind}) = [];
        choice_vector(nan_trials_in_spike_trains{ind}) = [];        
        cvPartition = cvpartition(stimulus_vector, 'KFold', nfolds); % Should these indices be the same for the analyses of S,C,II?
        avgSI_popqe_partitions_pop_size = zeros(length(pop_size_ALL), sample_combinations_pop_size);
        avgCI_popqe_partitions_pop_size = zeros(length(pop_size_ALL), sample_combinations_pop_size);
        avgII_popqe_partitions_pop_size = zeros(length(pop_size_ALL), sample_combinations_pop_size);
        avgSI_popnaive_partitions_pop_size = zeros(length(pop_size_ALL), sample_combinations_pop_size);
        avgCI_popnaive_partitions_pop_size = zeros(length(pop_size_ALL), sample_combinations_pop_size);
        avgII_popnaive_partitions_pop_size = zeros(length(pop_size_ALL), sample_combinations_pop_size);

        for pop_size_idx = 1:length(pop_size_ALL)
            fprintf('Pop. size %d is sent at time %s\n', pop_size_ALL(pop_size_idx), datestr(now,'HH:MM:SS.FFF'))
            allset_elements = nchoosek(1:20,pop_size_ALL(pop_size_idx));
            idx_selected_combinations = randperm(length(allset_elements));
            if (length(idx_selected_combinations)<sample_combinations_pop_size) && pop_size_ALL(pop_size_idx)<20
                select_allset_elements = allset_elements(idx_selected_combinations,:);
            elseif pop_size_ALL(pop_size_idx)==20
                select_allset_elements = allset_elements(idx_selected_combinations);
            else
                select_allset_elements = allset_elements(idx_selected_combinations(1:sample_combinations_pop_size),:);
            end
            avgSI_popqe_partitions = zeros(1,sample_combinations_pop_size);
            avgCI_popqe_partitions = zeros(1,sample_combinations_pop_size);
            avgII_popqe_partitions = zeros(1,sample_combinations_pop_size);
            avgSI_popnaive_partitions = zeros(1,sample_combinations_pop_size);
            avgCI_popnaive_partitions = zeros(1,sample_combinations_pop_size);
            avgII_popnaive_partitions = zeros(1,sample_combinations_pop_size);

            number_of_partitions = cvPartition.NumTestSets;
            pop_size = pop_size_ALL(pop_size_idx);
            % test_indexes = cvPartition.test;
            for sample = 1:min(sample_combinations_pop_size, length(idx_selected_combinations))
                data = spike_train_sig_cells_peak_info_ind{ind}(:,select_allset_elements(sample, 1:pop_size));
                data(nan_trials_in_spike_trains{ind},:) = [];
                SI_popqe_partitions = zeros(1,number_of_partitions);
                CI_popqe_partitions = zeros(1,number_of_partitions);
                II_popqe_partitions = zeros(1,number_of_partitions);
                SI_popnaive_partitions = zeros(1,number_of_partitions);
                CI_popnaive_partitions = zeros(1,number_of_partitions);
                II_popnaive_partitions = zeros(1,number_of_partitions);
                for l = 1:number_of_partitions
                    test_idxs = find(cvPartition.test(l));

                    stim_choice = map_Nd_array_to_1d([stimulus_vector; choice_vector]); % This function will pool stimuli and choices
                    % [stimchoice_vector_decoded, ~, ~, ~, ~] = svm_pipeline(double(data),stim_choice',test_idxs, opts); 
                    [stimulus_vector_decoded, ~, ~, ~, ~] = svm_pipeline(double(data),stimulus_vector,test_idxs, opts); 
                    [choice_vector_decoded, ~, ~, ~, ~] = svm_pipeline(double(data),choice_vector, test_idxs, opts);
                    stimchoice_vector_decoded = map_Nd_array_to_1d([stimulus_vector_decoded'; choice_vector_decoded']);

                    % Decode the stimuli and compute the information of the
                    % decoded stimuli compared to the actual stimuli
                    % [stimulus_vector_decoded, ~, ~, ~, ~] = svm_pipeline(double(data),stimulus_vector,test_idxs, opts); 
                    [SI_popnaive] = information(stimchoice_vector_decoded, stimulus_vector(test_idxs), MI_opts, {'I'});
                    [SI_popqe] = information(stimchoice_vector_decoded, stimulus_vector(test_idxs), MIqe_opts, {'I'});
%                     [SI_pop] = information(stimchoice_vector_decoded', stimulus_vector(test_idxs)', MI_opts, {'I'});
                    SI_popqe_partitions(l) = SI_popqe{1};
                    SI_popnaive_partitions(l) = SI_popnaive{1};
    
                    % Decode the choices and compute the information of the
                    % decoded choices compared to the actual choices
                    % [choice_vector_decoded, ~, ~, ~, ~] = svm_pipeline(double(data),choice_vector, test_idxs, opts);
                    [CI_popnaive] = information(stimchoice_vector_decoded, choice_vector(test_idxs), MI_opts, {'I'});
                    [CI_popqe] = information(stimchoice_vector_decoded, choice_vector(test_idxs), MIqe_opts, {'I'});
                    
%                     [CI_pop] = information(stimchoice_vector_decoded', choice_vector(test_idxs)', MI_opts, {'I'});
                    CI_popqe_partitions(l) = CI_popqe{1};
                    CI_popnaive_partitions(l) = CI_popnaive{1};
    
                    % Decode BOTH stimuli and choices and compute the
                    % intersection information of the decoded stimuli AND
                    % with the actual stimuli and choices
%                     stim_choice = map_Nd_array_to_1d([stimulus_vector'; choice_vector']); % This function will pool stimuli and choices
%                     [stimchoice_vector_decoded, ~, ~, ~, ~] = svm_pipeline(double(data),stim_choice', test_idxs, opts);               
                    [II_pop_qe, II_pop_naive] = II(stimulus_vector(test_idxs), stimchoice_vector_decoded, choice_vector(test_idxs), II_opts);

                    % % % % % % % % % % II_pop = II_out.biased.value;
                    II_popqe_partitions(l) = II_pop_qe;                               
                    II_popnaive_partitions(l) = II_pop_naive;
                end
                avgSI_popqe_partitions(sample) = mean(SI_popqe_partitions);
                avgCI_popqe_partitions(sample) = mean(CI_popqe_partitions);
                avgII_popqe_partitions(sample) = mean(II_popqe_partitions);
                avgSI_popnaive_partitions(sample) = mean(SI_popnaive_partitions);
                avgCI_popnaive_partitions(sample) = mean(CI_popnaive_partitions);
                avgII_popnaive_partitions(sample) = mean(II_popnaive_partitions);
            end

            avgSI_popqe_partitions_pop_size(pop_size_idx,:) = avgSI_popqe_partitions;
            avgCI_popqe_partitions_pop_size(pop_size_idx,:) = avgCI_popqe_partitions;
            avgII_popqe_partitions_pop_size(pop_size_idx,:) = avgII_popqe_partitions;
            avgSI_popnaive_partitions_pop_size(pop_size_idx,:) = avgSI_popnaive_partitions;
            avgCI_popnaive_partitions_pop_size(pop_size_idx,:) = avgCI_popnaive_partitions;
            avgII_popnaive_partitions_pop_size(pop_size_idx,:) = avgII_popnaive_partitions;
        end
        Nsessions = [Nsessions, ind];
        avgSI_popqe_partitions_pop_size_ind(ind,:,:) = avgSI_popqe_partitions_pop_size;
        avgCI_popqe_partitions_pop_size_ind(ind,:,:) = avgCI_popqe_partitions_pop_size;
        avgII_popqe_partitions_pop_size_ind(ind,:,:) = avgII_popqe_partitions_pop_size;
        avgSI_popnaive_partitions_pop_size_ind(ind,:,:) = avgSI_popnaive_partitions_pop_size;
        avgCI_popnaive_partitions_pop_size_ind(ind,:,:) = avgCI_popnaive_partitions_pop_size;
        avgII_popnaive_partitions_pop_size_ind(ind,:,:) = avgII_popnaive_partitions_pop_size;
    end
end
toc;
save('Results\current_analysis\info_analysis_with_population_sizekk.mat', ...
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
end