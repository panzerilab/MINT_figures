clc, clear all
% Load the stimuli and choices across trials
load('Data_A1/stim_choice.mat');
deconv_type = 1; type_data = 2; bin_size = 10;
load(['Data_A1/significant_neurons_sliding_window_',num2str(bin_size),'tp.mat'],'significant_neurons');
load('Data_A1/DFF_spktrain_alltrials.mat','spike_train_alltrials');
load('Data_A1/preprocessed_spk_input_data.mat','spike_train_sig_cells_peak_info_ind');

nSessions = 34;
nShuff = 2; 
nan_trials_in_spike_trains = cell(1,nSessions);
for ind = 1:nSessions
    if  length(significant_neurons{ind}) > 20
        C            = stim_choice(ind,type_data).c_vec;
        spike_train = spike_train_alltrials(ind,type_data,deconv_type).spike_train_alltrials(:,C' ~= 0,:);
        nan_trials_in_spike_trains{ind} = find((sum(sum(isnan(spike_train),1),3)));
    end
end

% -------------------------------------------------------------------------
% Define optional Inputs for MI and PID
% ------------------------------------------------------------------------

opts_MI.bias = 'shuffSub';
opts_MI.shuff = 30;
opts_MI.bin_method = {'none', 'none'};
opts_MI.supressWarnings = true;

opts_PID.bias = 'shuffSub';
opts_PID.shuff = 30;
opts_PID.bin_method = {'none', 'none'};
opts_PID.supressWarnings = true;

SVM_opts.cv = {'KFold', 2};
SVM_opts.libsvm = false;
optimization.optim_cv = {'KFold', 2};
optimization.parallel = true;
SVM_opts.optim_opts = optimization;

reqOutputs_MI = {'I(A;B)', 'Ilin(A;B)', 'Iss(A)', 'Ici(A;B)', 'Icd(A;B)'};
reqOutputs_PID = {'Joint','PID_atoms'};
reqOutputs_SVM = {'labels', 'testIdx'};
reqOutputs_popMI = {'I(A;B)'};


% Info breakdown opts
info_bdw_terms = {'Joint','ILIN','ISS','ICI','ICD'};
bdw_bins = 5;% number of bins used to discretize spike counts for information breakdown calculations
field_names = {};
for ind = 1:nSessions
    if  length(significant_neurons{ind}) > 20
        field_name = sprintf('session_%d', ind);
        field_names = [field_names, {field_name}];
        disp(field_name);
        S = stim_choice(ind,type_data).stimulus_vector;
        S(nan_trials_in_spike_trains{ind}) = [];
        R_all = spike_train_sig_cells_peak_info_ind{ind};
        R_all(nan_trials_in_spike_trains{ind},:) = [];
        R = ((R_all > 0) + (R_all > 1))';
        [nNeurons, nTrials] = size(R);
        nPairs = nNeurons*(nNeurons-1)/2;
        % -------------------------------------------------------------------------
        % Step 2: Compute PID and Information Breakdown
        % -------------------------------------------------------------------------    
        pairlist = nchoosek(1:nNeurons,2);
        MI_v = cell(length(pairlist), 5);
        PID_v = cell(length(pairlist), 5);       
        for pairi = 1:length(pairlist)
            cell1 = pairlist(pairi,1);
            cell2 = pairlist(pairi,2);
            resp1 = R(cell1,:);
            resp2 = R(cell2,:);
            resp1(resp1>bdw_bins -1) = bdw_bins -1;
            resp2(resp2>bdw_bins -1) = bdw_bins -1;
            jointResp = [resp1;resp2];            
            MI_v(pairi, :) = MI({jointResp,S}, reqOutputs_MI, opts_MI);            
            PID_v(pairi,:) = PID({resp1, resp2,S}, reqOutputs_PID, opts_PID); %order = 'Syn', 'Red', 'Unq1', 'Unq2'
        end
        PID_result.(field_name) = PID_v;
        for bdwIdx = 1:numel(info_bdw_terms)
            bdwLab = info_bdw_terms{bdwIdx};
            for pairi = 1:length(pairlist)
                MI_breakdown.(field_name).(bdwLab)(pairi) = MI_v{pairi,bdwIdx};
            end
        end              
        % -------------------------------------------------------------------------
        % Step 3: Compute Population Information with SVM Decoder
        % % -------------------------------------------------------------------------
        SVM_opts.svm_family = 'linear';
        SVM_out = svm_wrapper({R,S}, reqOutputs_SVM, SVM_opts);
        S_p = SVM_out{1};
        testIdx = SVM_out{2};
        for l = 1:length(testIdx)
            test_index = cell2mat(testIdx(l));
            MI_partitions(l) = cell2mat(MI({S_p(test_index), S(test_index)},reqOutputs_popMI,opts_MI));
        end 
        MI_pop.(field_name).linear = MI_partitions;


        SVM_opts.svm_family = 'RBF';
        SVM_out = svm_wrapper({R,S}, reqOutputs_SVM, SVM_opts);
        S_p = SVM_out{1};
        testIdx = SVM_out{2};
        for l = 1:length(testIdx)
            test_index = cell2mat(testIdx(l));
            MI_partitions(l) = cell2mat(MI({S_p(test_index), S(test_index)},reqOutputs_popMI,opts_MI));
        end 
        MI_pop.(field_name).RBF = MI_partitions;      
        MI_partitions_rbf = [];
        MI_partitions_lin = [];
        % Perform shufflung and compute MI for shuffled data
        for shIdx = 1:nShuff
            inputs_sh = shuffle({R, S}, {'A_B'});
            RSh = inputs_sh{1};
            SVM_opts.svm_family = 'linear';
            SVM_out = svm_wrapper({RSh,S}, reqOutputs_SVM, SVM_opts);
            S_p = SVM_out{1};
            testIdx = SVM_out{2};
            for l = 1:length(testIdx)
                test_index = cell2mat(testIdx(l));
                MI_partitions(l) = cell2mat(MI({S_p(test_index), S(test_index)},reqOutputs_popMI,opts_MI));
            end
             MI_partitions_lin = [MI_partitions_lin, MI_partitions];

            SVM_opts.svm_family = 'RBF';
            SVM_out = svm_wrapper({RSh,S}, reqOutputs_SVM, SVM_opts);
            S_p = SVM_out{1};
            testIdx = SVM_out{2};
            for l = 1:length(testIdx)
                test_index = cell2mat(testIdx(l));
                MI_partitions(l) = cell2mat(MI({S_p(test_index), S(test_index)},reqOutputs_popMI,opts_MI));
            end
            MI_partitions_rbf = [MI_partitions_rbf, MI_partitions];
        end 
        MISh_pop.(field_name).RBF = MI_partitions_lin;
        MISh_pop.(field_name).linear = MI_partitions_rbf;
    end
end
if ~exist('Figure2', 'dir')
    mkdir('Figure2');

end
if ~exist('Figure2/Results', 'dir')
    mkdir('Figure2/Results');
end

filename = sprintf('Figure2/Results/A1_Figure2.mat');
save(filename);
