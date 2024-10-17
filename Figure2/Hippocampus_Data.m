% -------------------------------------------------------------------------
% Tutorial: Hippocampus_data
% -------------------------------------------------------------------------
clear; clc
rng(0);
load('Figure2/Hippocampus_data.mat');

field_names = fieldnames(Hippocampus_data);

% -------------------------------------------------------------------------
% Define optional Inputs for MI and PID
% ------------------------------------------------------------------------
nShuff = 2;
opts_MI.bias = 'qe';
opts_MI.xtrp = 10;
opts_MI.bin_method = {'none', 'none'};
opts_MI.supressWarnings = true;

opts_PID.bias = 'shuffSub';
opts_PID.shuff = 30;
opts_PID.bin_methodX = 'none';
opts_PID.bin_methodY = 'none';
opts_PID.supressWarnings = true;

SVM_opts.cv = {'KFold', 5};
SVM_opts.libsvm = false;
optimization.optim_cv = {'KFold', 5};
SVM_opts.optim_opts = optimization;

outputs_MI = {'I(A;B)', 'Ilin(A;B)', 'Iss(A)', 'Ici(A;B)', 'Icd(A;B)'};
outputs_PID = {'Joint','PID_atoms'};
outputs_SVM = {'labels'};
outputs_popMI = {'I(A;B)'};

% Info breakdown opts
info_bdw_terms = {'Joint','ILIN','ISS','ICI','ICD'};
bdw_bins = 5;% number of bins used to discretize spike counts for information breakdown calculations
% -------------------------------------------------------------------------
% Step 1: Load Subject data
% -------------------------------------------------------------------------
for i = 1:length(field_names)
    % Extract data
    field_name = field_names{i, 1};
    disp(['Subject ',field_name])
    S = Hippocampus_data.(field_name).Position;
    R = Hippocampus_data.(field_name).R_neuron_binned;
    [nNeurons, nTrials] = size(R);
    nPairs = nNeurons*(nNeurons-1)/2;
    % -------------------------------------------------------------------------
    % Step 2: Compute PID and Information Breakdown
    % -------------------------------------------------------------------------
    pairlist = nchoosek(1:nNeurons,2);
    MI_v = cell(length(pairlist), 5);
    PID_v = cell(length(pairlist), 5);
    parfor pairi = 1:length(pairlist)
        cell1 = pairlist(pairi,1);
        cell2 = pairlist(pairi,2);
        resp1 = R(cell1,:);
        resp2 = R(cell2,:);
        resp1(resp1>bdw_bins -1) = bdw_bins -1;
        resp2(resp2>bdw_bins -1) = bdw_bins -1;
        jointResp = [resp1;resp2];
        MI_v(pairi, :) = MI({jointResp,S}, outputs_MI, opts_MI);
        PID_v(pairi,:) = PID({resp1, resp2,S}, outputs_PID, opts_PID); %order = 'Syn', 'Red', 'Unq1', 'Unq2'
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
        % -------------------------------------------------------------------------
        SVM_opts.svm_family = 'linear';
        SVM_out = SVM({R,S}, outputs_SVM, SVM_opts);
        S_p = SVM_out{1};
        testIdx = SVM_out{2};
        for l = 1:length(testIdx)
            test_index = cell2mat(testIdx(l));
            MI_partitions(l) = cell2mat(MI({S_p(test_index), S(test_index)},outputs_popMI,opts_MI));
        end 
        MI_pop.(field_name).linear = MI_partitions;
        

        SVM_opts.svm_family = 'RBF';
        SVM_out = SVM({R,S}, outputs_SVM, SVM_opts);
        S_p = SVM_out{1};
        testIdx = SVM_out{2};
        for l = 1:length(testIdx)
            test_index = cell2mat(testIdx(l));
            MI_partitions(l) = cell2mat(MI({S_p(test_index), S(test_index)},outputs_popMI,opts_MI));
        end 
        MI_pop.(field_name).RBF = MI_partitions;      

        % Perform shufflung and compute MI for shuffled data
        for shIdx = 1:nShuff
            RSh = cell2mat(shuffle({R}));

            SVM_opts.svm_family = 'linear';
            SVM_out = SVM({RSh,S}, outputs_SVM, SVM_opts);
            S_p = SVM_out{1};
            testIdx = SVM_out{2};
            for l = 1:length(testIdx)
                test_index = cell2mat(testIdx(l));
                MI_partitions(l) = cell2mat(MI({S_p(test_index), S(test_index)},outputs_popMI,opts_MI));
            end
            MISh_pop.(field_name).linear = MI_partitions;

            SVM_opts.svm_family = 'RBF';
            SVM_out = SVM({RSh,S}, outputs_SVM, SVM_opts);
            S_p = SVM_out{1};
            testIdx = SVM_out{2};
            for l = 1:length(testIdx)
                test_index = cell2mat(testIdx(l));
                MI_partitions(l) = cell2mat(MI({S_p(test_index), S(test_index)},outputs_popMI,opts_MI));
            end
            MISh_pop.(field_name).RBF = MI_partitions;
        end 
end
if ~exist('Figure2/Results', 'dir')
    mkdir('Figure2/Results');
end
filename = sprintf('Figure2/Results/Hippocampus_Results.mat');
save(filename);