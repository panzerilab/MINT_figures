% -------------------------------------------------------------------------
% Tutorial: Population Information Analysis
% -------------------------------------------------------------------------
% In this tutorial, we simulate and analyze population responses under
% different correlation scenarios. We perform information breakdown
% analysis of neuron pairs and SVM analysis of the whole population.

% -------------------------------------------------------------------------
% Step 1: Initialize Parameters and Options
% -------------------------------------------------------------------------

clear; clc
rng(0);

nNeurons = 20; 
nPairs = nNeurons*(nNeurons-1)/2;
nTrials_per_stim = 200; 
nShuff = 2; 
nReps = 10; 
tPoints = 500; 

scenario_labels = {'limitNoise','onlyNoise'};
% -------------------------------------------------------------------------
% Define optional Inputs for MI and PID
% ------------------------------------------------------------------------

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

outputs_MI = {'I(A;B)', 'Ilin(A;B)', 'Iss(A)', 'Ici(A;B)', 'Icd(A;B)'};
outputs_PID = {'Joint','PID_atoms'};
outputs_SVM = {'labels', 'testIdx'};
outputs_popMI = {'I(A;B)'};


% Info breakdown opts
info_bdw_terms = {'Joint','ILIN','ISS','ICI','ICD'};
bdw_bins = 5;% number of bins used to discretize spike counts for information breakdown calculations

% -------------------------------------------------------------------------
% Simulation and Analysis
% ------------------------------------------------------------------------
for repIdx = 1:nReps
    disp(['Repetition number ', num2str(repIdx)])
    % -------------------------------------------------------------------------
    % Step 2: Run Simulations and Analye the Results for Each Scenario
    % -------------------------------------------------------------------------

    % Generate Stimulus
    S = [1*ones(1,nTrials_per_stim),2*ones(1,nTrials_per_stim)];
    % Here we perform target-fixed shuffling of the responses
    shIdxs = nan(nShuff,nNeurons,2*nTrials_per_stim);

    for shIdx = 1:nShuff
        for stimVal = 1:2
            for cellIdx = 1:nNeurons
                stim_idxs = find(S == stimVal);
                rand_idxs = randperm(numel(stim_idxs));
                val_idxs = stim_idxs(rand_idxs);
                shIdxs(shIdx,cellIdx,stim_idxs) = val_idxs;
            end
        end
    end
    % Loop through each scenario
    for scIdx = 1:numel(scenario_labels)
        scLab = scenario_labels{scIdx};

        % ---------------------------------------------------------------------
        % Step 2.1: Simulate Spike Trains for each scenario
        % ---------------------------------------------------------------------

        switch scLab
            case 'limitNoise' 
                % Scenario 2: info limiting correlations
                disp(['Simulating ',scenario_labels{1},' scenario'])
                lambda_1 = 0.8;
                lambda_2 = 1.9;
                lambda_noise_1 = 0.2;
                lambda_noise_2 = 0.1;
                noise = true; 
                [R_1, R_2, R] = simulate_spike_trains(nNeurons, tPoints, nTrials_per_stim, lambda_1, lambda_2, lambda_noise_1, lambda_noise_2, scLab);

            case 'onlyNoise'
                disp(['Simulating ',scenario_labels{2},' scenario'])
                lambda_1 = 1;
                lambda_2 = 2;
                lambda_noise_1 = 1;
                lambda_noise_2 = 0;
                noise = true;
                [R_1, R_2, R] = simulate_spike_trains(nNeurons, tPoints, nTrials_per_stim, lambda_1, lambda_2, lambda_noise_1, lambda_noise_2,scLab);
        end

        % -------------------------------------------------------------------------
        % Step 2: Compute PID and Information Breakdown
        % -------------------------------------------------------------------------
        % pairlist = nchoosek(1:nNeurons,2);
        % MI_v = cell(length(pairlist), 5);
        % PID_v = cell(length(pairlist), 5);       
        % parfor pairi = 1:length(pairlist)
        %     cell1 = pairlist(pairi,1);
        %     cell2 = pairlist(pairi,2);
        %     resp1 = R(cell1,:);
        %     resp2 = R(cell2,:);
        %     resp1(resp1>bdw_bins -1) = bdw_bins -1;
        %     resp2(resp2>bdw_bins -1) = bdw_bins -1;
        %     jointResp = [resp1;resp2];            
        %     MI_v(pairi, :) = MI({jointResp,S}, outputs_MI, opts_MI);            
        %     PID_v(pairi,:) = PID({resp1, resp2,S}, outputs_PID, opts_PID); %order = 'Syn', 'Red', 'Unq1', 'Unq2'
        % end
        % PID_result.(scLab)(:,:,repIdx) = PID_v;
        % for bdwIdx = 1:numel(info_bdw_terms)
        %     bdwLab = info_bdw_terms{bdwIdx};
        %     for pairi = 1:length(pairlist)
        %         MI_breakdown.(scLab).(bdwLab)(repIdx,pairi) = MI_v{pairi,bdwIdx};
        %     end
        % end        
        % % -------------------------------------------------------------------------
        % Step 2.3: Compute Population Information with SVM Decoder
        % -------------------------------------------------------------------------
        SVM_opts.svm_family = 'linear';
        SVM_out = SVM({R,S}, outputs_SVM, SVM_opts);
        S_p = SVM_out{1};
        testIdx = SVM_out{2};
        for l = 1:length(testIdx)
            test_index = cell2mat(testIdx(l));
            MI_partitions(l) = cell2mat(MI({S_p(test_index), S(test_index)},outputs_popMI,opts_MI));
        end 
        MI_pop.(scLab)(repIdx).linear = MI_partitions;
        

        SVM_opts.svm_family = 'RBF';
        SVM_out = SVM({R,S}, outputs_SVM, SVM_opts);
        S_p = SVM_out{1};
        testIdx = SVM_out{2};
        for l = 1:length(testIdx)
            test_index = cell2mat(testIdx(l));
            MI_partitions(l) = cell2mat(MI({S_p(test_index), S(test_index)},outputs_popMI,opts_MI));
        end 
        MI_pop.(scLab)(repIdx).RBF = MI_partitions;      

        % Perform shufflung and compute MI for shuffled data
        for shIdx = 1:nShuff
            RSh = cell2mat(shuffle({R, S}), {'A_B'}); 

            SVM_opts.svm_family = 'linear';
            SVM_out = SVM({RSh,S}, outputs_SVM, SVM_opts);
            S_p = SVM_out{1};
            testIdx = SVM_out{2};
            for l = 1:length(testIdx)
                test_index = cell2mat(testIdx(l));
                MI_partitions(l) = cell2mat(MI({S_p(test_index), S(test_index)},outputs_popMI,opts_MI));
            end
            MISh_pop.(scLab)(repIdx).linear = MI_partitions;

            SVM_opts.svm_family = 'RBF';
            SVM_out = SVM({RSh,S}, outputs_SVM, SVM_opts);
            S_p = SVM_out{1};
            testIdx = SVM_out{2};
            for l = 1:length(testIdx)
                test_index = cell2mat(testIdx(l));
                MI_partitions(l) = cell2mat(MI({S_p(test_index), S(test_index)},outputs_popMI,opts_MI));
            end
            MISh_pop.(scLab)(repIdx).RBF = MI_partitions;
        end
    end
end
if ~exist('Figure2/Results', 'dir')
    mkdir('Figure2/Results');
end
filename = sprintf('Figure2/Results/Simulation_Results_%d_%d.mat', nTrials_per_stim, nReps);
save(filename);
% -------------------------------------------------------------------------
% Helper Function: Simulate Spike Trains
% ------------------------------------------------------------------------
function [R_1, R_2, R_count_tot] = simulate_spike_trains(nNeurons, tPoints, nTrials_per_stim, lambda_1, lambda_2, lambda_noise_1,lambda_noise_2, scIdx)
R_1 = zeros(nNeurons,tPoints,nTrials_per_stim);
R_2 = zeros(nNeurons,tPoints,nTrials_per_stim);

if strcmp(scIdx, 'noNoise')
    for cellIdx = 1:nNeurons
        for trialIdx = 1:nTrials_per_stim
            [spikes1] = poisson_spike_gen(1:tPoints, lambda_1/tPoints, 0);
            [spikes2] = poisson_spike_gen(1:tPoints, lambda_2/tPoints, 0);
            R_1(cellIdx,:,trialIdx) = spikes1;
            R_2(cellIdx,:,trialIdx) = spikes2;
        end
    end
else
    for trialIdx = 1:nTrials_per_stim
        [spikes1_noise] = poisson_spike_gen(1:tPoints, lambda_noise_1/tPoints, 0);
        if strcmp(scIdx, 'limitNoise')
            [spikes2_noise] = poisson_spike_gen(1:tPoints, lambda_noise_2/tPoints, 0);
        end
        for cellIdx = 1:nNeurons
            [spikes1] = poisson_spike_gen(1:tPoints, lambda_1/tPoints, 0);
            [spikes2] = poisson_spike_gen(1:tPoints, lambda_2/tPoints, 0);
            R_1(cellIdx,:,trialIdx) = spikes1 + spikes1_noise;
            if strcmp(scIdx, 'limitNoise')
                R_2(cellIdx,:,trialIdx) = spikes2 + spikes2_noise;
            else
                R_2(cellIdx,:,trialIdx) = spikes2;
            end
        end
    end   
end
% Compute single-trial features
R_1_count = squeeze(sum(R_1,2));
R_2_count = squeeze(sum(R_2,2));
R_count_tot = [R_1_count,R_2_count];
end



