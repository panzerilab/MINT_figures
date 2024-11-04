% -------------------------------------------------------------------------
% Tutorial: Combined Usage of Information Theoretic tools
% -------------------------------------------------------------------------
% DESCRIPTION 

% -------------------------------------------------------------------------
% Step 1: Initialize Parameters and Options
% -------------------------------------------------------------------------

clc, clear all;
rng(0);

nShuff = 100;
nReps = 10; 
nPopulations = 4;
nSubpopulations = 2;
nTrials = 200;
nTimepoints = 30;
delay = 5;
stim_window = [3 12];
stim_window_2 = [stim_window(1)+12 stim_window(2)+12]; 
null_samples = 300;
delay_shifts = 1:19;

% Initialize opts
MI_opts.bias = 'shuffSub';
MI_opts.shuff = 30;
MI_opts.parallel = 0;
MI_opts.n_bins = {3};
MI_opts.supressWarnings = true;
MI_opts.bin_method = {'eqpop','none'};

TE_opts.bias = 'shuffSub';
TE_opts.shuff = 30;
TE_opts.n_bins = {3}; 
TE_opts.bin_method = {'eqpop','eqpop'};
TE_opts.tau = {-delay};
TE_opts.supressWarnings = true;
TE_opts.computeNulldist = true;
TE_opts.singleTimepoint = true;
TE_opts.n_samples = nShuff;

FIT_opts.bias = 'shuffSub';
FIT_opts.shuff = 30;
FIT_opts.n_bins = {3}; 
FIT_opts.computeNulldist = true;
FIT_opts.n_samples = nShuff;
FIT_opts.shuffling = {'A','A_C'};
FIT_opts.dim_shuffle = {'Trials'};
FIT_opts.supressWarnings = true;
FIT_opts.bin_method = {'eqpop','eqpop', 'none'};
FIT_opts.singleTimepoint = true;

FIT12_1_shift = zeros(length(delay_shifts),(nTimepoints-10),nReps);
FIT12_2_shift = zeros(length(delay_shifts),(nTimepoints-10),nReps);
FIT21_2_shift = zeros(length(delay_shifts),(nTimepoints-10),nReps);
FIT1 = zeros(nPopulations, nPopulations, nReps);
FIT2 = zeros(nPopulations, nPopulations, nReps);
TE_v = zeros(nPopulations, nPopulations, nReps);

for repIdx = 1:nReps
    disp(['Repetition number ', num2str(repIdx)]);
    X = zeros(nPopulations, nSubpopulations, nTimepoints, nTrials);
    S1 = randi(2,1,nTrials);
    S2 = randi(2,1,nTrials);
    epsX = 0.5;
    for t=1:nTimepoints
        X(:,:,t,:) = epsX*randn(nPopulations, nSubpopulations,1, nTrials);
        if (t >stim_window(1))  && (t <=  stim_window(2))
            X(1,1,t,:) = squeeze(X(1,1,t,:)) +epsX*randn+2*(S1-1.5)';
            X(2,2,t,:) = squeeze(X(2,2,t,:)) +epsX*randn+2*(S2-1.5)';
        end
        if (t > stim_window_2(1))  && (t <=  stim_window_2(2))
            X(4,1,t,:) = squeeze(X(4,1,t,:)) +epsX*randn+2*(S1-1.5)';
        end
        if t > delay
            X(1,2,t,:) = X(1,2,t,:) + X(2,2,t-delay,:);
            X(2,1,t,:) = X(2,1,t,:) + X(1,1,t-delay,:);
            X(3,1,t,:) = X(3,1,t,:) + X(1,1,t-delay,:);
            X(4,2,t,:) = X(4,2,t,:) + 3*X(3,2,t-delay,:); % noise transfer
        end
    end
   
    X = squeeze(sum(X,2));

    for neuron = 1:nPopulations
        infoS1(neuron, :, repIdx) = cell2mat(MI({X(neuron,:,:), S1}, {'I(A;B)'},MI_opts));
        infoS2(neuron, :, repIdx) = cell2mat(MI({X(neuron,:,:), S2}, {'I(A;B)'},MI_opts));
    end

    for d = 1:length(delay_shifts)
        delay_shift = delay_shifts(d);
        FIT_opts.tau = {delay_shift};
        FITshuff_12_1_simple_slice = zeros(nTimepoints-10, FIT_opts.n_samples);
        FITshuff_12_1_cond_slice = zeros(nTimepoints-10, FIT_opts.n_samples);
        FITshuff_12_2_simple_slice = zeros(nTimepoints-10, FIT_opts.n_samples);
        FITshuff_12_2_cond_slice = zeros(nTimepoints-10, FIT_opts.n_samples);
        FITshuff_21_2_simple_slice = zeros(nTimepoints-10, FIT_opts.n_samples);
        FITshuff_21_2_cond_slice = zeros(nTimepoints-10, FIT_opts.n_samples);
        parfor t=1:(nTimepoints-10)
            FIT_localopts = FIT_opts;
            FIT_localopts.tpres = {t};
            if t-delay_shift > 0
                [FIT_1_tmp, ~, FIT_1_nullDist] = FIT({X(1,:,:),X(2,:,:),S1}, {'FIT(A->B;S)'},FIT_localopts);               
                FIT12_1_shift(d, t, repIdx) = FIT_1_tmp{1};                              
                FITshuff_12_1_simple = cell2mat(FIT_1_nullDist{1});
                FITshuff_12_1_cond = cell2mat(FIT_1_nullDist{2});                        
                FITshuff_12_1_simple_slice(t,:) =  FITshuff_12_1_simple(:); 
                FITshuff_12_1_cond_slice(t,:) =  FITshuff_12_1_cond(:);              

                [FIT_2_tmp, ~, FIT_2_nullDist] = FIT({X(1,:,:),X(2,:,:),S2}, {'FIT(A->B;S)', 'FIT(B->A;S)'},FIT_localopts);
                FIT12_2_shift(d, t, repIdx) = FIT_2_tmp{1};
                FIT21_2_shift(d, t, repIdx) = FIT_2_tmp{2};                          
                FITshuff_2_simple = cell2mat(FIT_2_nullDist{1});
                FITshuff_2_cond = cell2mat(FIT_2_nullDist{2});      
                FITshuff_12_2_simple_slice(t,:) = FITshuff_2_simple(:,1); 
                FITshuff_12_2_cond_slice(t,:) = FITshuff_2_cond(:,1);
                FITshuff_21_2_simple_slice(t,:) = FITshuff_2_simple(:,2); 
                FITshuff_21_2_cond_slice(t,:) = FITshuff_2_cond(:,2);               
            end
        end
        FIT12_1_sh.simple(d,:,repIdx,:) = FITshuff_12_1_simple_slice;
        FIT12_1_sh.cond(d,:,repIdx,:) = FITshuff_12_1_cond_slice;        
        FIT12_2_sh.simple(d,:,repIdx,:) = FITshuff_12_2_simple_slice;
        FIT12_2_sh.cond(d,:,repIdx,:) = FITshuff_12_2_cond_slice;        
        FIT21_2_sh.simple(d,:,repIdx,:) = FITshuff_21_2_simple_slice;
        FIT21_ 2_sh.cond(d,:,repIdx,:) = FITshuff_21_2_cond_slice;
    end
    
    % compute FIT and TE for all pairs of populations
    FIT_opts.tau = {delay};
    FIT_opts.tpres = {7 + delay};
    TE_opts.tau = {delay};
    TE_opts.tpres = {7 + delay};
    for pop1 = 1:nPopulations
       for pop2 = pop1:nPopulations
            if pop1 ~= pop2
                 [FIT_values, ~, FIT_nullDist] = FIT({X(pop1,:,:),X(pop2,:,:),S1},{'FIT(A->B;S)', 'FIT(B->A;S)'},FIT_opts);               
                  FIT1(pop1, pop2, repIdx) = FIT_values{1};
                  FIT1(pop2, pop1, repIdx) = FIT_values{2};
                  FITshuff_simple = FIT_nullDist.A;
                  FITshuff_cond = FIT_nullDist.A_C;
                  FIT1_sh.simple(pop1, pop2, repIdx, :) =  FITshuff_simple(:,1);
                  FIT1_sh.simple(pop2, pop1, repIdx, :) =  FITshuff_simple(:,2);
                  FIT1_sh.conditioned(pop1, pop2, repIdx, :) =  FITshuff_cond(:,1);
                  FIT1_sh.conditioned(pop2, pop1, repIdx, :) =  FITshuff_cond(:,2);

                  [FIT_values, ~, FIT_nullDist] = FIT({X(pop1,:,:),X(pop2,:,:),S2},{'FIT(A->B;S)', 'FIT(B->A;S)'},FIT_opts);
                  FIT2(pop1, pop2, repIdx) = FIT_values{1};
                  FIT2(pop2, pop1, repIdx) = FIT_values{2};
                  FITshuff_simple = FIT_nullDist.A;
                  FITshuff_cond = FIT_nullDist.A_C;
                  FIT2_sh.simple(pop1, pop2, repIdx, :) =  FITshuff_simple(:,1);
                  FIT2_sh.simple(pop2, pop1, repIdx, :) =  FITshuff_simple(:,2); 
                  FIT2_sh.conditioned(pop1, pop2, repIdx, :) =  FITshuff_cond(:,1);
                  FIT2_sh.conditioned(pop2, pop1, repIdx, :) =  FITshuff_cond(:,2);

                  [TE_values, ~, TE_nullDist] =TE({X(pop1,:,:),X(pop2,:,:)},{'TE(A->B)','TE(B->A)'}, TE_opts);
                  TE_v(pop1, pop2, repIdx) = TE_values{1};
                  TE_v(pop2, pop1, repIdx) = TE_values{2};
                  TEshuff_simple = TE_nullDist;
                  TE_sh.simple(pop1, pop2, repIdx, :) =  TEshuff_simple(:,1);
                  TE_sh.simple(pop2, pop1, repIdx, :) =  TEshuff_simple(:,2);                
            end
        end
    end     
end 
if ~exist('Results', 'dir')
    mkdir('Results');
end
currentDateTime = datetime('now', 'Format', 'dd_MM_yyyy');
filename = sprintf('Results/Simulation_Figure4_%s.mat', currentDateTime);
save(filename);


