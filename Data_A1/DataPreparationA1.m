clc, clear all;
% -------------------------------------------------------------------------
% This code is preparing the data of A1 for the later analysis in Figure 2
% and 3 of the MINT toolbox paper the data can be downloaded here:
% https://drum.lib.umd.edu/items/30d43732-7149-4726-a860-0ae3d210b2ae
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Step 1: Load the data
% -------------------------------------------------------------------------
data = dir(['Data/*.mat']);
for ind  =  1:length(data)
    data_mouse = load(['Data/',data(ind).name]); %%% Load all the data
    Behavior = data_mouse.Behavior;              %%% Load only trials during the task
    % Un-filtered DF
    DFF_corrected_Beh = data_mouse.Behavior.Fluorescence.DeltaF;
    %%%%%%%%%%%%%%%%%%
    % Butterworth filter DFF
    % Define Butterworth filter parameters
    order = 2;  cutoff_frequency = 2;
    nyquist_frequency = data_mouse.Behavior.Fluorescence.fps / 2;
    normalized_cutoff_frequency = cutoff_frequency / nyquist_frequency;
    [b, a] = butter(order, normalized_cutoff_frequency,'low');
    DFF_smoothed_Beh = nan(size(DFF_corrected_Beh));
    for i = 1:size(DFF_corrected_Beh,2)
        for j = 1:size(DFF_corrected_Beh,3)
            if sum(isnan(data_mouse.Behavior.Fluorescence.DeltaF(:,i,j)))==0
                DFF_smoothed_Beh(:,i,j) = filtfilt(b, a, data_mouse.Behavior.Fluorescence.DeltaF(:,i,j));
            else
                DFF_smoothed_Beh(:,i,j) = nan;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%
    brightCells             = Behavior.Fluorescence.brightCells'==1;            %%% Extract brightcells
    n_responsive_neurons    = size(DFF_corrected_Beh(:,:,brightCells),3);       %%% Responsive neurons
    timepoints              = size(DFF_corrected_Beh(:,:,brightCells),1);       %%% Size of timepoints
    name_exp_mice           = Behavior.Animal;                                  %%% Name of the individual
    cell_location_mice      = Behavior.Fluorescence.cellLocation;               %%% XY coordinates of all neurons per experiment
    stimonset               = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'StimOnset'))));          %%% A vector of time frames when a stimulus occured in different trials (at tf=31)
    StimFrequency           = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'StimFrequency'))));      %%% A vector of frequencies in different trials (at tf=31)
    early_trials            = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'Early'))));              %%% A vector of early trials (at tf=31, 1==early trials, 0==otherwise)
    early_FA_trials         = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'EarlyFalseAlarm')))); %%% A vector of early false alarms (at tf=31, 1==early false alarms, 0==otherwise)
    early_hit_trials        = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'EarlyHit')))); %%% A vector of early hits (at tf=31, 1==early hits, 0==otherwise)
    hit_trials              = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'Hit')))); %%% A vector of hit trials (at tf=31, 1==hit trials, 0==otherwise)
    miss_trials             = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'Miss')))); %%% A vector of miss trials (at tf=31, 1==miss trials, 0==otherwise)
    FA_trials               = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'FalseAlarm')))); %%% A vector of false alarms (at tf=31, 1==false alarms, 0==otherwise)
    CR_trials               = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'CorrectReject')))); %%% A vector of correct rejection trials (at tf=31, 1==correct rejection trials, 0==otherwise)
    FirstResponse           = squeeze(data_mouse.Behavior.Tags(:,:,find(strcmpi(data_mouse.Behavior.TagNames(:,2),'FirstResponse')))');        %find where the first lick occured relative to the image frame timing
    save(['Data/Data_processed/data_processed_',char(data(ind).name)],...
        'DFF_corrected_Beh','DFF_smoothed_Beh',...
        'stimonset','StimFrequency','early_trials','early_FA_trials', 'early_hit_trials',...
        'hit_trials','miss_trials','FA_trials','CR_trials','brightCells', 'name_exp_mice',...
        'cell_location_mice','n_responsive_neurons','timepoints','FirstResponse');    %%% Save (what I need) to Data/data_for_information_analysis folder
end
clear
%% ------------------------------------------------------------------------
% Step 2:  DFF deconvolution (only brightcells; remove NaN trials) 
% -------------------------------------------------------------------------
% You need CAImAn toolbox for this step
processed_original_data = dir(['Data/Data_processed/*.mat']);
alldata = cell(1,length(processed_original_data));
for ind = 1:length(processed_original_data)
    alldata{ind} = load(['Data/Data_processed/', char(processed_original_data(ind).name)]);
end
type_data_name = {'corrected','smoothed'};
DFF_alltrials = struct;
spike_train_alltrials = struct;
deconv_type_name = {'foopsi_ar1'};
for deconv_type = 1
    deconv_type
    for ind = 1:length(processed_original_data)
        ind
        for type_data = 1:length(type_data_name)
            DFF_data_Beh = {alldata{ind}.DFF_corrected_Beh, alldata{ind}.DFF_smoothed_Beh};
            nan_trials = sum(sum(isnan(DFF_data_Beh{type_data}(:,:,alldata{ind}.brightCells == 1))),3);
            DFF_alltrials(ind,type_data).DFF_alltrials = DFF_data_Beh{type_data}(:, :, alldata{ind}.brightCells == 1);
            T = size(DFF_alltrials(ind,type_data).DFF_alltrials,1);                     %%% Timepoints
            N = size(DFF_alltrials(ind,type_data).DFF_alltrials,3);                     %%% Neurons
            DFF_alltrials(ind,type_data).name_exp_mice = alldata{ind}.name_exp_mice;
            DFF_alltrials(ind,type_data).N = N;
            DFF_alltrials(ind,type_data).type_data_name = type_data_name(type_data);
            spike_train_alltrials(ind,type_data,deconv_type).spike_train_alltrials = deconv_DFF(DFF_alltrials(ind,type_data).DFF_alltrials, deconv_type_name{deconv_type});
            spike_train_alltrials(ind,type_data,deconv_type).name_exp_mice = alldata{ind}.name_exp_mice;
            spike_train_alltrials(ind,type_data,deconv_type).N = N;
            spike_train_alltrials(ind,type_data,deconv_type).deconv_type = deconv_type_name{deconv_type};
        end
    end
end
save('DFF_spktrain_alltrials.mat','DFF_alltrials','spike_train_alltrials','-v7.3');

%% ------------------------------------------------------------------------
% Step 3: extract Stimulus and Choice
% -------------------------------------------------------------------------
processed_original_data = dir(['Data/Data_processed/*.mat']);  %%% Consider all experiments and mice
type_data_name = {'corrected','smoothed'};
stim_choice = struct;
for ind = 1:length(processed_original_data)  %%% Generalize with a parfor here (Load issue not save)
    ind
    for type_data = 1:length(type_data_name)
        load(['Data/Data_processed/', char(processed_original_data(ind).name)]);                        %%% Load file with Fluorescence responses, stimuli and choices
        stim_choice(ind,type_data).CellLoc = cell_location_mice(brightCells == 1,:);
        stim_choice(ind,type_data).name_exp_mice = name_exp_mice;
        DFF_data_Beh = {DFF_corrected_Beh, DFF_smoothed_Beh};
        DeltaF_trials_stim  = DFF_data_Beh{type_data}(:, :, brightCells == 1);                          %%% Consider only responsive cells at a given trial
        stim_choice(ind,type_data).nan_trials  = sum(sum(isnan(DeltaF_trials_stim)),3);
        stim_choice(ind,type_data).hits_vec  = hit_trials(31,:);                                        %%% Hits = c1
        stim_choice(ind,type_data).FA_vec = FA_trials(31,:);                                            %%% FA = c1
        stim_choice(ind,type_data).early_vec = early_trials(31,:);                                      %%% FA = c1
        stim_choice(ind,type_data).early_hit_vec = early_hit_trials(31,:);                              %%% FA = c1
        stim_choice(ind,type_data).early_FA_vec = early_FA_trials(31,:);                                %%% FA = c1
        miss_vec = miss_trials(31,:);                                                                   %%% Miss = c2
        miss_vec(miss_vec == 1) = 2;                                                                    %%% Re-label miss = c2
        stim_choice(ind,type_data).miss_vec  = miss_vec;                                   
        CR_vec = CR_trials(31,:);                                       
        CR_vec(CR_vec == 1) = 2;                                                                        %%% Re-label CR = c2
        stim_choice(ind,type_data).CR_vec = CR_vec;
        c_vec = stim_choice(ind,type_data).hits_vec + stim_choice(ind,type_data).FA_vec + stim_choice(ind,type_data).miss_vec + stim_choice(ind,type_data).CR_vec; %%% choice vector (without the early trials) for the information toolbox
        stim_choice(ind,type_data).c_vec = c_vec;
        stim_choice(ind,type_data).choice_vector = c_vec(c_vec'~=0);                                    %%% Select those neurons that have made a choice (H,M,FA,CR) (no early trials)
        s_vec = StimFrequency(31,:);                        
        s_vec(s_vec<11000) = 1;                                                                         %%% Re-label low freqs  = s1
        s_vec(s_vec>11000) = 2;                                                                         %%% Re-label high freqs = s2
        stim_choice(ind,type_data).stimulus_vector = s_vec(c_vec'~=0);                                  %%% Select those stimuli that correspond to the choice (H,M,FA,CR) (no early trials)
        stim_choice(ind,type_data).considered_trials = sum(c_vec'~=0);                                  %%% Count number of trials with a choice (H,M,FA,CR)
        stim_choice(ind, type_data).FirstResponse = FirstResponse;
    end
end
save('stim_choice.mat','stim_choice');

%% ------------------------------------------------------------------------
% Step 3: Load Intersection Information peak
% -------------------------------------------------------------------------
load('DFF_spktrain_alltrials.mat','spike_train_alltrials');

GC_type  = {'H','M','C','F'}; time_scale_GC = 1; %%%time_scale_GC = 1: short time-scales;
load([pwd,'/GCNetOutputs_10tp_earlyIIpeak']);
deconv_type = 1; type_data = 2; bin_size = 10; type_sig_cells = 'Sig_cells';
load([pwd,'/stim_choice.mat']);
load([pwd,'/Neurons_class_GC_sliding_window_',char(type_sig_cells),'_',num2str(bin_size),'tp.mat'],'II_CLASS_PEAK_value_mice');
load([pwd,'/significant_neurons_sliding_window_',num2str(bin_size),'tp.mat'],'significant_neurons');
spike_train_sig_cells_peak_info_ind = cell(1,length(II_CLASS_PEAK_value_mice));
for ind   = 1:length(II_CLASS_PEAK_value_mice)
    if  length(significant_neurons{ind}) > 20
        GC_neurons  = GCNets.(GC_type{1}){time_scale_GC,ind}.cellids;
        c_vec = stim_choice(ind,type_data).c_vec;
        n_trials  = length(find(c_vec'~=0));
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