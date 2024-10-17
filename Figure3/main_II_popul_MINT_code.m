% Load the toolbox and the needed folders
clear all; clc; close all; path(pathdef); 
addpath(genpath('CaImAn-MATLAB-master/')); addpath(genpath('MINT-master/'));

%% Load the raw data
tic; data_for_information_analysis; toc;
tic; generate_spike_train_data; toc;
extract_trials_stim_choice

%% Load the processed data of the spiking activity in significant neurons at their II peaks
load(['Data/spk_trains/DFF_spktrain_alltrials.mat'],'spike_train_alltrials'); 
[spike_train_sig_cells_peak_info_ind] = load_data_II_cells_spk_DFF_alltrials_toneD(spike_train_alltrials);
save([pwd, '/Results/current_analysis/preprocessed_spk_input_data.mat'],'spike_train_sig_cells_peak_info_ind');

%% Decode the stimuli and choice based on the population spiking activity
load(['Data/spk_trains/DFF_spktrain_alltrials.mat'],'spike_train_alltrials'); 
load([pwd, '/Results/current_analysis/preprocessed_spk_input_data.mat'],'spike_train_sig_cells_peak_info_ind');
decode_stim_choice_popul_spk_MINT(spike_train_sig_cells_peak_info_ind, spike_train_alltrials)

%% Plot the results
info_scaling_plot_MINT;