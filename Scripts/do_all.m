%% Generate the data for the paper
%{
% generate the transition data
% initialized using the spectral method
% for random initialization, modify corresponding scripts
generate_transition_data_gaussian; % gaussian model
generate_transition_data_unitary; % unitary model
generate_transition_data_cdp; % coded diffraction 1D
generate_transition_data_cdp2d; % coded diffraction 2D

% generate timing data
generate_timing_data_gaussian; % gaussian model
generate_timing_data_cdp; % coded diffraction 1D
generate_timing_data_cdp2d; % coded diffraction 2D

% generate additive noise data
generate_noise_data_gaussian; % gaussian model
generate_noise_data_cdp; % coded diffraction 1D
generate_noise_data_cdp2d; % coded diffraction 2D

% generate molecule and galaxy data
generate_molecule_data_cdp2d; % molecule
generate_galaxy_data_cdp2d; % galaxy

% generate comparison data betweew Kaczmarz and relaxed Kaczmarz
generate_data_kacz_vs_relaxedkacz;
%}

%% Process all the data
% the .txt data files must be on the path set in procee_all_data.m script
%process_all_data

%% Plots all the figures
% the .mat files must be on the path set in the plot scripts
% all the plots will be saved to ../Plots
make_plots_transition_all_algs; % Figure 1
make_plots_transition_diff_init; % Figure 2
make_plots_robust_to_noise; % Figure 3
make_plots_kacz_vs_relaxedkacz; % Figure 4

%% Tex files for the tables
% the .tex files for the tables are located in ../Plots
% the corresponding data files are located in ../Timing and ../Image
