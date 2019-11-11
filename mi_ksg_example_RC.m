%% THINK ABOUT THIS
% 1. Upload data and instantiate data object
% 2. Format data (i.e., phase, timing, etc.)
% 3. Run analysis class
% 4. Make plots
% 5. Repeat 2 - 4 for different variable/formatting combinations

% For auditing purposes, we probably need to include in the analysis
% object, a single struct that copies the parameters of data parameters adn
% stores it for good will

% It might actually be most advantageous to create another class.... for
% storing results objects --> this would include, sims objects and
% parameters.

%% Upload Data and Instantiate Data Object

% Clear variables from workspace
clear

% Load spiking data
load('bl21lb21_171218_dtvw-0002.mat')

% Extract spike train of interest and convert times to ms
ts = spikedata.ts(:,1)*1000;

% Clear irrelevant data from workspace
clearvars -except ts

% Note- we will need to add a line of code to upload corresponding pressure
% data. 

% Load pressure data (THIS DATA DOES NOT CORRESPOND TO THE SPIKES)
%load('bl21lb21_171218_125431_300s_CH1-2-3-4-5-6-7-8_Fs30000_Filt300-7500-analog.mat')

% Extract pressure wave
%pressure = wav_tmp(9,1:end);

% Instantiate data object
neural_data = mi_data(30000,30000);

% Add spiking data to data object
neural_data.add_spikes(ts, 'n_timebase', 'phase');

%% Add pressure wave to object
neural_data.add_behavior(pressure);


%% Instantiate Analysis Class

% This is the level where you would switch which type of analysis you want
% to run or 
% Instantiate analysis object(s)
% The subclass below calculates consecutive ISI MI.
MI = calc_isi_isi(neural_data, [1], 1, 200);

%% DATA ANALYSIS SAMPLE SCRIPT
close('all');

% Load spike times
% BC: 
% load('bl21lb21_171218_spikedata-0002-CH5.mat');

% RC:
clear
load('bl21lb21_171218_dtvw-0002.mat')
ts = spikedata.ts(:,1)*1000;
clearvars -except ts


% Load pressure data
load('bl21lb21_171218_125431_300s_CH1-2-3-4-5-6-7-8_Fs30000_Filt300-7500-analog.mat')
pressure = wav_tmp(9,1:end);

% Instantiate data object
neural_data = mi_data(30000,30000);

% Add neural data
add_spikes(neural_data, ts);

% Add cyclic pressure data
testObj.add_behavior(pressure);




%% Instantiate analysis object(s)
MI_isi = calc_isi_isi(neural_data, [1], 1, 200);

buildMIs(MI_isi);

% --> update k-values, etc.

MI_isi.arrMIcore{1}.opt_k = 0;
MI_isi.arrMIcore{1}.k_values = 3:8;

% Run mutual information calculations
calcMIs(MI_isi);


%% MAKE DATA PLOTS
core1 = MI_isi.arrMIcore{1};
viz = mi_ksg_viz;

make_ksg_graph(viz, core1);
plot_k_dependence(viz, core1);

for i=1:length(core1.k_values)
    plot_data_fraction(viz, core1, core1.k_values(i));
end


a = core1.k_values;
b = 1:10;

m = cell2mat(core1.mi_data(:,1));
c = reshape(m, length(b), length(a));

figure();
surf(a, b, c, 'FaceAlpha', 0.5);
hold on;
plot3(core1.k_values, ones(1,length(core1.k_values)), c(1,:), 'k-', 'LineWidth', 5);

xlabel('k-value');
ylabel('Data Fraction (1/N)');
zlabel('MI (bits)');

%% COUNT-COUNT EXAMPLE - working as of 20190129
% RC:
clear
load('bl21lb21_171218_dtvw-0002.mat')
ts1 = spikedata.ts(:,1)*1000;

load('bl21lb21_171218_dtvw-0003.mat')
ts2 = spikedata.ts(:,1)*1000;

load('bl21lb21_171218_dtvw-pressure.mat')
ptimes = spikedata.ts(:,1)*1000;

clearvars -except ts1 ts2 ptimes

% Instantiate data object
neural_data = mi_data(30000,30000);
add_spikes(neural_data, ts1);
add_spikes(neural_data, ts2);

% Add pressure times
neural_data.cycleTimes = ptimes;

clearvars -except neural_data


% Instantiate analysis object(s)
MI_cc = calc_count_count(neural_data, [1 2]);

buildMIs(MI_cc, MI_cc.verbose);

calcMIs(MI_cc);



%%

if size(core1.mi_data,1) > length(core1.k_values)
    a = core1.k_values;
    b = 1:10;

    m = cell2mat(core1.mi_data(:,1));
    c = reshape(m, length(b), length(a));

    figure();
    % plot mi surface
    surf(a, b, c, 'FaceAlpha', 0.5);
    hold on;
    
    % plot mi along data fracs = 1
    plot3(core1.k_values, ones(1,length(core1.k_values)), c(1,:), 'k-', 'LineWidth', 5);

    % plot optimal k-value
    k_ix = core1.k_values == core1.opt_k;
    k_ixs = [core1.mi_data{:,4}] == core1.opt_k;
    plot3(ones(1,length(b))*core1.opt_k, 1:10, c(:,k_ix), 'r-', 'LIneWidth', 5);
    
    xlabel('k-value');
    ylabel('Data Fraction (1/N)');
    zlabel('MI (bits)');
end
