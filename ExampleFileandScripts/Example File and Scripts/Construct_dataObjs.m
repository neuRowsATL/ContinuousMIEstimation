%% LOAD DATA
load('20200422_bl21lb21_example_data.mat')

%% SET RELEVANT VARIABLES
% Note- for your data, you make have different info for neurons and
% pressure data
standard_fileInfo = '20200422_bl21lb21_example_data';

% Set behav sample frequency
bFs = 30000;

% Set neural sample frequency
nFs = 30000;

behav_fileInfo = strcat('example script/', datestr(date, 'mmddyyyy'));

neural_fileInfo = strcat('example script/', datestr(date, 'mmddyyyy'));

verbose_level = 5;

%% MAKE PRESSURE DATA OBJECT

% Construct pressure behavioral object
b = mi_data_pressure(behav_fileInfo, 'verbose', verbose_level);

% Add cycle times to pressure object
b.add_cycleTimes(final_cycle_times', strcat(standard_fileInfo, datestr(date, 'mmddyyyy')), bFs);

%% MAKE NEURAL DATA OBJECT

% Set up info strings for each unit
str_UnitG = strcat(standard_fileInfo, '/unitG');
str_UnitD = strcat(standard_fileInfo, '/unitD');

% Construct neural data object
d = mi_data_neural(neural_fileInfo, 'verbose', verbose_level);

% Add spikes to neural data object
d.add_spikes(unitG, str_UnitG, nFs, 'unitG');
d.add_spikes(unitD, str_UnitD, nFs, 'unitD');

%% Save Data Objects
save(strcat(standard_fileInfo, datestr(date, 'mmddyyyy'), 'dataObjs.mat'), 'd', 'b')

