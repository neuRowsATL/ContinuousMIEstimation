%% Load Data objects
% NOTE- INSERT NAME OF DATA OBJECT
load('2020bl21lb21_example_data07292020dataObjs.mat')

save_str = '20200422_bl21lb21_example_analysis';

%Set analysis variables of interest
unit1_name = 'unitD';
unit2_name = 'unitG';

% NOTE- SAVE A COPY OF THE SCRIPT WITH A NEW DATA AND UNIT NAMES AFTER
% REVISING THE VARIABLES ABOVE

verbose_level = 5;


%% Construct CC analysis object and run estimates

% Construct mi_analysis object
a_cc = calc_count_count(d, b, {unit1_name , unit2_name}, 'verbose', verbose_level);

a_cc.buildMIs();

a_cc.calcMIs();

a_cc.getMIs();

save(strcat(save_str, datestr(date, 'mmddyyyy'), unit1_name, unit2_name, '_analysis_c1c2.mat'), 'a_cc')


%% Construct TC analysis object and run estimates

% Construct mi_analysis object
a_tc = calc_timing_count(d, b, {unit1_name , unit2_name}, 'verbose', verbose_level);

a_tc.buildMIs();

a_tc.calcMIs();

a_tc.getMIs();

save(strcat(save_str, datestr(date, 'mmddyyyy'), unit1_name, unit2_name, '_analysis_t1c2.mat')', 'a_tc')

%% Construct TC analysis object and run estimates

% Construct mi_analysis object
a_tc = calc_timing_count(d, b, {unit2_name , unit1_name}, 'verbose', verbose_level);

a_tc.buildMIs();

a_tc.calcMIs();

a_tc.getMIs();

save(strcat(save_str, datestr(date, 'mmddyyyy'), unit1_name, unit2_name, '_analysis_t2c1.mat'), 'a_tc')

%% Construct TT analysis object and run estimates

% Construct mi_analysis object
a_tt = calc_timing_timing(d, b, {unit1_name , unit2_name}, 'verbose', 5);

a_tt.buildMIs();

a_tt.calcMIs();

a_tt.getMIs();

save(strcat(save_str, datestr(date, 'mmddyyyy'), unit1_name, unit2_name, '_analysis_t1t2.mat'), 'a_tt')