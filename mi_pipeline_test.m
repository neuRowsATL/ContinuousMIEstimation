%%
% This script will systematically test the mutual information analysis
% pipeline in order to identify errors and bugs.
%
% This script requires data files in the TestData folder
%%
% Instantiate data objects

clear all
close('all')


global_errs = {};

diary mi_pipeline_test_diary.txt
diary on

with_plots = true;

%% RUN MI_DATA
load('TestData/20191018_bl21lb21_171218_spikes.mat');
unit1 = unit1*1000.;
unit2 = unit2*1000.;
unit3 = unit3*1000.;

str_unit1 = 'TestData/20191018_bl21lb21_171218_spikes.mat/unit1';
str_unit2 = 'TestData/20191018_bl21lb21_171218_spikes.mat/unit2';
str_unit3 = 'TestData/20191018_bl21lb21_171218_spikes.mat/unit3';

load('TestData/20191018_bl21lb21_171218_cycles.mat');
cycle_times = [cycle_times(1:end-1)' cycle_times(2:end)']; % Needs to be N x 2 matrix of [on off] x N

str_cycles = 'TestData/20191018_bl21lb21_171218_cycles.mat';

%%
try
    disp([newline newline]);
    clear d;
    disp([newline '===== ===== ===== ===== =====']);
    disp(['RUNNING: mi_data() with ID only' newline newline]);
    d = mi_data('test');
    
    
    % CHECK OBJECT FOR INSTANTIATION CONSISTENCY
    success = (['----- ----- ----- ----- -----' newline 'DEBUG REPORT:' newline]);
    
    add_data(d, unit1, str_unit1, 30000);
    
    success = [success newline 'Imported: data'];
    if ~all(size(d.data.noname.data) == size(unit1)); success = [success '  >> FAILED']; end
    
    
    success = [success newline 'Pulled: data'];
    if ~all(size(get_data(d)) == size(unit1)); success = [success ' >> FAILED']; end
    
    disp(success);
catch e
    e
    global_errs{end+1} = {'Intantiating mi_data with ID only'};
    disp([newline 'ERROR: Unable to instantiate mi_data with ID only']);
end
%% - FOR RUNNING IN EMACS
try
    disp([newline newline]);
    clear d;
    disp([newline '===== ===== ===== ===== =====']);
    disp(['RUNNING: mi_data()' newline newline]);
    d = mi_data('test', 'verbose', 5);

    add_data(d, unit1, str_unit1, 30000, 'unit1');
    add_data(d, unit2, str_unit2, 30000, 'unit2');
    add_data(d, unit3, str_unit3, 30000, 'unit3');


    % CHECK OBJECT FOR INSTANTIATION CONSISTENCY
    success = (['----- ----- ----- ----- -----' newline 'SUCCESSFUL:' newline]);

    % Check that ID matches what was set
    success = [success newline 'Assigned: ID'];
    if ~strcmp(d.ID, 'test'); success = [success '>> FAILED']; end
    
    % Check for correct Fs 
    success = [success newline 'Assigned: Fs'];
    if ~d.Fs == 30000;  success = [success '>> FAILED']; end
    
    % Check for correct verbose
    success = [success newline 'Assigned: verbose'];
    if ~d.verbose == 5;  success = [success '>> FAILED']; end
    
    % Check for correct data names
    success = [success newline 'Assigned: data'];
    if ~isfield(d.data, 'unit1') || ~isfield(d.data, 'unit2') || ~isfield(d.data, 'unit3');  success = [success '>> FAILED']; end

    % Check for correct data sizes
    success = [success newline 'Imported: data'];
    if any(size(d.data.unit1.data) ~= size(unit1)) || ...
            any(size(d.data.unit2.data) ~= size(unit2)) || ...
            any(size(d.data.unit3.data) ~= size(unit3))
         success = [success '>> FAILED'];
    end
    
    % Check for correct data info
    success = [success newline 'Imported: info'];
    if ~strcmp(d.data.unit1.info,str_unit1) || ...
            ~strcmp(d.data.unit2.info,str_unit2) || ...
            ~strcmp(d.data.unit3.info,str_unit3)
        success = [success '>> FAILED'];
    end
    
    % Check for correct data from get_data
    success = [success newline 'Pulled: data'];    
    if any(size(get_data(d,'unit1')) ~= size(unit1)) || ...
            any(size(get_data(d,'unit2')) ~= size(unit2)) || ...
            any(size(get_data(d,'unit3')) ~= size(unit3))
        success = [success '>> FAILED'];
    end
    
    disp(success);
catch e
    e
    global_errs{end+1} = {'Instantiating mi_data with ID and verbose'};
    % Not possible to proceed without mi_data class
    
    disp([newline newline '===== ===== ===== ===== =====' newline 'GLOBAL ERRORS' newline]);
    for i=1:length(global_errs)
        disp(global_errs{i});
    end
    disp(['----- ----- ----- ----- -----' newline]);
    
    error('FATAL ERROR: Unable to construct mi_data objects');
end


%%
try
    clear d
    disp([newline '===== ===== ===== ===== =====']);
    disp(['RUNNING: mi_data_neural()' newline newline]);

    d = mi_data_neural('test', 'verbose', 5);

    add_spikes(d, unit1, str_unit1, 30000, 'unit1');


    
   % CHECK OBJECT FOR INSTANTIATION CONSISTENCY
success = (['----- ----- ----- ----- -----' newline 'SUCCESSFUL:' newline]);

    % Check for correct ID
    success = [success newline 'Assigned: ID'];
    if ~strcmp(d.ID, 'test'); success = [success '>> FAILED']; end
    
    % Check for correct Fs
    success = [success newline 'Assigned: Fs'];
    if ~d.Fs == 30000;  end
    
    % Check for correct verbose
    success = [success newline 'Assigned: verbose'];
    if ~d.verbose == 5; success = [success '>> FAILED']; end
    
    % Check for correct data name
    success = [success newline 'Assigned: data'];
    if ~isfield(d.data, 'unit1'); success = [success '>> FAILED']; end

    % Check for correct data size
    success = [success newline 'Imported: data'];
    if any(size(d.data.unit1.data) ~= size(unit1)); success = [success '>> FAILED']; end

    % Check for correct data info
    success = [success newline 'Imported: info'];
    if ~strcmp(d.data.unit1.info,str_unit1); success = [success '>> FAILED']; end

    % Check for correct data from get_data
    success = [success newline 'Pulled: data'];
    if any(size(get_data(d,'unit1')) ~= size(unit1)); success = [success '>> FAILED']; end
    
    % Check for correct data from get_spikes, 'raw'
    success = [success newline 'Pulled: raw data'];
    if any(size(get_spikes(d, 'format', 'raw', 'name', 'unit1')) ~= size(get_data(d,'unit1'))); success = [success '>> FAILED']; end

       % VALIDATION CHECK - Check by plotting and functions
    % 1. Number of non-nan values = # spikes - # spikes outside of cycles
    % 2. Save validated matrix to check against
    
    c1 = get_count(d, cycle_times, 'unit1');
    
    % Check for correct data from get_count
    success = [success newline 'Pulled: count data'];
    if sum(c1) ~= (sum(~isnan(unit1)) - sum(unit1 < cycle_times(1,1) | unit1 > cycle_times(end,2))); success = [success '>> FAILED']; end
    
    c2 = get_spikes(d, 'format', 'count', 'cycleTimes', cycle_times, 'name', 'unit1');

    % Check for correct data from get_spikes, count
    success = [success newline 'Pulled: count data'];
    if sum(c2) ~= (sum(~isnan(unit1)) - (sum(unit1 < cycle_times(1,1) | unit1 > cycle_times(end,2)))); success = [success '>> FAILED']; end

    % Check for matching data from get_count and from get_spikes, 'count'
    success = [success newline 'Matched: spike counts'];
    if c1 ~= c2; success = [success '>> FAILED']; end
    
    if with_plots
        figure();
        plot(c1);
        hold on;
        plot(c2);
        xlabel('Cycle Index');
        ylabel('Spike Count');
        title('Comparison of Spike Counts');
    end
    
    % VALIDATION CHECKS
    % 1. Number of non-nan values = # spikes - # spikes outside of cycles
    % 2. Save validated matrix to check against
    
    t1 = get_timing(d, cycle_times, 'timeBase', 'time', 'name', 'unit1');
    
    % Check for correct data from get_timing
    success = [success newline 'Pulled: count data'];
    if sum(sum(~isnan(t1))) ~= (sum(~isnan(unit1)) - sum(unit1 < cycle_times(1,1) | unit1 > cycle_times(end,2))); success = [success '>> FAILED']; end
    
    t2 = get_spikes(d, 'format', 'timing', 'cycleTimes', cycle_times, 'timeBase', 'time', 'name', 'unit1');
    
    % Check for correct data from get_spikes, 'timing'
    success = [success newline 'Pulled: count data'];
    if sum(sum(~isnan(t2))) ~= (sum(~isnan(unit1)) - sum(unit1 < cycle_times(1,1) | unit1 > cycle_times(end,2))); success = [success '>> FAILED']; end

    % Check for matching data from get_timing and get_spikes, 'timing'
    success = [success newline 'Matched: spike timing (time)'];
    if ~isequalwithequalnans(t1,t2); errs = [errs newline 'Matching: spike timing (time)']; end
    
    if with_plots
        figure();
        scatter(reshape(t1',1,[]), reshape(repmat([1:size(t1,1)]',1,size(t1,2))',1,[]), 'b.');
        hold on;
        scatter(reshape(t2',1,[]), reshape(repmat([1:size(t2,1)]',1,size(t2,2))',1,[]), 'b.');
        xlabel('Time (ms)');
        ylabel('Cycle Index');
        title('Comparison of spike timing (time)');
    end
    
    t1 = get_timing(d, cycle_times, 'timeBase', 'phase', 'name', 'unit1');
    t2 = get_spikes(d, 'format', 'timing', ...
        'cycleTimes', cycle_times, ...
        'timeBase', 'phase', ...
        'name', 'unit1');
    
    % Check for matching data from get_timing and get_spikes, 'timing',
    % 'phase'
    success = [success newline 'Matched: spike timing (phase)'];
    if ~isequalwithequalnans(t1,t2); errs = [errs newline 'Matching: spike timing (phase)']; end
    
    if with_plots
        figure();
        scatter(reshape(t1',1,[]), reshape(repmat([1:size(t1,1)]',1,size(t1,2))',1,[]), 'b.');
        hold on;
        scatter(reshape(t2',1,[]), reshape(repmat([1:size(t2,1)]',1,size(t2,2))',1,[]), 'b.');
        xlabel('Phase (rad)');
        ylabel('Cycle Index');
        title('Comparison of spike timing (phase)');
    end
catch e
    e
    global_errs{end+1} = {'Instantiating mi_data_neural'};
    disp([newline 'ERROR: Unable to instantiate mi_data_neural']);
end
    % PRINT RESULTS FROM CHECKS    
    disp(success);

    
%%

% In this script, we are just checking that the behavior class is
% functioning and that we can load cycleTimes into the class. For this
% reason, we arbitrarily set the d.data.cycleTimes.data = [unit1' unit1']. 

% try
%     disp([newline newline]);
%     clear d;
%     d = mi_data_behavior('test');
%     clear d;
% catch
%     global_errs{end+1} = {'Intantiating mi_data_behavior with ID only'};
%     disp([newline 'ERROR: Unable to instantiate mi_data_behavior with ID only']);
% end

try
    disp([newline newline]);
    clear d
    disp([newline '===== ===== ===== ===== =====']);
    disp(['RUNNING: mi_data_behavior()' newline newline]);
    d = mi_data_behavior('test', 'verbose', 5);

    add_cycleTimes(d, [unit1' unit1'], str_unit1, 30000);

    % CHECK OBJECT FOR INSTANTIATION CONSISTENCY
    success = (['----- ----- ----- ----- -----' newline 'SUCCESSFUL:' newline]);

    % Check for correct ID
    success = [success newline 'Assigned: ID'];
    if ~strcmp(d.ID, 'test'); success = [success '>> FAILED']; end
    
    % Check for correct Fs
    success = [success newline 'Assigned: Fs'];
    if ~d.Fs == 30000; success = [success '>> FAILED']; end
    
    % Check for correct verbose
    success = [success newline 'Assigned: verbose'];
    if ~d.verbose == 5; success = [success '>> FAILED']; end
    
     % Check for correct data name
    success = [success newline 'Assigned: data'];
    if ~isfield(d.data, 'cycleTimes'); success = [success '>> FAILED']; end

    % Check for correct data size
    success = [success newline 'Imported: data'];
    if any(size(d.data.cycleTimes.data) ~= [size(unit1,2), 2]); success = [success '>> FAILED']; end

    % Check for correct data info
    success = [success newline 'Imported: info'];
    if ~strcmp(d.data.cycleTimes.info,str_unit1); success = [success '>> FAILED']; end


    % Check for correct data from get_cycleTimes
    success = [success newline 'Pulled: cycleTimes'];
    if any(size(get_cycleTimes(d)) ~= [size(unit1,2) 2]); success = [success '>> FAILED']; end
    
    disp(success);
catch e
    e
    global_errs{end+1} = {'Instantiating mi_data_behavior with ID and verbose'};
    % Not possible to proceed without mi_data class
    
    disp([newline newline '===== ===== ===== ===== =====' newline 'GLOBAL ERRORS' newline]);
    for i=1:length(global_errs)
        disp(global_errs{i});
    end
    disp(['----- ----- ----- ----- -----' newline]);
    
    error('FATAL ERROR: Unable to construct mi_data_behavior objects');
end


%% CHECK mi_data_pressure: phase

% BRYCE:
%fnames = dir('D:\EMG_Data\chung\for_analysis\bl21lb21_20171218\bl21lb21_trial1_ch1_ch16\*.rhd');
%fnames = {fnames.name};
%fpath = 'D:\EMG_Data\chung\for_analysis\bl21lb21_20171218\bl21lb21_trial1_ch1_ch16';

% % RACHEL_lab:
% fnames = dir('C:\Users\RBARKE2\Projects\MergingCode\ContinuousMIEstimation\TestData\bl21lb21_trial1_ch1_ch16\*.rhd');
% fnames = {fnames.name};
% fpath = 'C:\Users\RBARKE2\Projects\MergingCode\ContinuousMIEstimation\TestData\bl21lb21_trial1_ch1_ch16';

% RACHEL_mac:
fnames = dir('/Users/Rachel/ContinuousMIEstimation/TestData/bl21lb21_trial1_ch1_ch16/*.rhd');
fnames = {fnames.name};
fpath = '/Users/Rachel/ContinuousMIEstimation/TestData/bl21lb21_trial1_ch1_ch16';


try
    disp([newline newline]);
    clear d
    d = mi_data_pressure('test', 'verbose', 5);
    add_cycleTimes(d, cycle_times, str_cycles, 30000);
    set_data_files(d, fnames, fpath);
    
    build_behavior(d);
    disp([newline '===== ===== ===== ===== =====']);
    disp(['RUNNING: mi_data_behavior()' newline newline]);
    disp(['--> mi_data_pressure()' newline newline]);
    disp(['--> --> phase' newline newline]);


    % CHECK OBJECT FOR INSTANTIATION CONSISTENCY
    success = (['----- ----- ----- ----- -----' newline 'SUCCESSFUL:' newline]);

    % Check for correct ID
    success = [success newline 'Assigned: ID'];
    if ~strcmp(d.ID, 'test'); success = [success '>> FAILED']; end
    
    % Check for correct Fs
    success = [success newline 'Assigned: Fs'];
    if ~d.Fs == 30000; success = [success '>> FAILED']; end
    
    % Check for correct verbose
    success = [success newline 'Assigned: verbose'];
    if ~d.verbose == 5; success = [success '>> FAILED']; end
    
     % Check for correct data name
    success = [success newline 'Assigned: data'];
    if ~isfield(d.data, 'cycleTimes'); success = [success '>> FAILED']; end

    % Check for correct data size
    success = [success newline 'Imported: data'];
    if any(size(d.data.cycleTimes.data) ~= size(cycle_times)); success = [success '>> FAILED']; end

    % Check for correct data info
    success = [success newline 'Imported: info'];
    if ~strcmp(d.data.cycleTimes.info,str_cycles); success = [success '>> FAILED']; end


    % Check for correct data from get_cycleTimes
    success = [success newline 'Pulled: cycleTimes'];
    if any(size(get_cycleTimes(d)) ~= size(cycle_times)); success = [success '>> FAILED']; end
    
    % Check for correct strFldr:
    success = [success newline 'Assigned: strFldr'];
    if ~strcmp(d.strFldr, fpath); succes= [success '>> FAILED']; end
    
    % Check for correct arrFiles
    success = [success newline 'Assigned: arrFiles'];
    if ~isequal(d.arrFiles, fnames); success = [success '>> FAILED' ]; end
    
%     % NOTE: TEMPORARIlY THIS IS NOT WORKING BECAUSE WE ARE NOT USING ALL
%     THE FILES CURRENTLY. ALSO, THIS MAY NEVER WORK DEPENDING ON HOW WE
%     DECIDE TO OMIT CYCLES AND HOW WE DOCUMENT THAT. 

%     % Check that size of rawBehav is correct
%     success = [success newline 'Pulled: rawBehav']; 
%     if ~any(size(d.data.cycleTimes.data) == size(find(~cellfun('isempty',obj.rawBehav)))); end

    % Plot some random cycles that were specified
    % Identify the cycles that were specified
    if with_plots
        cycles = find(~cellfun('isempty', d.rawBehav));
        cyclesToPlot = round(linspace(1,length(cycles),5));
        F1 = figure();
        colors = {[1 0 0], [0 1 0], [0 0 1], [0 1 1], [1 0 1]};
        for i = 1:length(cyclesToPlot)
            plot(d.rawBehav{cycles(cyclesToPlot(i)),1}, 'color', colors{i})
            hold on
        end
    end

    % Get behavior in many different ways
    % PHASE, RAW
    b1 = get_behavior(d, 'phase', 'raw', pi/2 , pi, 11);
    if with_plots
        F2 = figure();
        for i = 1:length(cyclesToPlot)
            plot(b1(cycles(cyclesToPlot(i)),:), 'color', colors{i})
            hold on
        end
        title('Transformed Sample Cycles: Phase, Raw')

        figure(F1);
        for i = 1:length(cyclesToPlot)
            len = length(d.rawBehav{cycles(cyclesToPlot(i))});
            idx1 = round(len*((pi/2)/(2*pi)));
            idx2 = idx1 + round(len*(pi/(2*pi)));
            idxs = round(linspace(idx1,idx2, 11));
            plot(idxs,b1(cycles(cyclesToPlot(i)), :), 'Marker', 'o','MarkerSize', 10 , 'MarkerFaceColor', colors{i}, 'LineStyle', 'none', 'MarkerEdgeColor', [0 0 0])
        end
        title('Raw Sample Cycles with projected transform: Phase, Raw')
    end
    
    % PHASE, Residual
    b2 = get_behavior(d, 'phase', 'residual', pi/2 , pi, 11);
    
    % Check for correct size behavior residuals
    success = [success newline 'Pulled: behavior (phase, residual)'];
    m = mean(b1,1, 'omitnan');
    if sum(sum(~isnan(b1))) ~= sum(sum(~isnan(b2)))
        success = [success '>> FAILED' ];    
    % Check that raw equals residual plus mean. 
    elseif ~isequaln(b1,(b2 + m))
        success = [success '>> FAILED' ];     
    end
    
    % Again, plot some random cycles that were specified
    % Identify the cycles that were specified
    if with_plots
        cycles = find(~cellfun('isempty', d.rawBehav));
        cyclesToPlot = round(linspace(1,length(cycles),5));
        F3 = figure();
        colors = {[1 0 0], [0 1 0], [0 0 1], [0 1 1], [1 0 1]};
        for i = 1:length(cyclesToPlot)
            plot(d.rawBehav{cycles(cyclesToPlot(i)),1}, 'color', colors{i})
            hold on
        end
    end
    
    if with_plots
        F4 = figure();
        for i = 1:length(cyclesToPlot)
            plot(b2(cycles(cyclesToPlot(i)),:), 'color', colors{i})
            hold on
        end
        title('Transformed Sample Cycles: Phase, Residual')

        figure(F3);
        for i = 1:length(cyclesToPlot)
            len = length(d.rawBehav{cycles(cyclesToPlot(i))});
            idx1 = round(len*((pi/2)/(2*pi)));
            idx2 = idx1 + round(len*(pi/(2*pi)));
            idxs = round(linspace(idx1,idx2, 11));
            plot(idxs,(b2(cycles(cyclesToPlot(i)), :) + m), 'Marker', 'o','MarkerSize', 10 , 'MarkerFaceColor', colors{i}, 'LineStyle', 'none', 'MarkerEdgeColor', [0 0 0])
        end
        title('Raw Sample Cycles with projected transform: Phase, Residual')        
        
    end
    
    % PHASE, PCA
    b1 = get_behavior(d, 'phase', 'raw', pi/2, pi, 600);
    b2 = get_behavior(d, 'phase', 'pca', pi/2 , pi, 600, 2);
    
    % Run PCA on b1 to compare with b2.
    m = mean(b1, 1, 'omitnan');
    [coeff, score, latent] = pca(b1);
    
    % Check that scores 1 and 2 match b2. 
    success = [success newline 'Pulled: behavior (Phase, PCA)'];
    if ~isequaln(score(:, 1:2), b2); success = [success '>> FAILED' ]; end
    
    % Check that the error on the full projections is low
    b1_Full_proj = score * coeff' + repmat(m, size(b1,1), 1);
    b1_proj_err = (b1 - b1_Full_proj)./b1;
    err = max(max(b1_proj_err));
    
    % Check that error between data and rpojections is less than .01%
    % RC20191202: NOTE: The error threshold here may need to be changed. 
    success = [success newline 'Matched: all PCs projection with raw'];
    if err >= 1e-4; success = [success '>> FAILED' ]; end
    
    
    % Plot the projections of the two PCs of interest with the raw data.    
    % Identify the cycles that were specified
    if with_plots
        cycles = find(~cellfun('isempty', d.rawBehav));
        cyclesToPlot = round(linspace(1,length(cycles),5));
        F5 = figure();
        colors = {[1 0 0], [0 1 0], [0 0 1], [0 1 1], [1 0 1]};
        for i = 1:length(cyclesToPlot)
            plot(d.rawBehav{cycles(cyclesToPlot(i)),1}, 'color', colors{i})
            hold on
        end
    end
    
    if with_plots
        F6 = figure();
        for i = 1:length(cyclesToPlot)
            plot(b2(cycles(cyclesToPlot(i)),:), 'color', colors{i})
            hold on
        end
        title('Transformed Sample Cycles: Phase, PCA')

        figure(F5);
        coeff_t = coeff';
        b2_proj = b2 * coeff_t(1:2,:);
        m = mean(b1,1,'omitnan');
        colors2 = {};
        for i = 1:length(cyclesToPlot)
            len = length(d.rawBehav{cycles(cyclesToPlot(i))});
            idx1 = round(len*((pi/2)/(2*pi)));
            idx2 = idx1 + round(len*(pi/(2*pi)));
            idxs = round(linspace(idx1,idx2, 600));
            colors2{i} = colors{i}*.7;
            plot(idxs,(b2_proj(cycles(cyclesToPlot(i)), :) + m),'LineWidth', 2 , 'LineStyle', ':', 'color', colors2{i})
        end
        title('Raw Sample Cycles with projected transform: Phase, PCA')        
        
    end
    disp(success);
catch e
    e
    global_errs{end+1} = {'Instantiating mi_data_behavior with ID and verbose'};
    % Not possible to proceed without mi_data class
    
    disp([newline newline '===== ===== ===== ===== =====' newline 'GLOBAL ERRORS' newline]);
    for i=1:length(global_errs)
        disp(global_errs{i});
    end
    disp(['----- ----- ----- ----- -----' newline]);
    
    error('FATAL ERROR: Unable to construct mi_data_behavior objects');
end

%% CHECK mi_data_pressure: time

% BRYCE:
%fnames = dir('D:\EMG_Data\chung\for_analysis\bl21lb21_20171218\bl21lb21_trial1_ch1_ch16\*.rhd');
%fnames = {fnames.name};
%fpath = 'D:\EMG_Data\chung\for_analysis\bl21lb21_20171218\bl21lb21_trial1_ch1_ch16';

% % RACHEL_lab:
% fnames = dir('C:\Users\RBARKE2\Projects\MergingCode\ContinuousMIEstimation\TestData\bl21lb21_trial1_ch1_ch16\*.rhd');
% fnames = {fnames.name};
% fpath = 'C:\Users\RBARKE2\Projects\MergingCode\ContinuousMIEstimation\TestData\bl21lb21_trial1_ch1_ch16';

%RACHEL_mac:
fnames = dir('/Users/Rachel/ContinuousMIEstimation/TestData/bl21lb21_trial1_ch1_ch16/*.rhd');
fnames = {fnames.name};
fpath = '/Users/Rachel/ContinuousMIEstimation/TestData/bl21lb21_trial1_ch1_ch16';


try
    disp([newline newline]);
    clear d
    d = mi_data_pressure('test', 'verbose', 5);
    add_cycleTimes(d, cycle_times, str_cycles, 30000);
    set_data_files(d, fnames, fpath);
    
    build_behavior(d);
    disp([newline '===== ===== ===== ===== =====']);
    disp(['RUNNING: mi_data_behavior()' newline newline]);
        disp(['--> mi_data_pressure()' newline newline]);
    disp(['--> --> time' newline newline]);


    % CHECK OBJECT FOR INSTANTIATION CONSISTENCY
    success = (['----- ----- ----- ----- -----' newline 'SUCCESSFUL:' newline]);

    % Check for correct ID
    success = [success newline 'Assigned: ID'];
    if ~strcmp(d.ID, 'test'); success = [success '>> FAILED']; end
    
    % Check for correct Fs
    success = [success newline 'Assigned: Fs'];
    if ~d.Fs == 30000; success = [success '>> FAILED']; end
    
    % Check for correct verbose
    success = [success newline 'Assigned: verbose'];
    if ~d.verbose == 5; success = [success '>> FAILED']; end
    
     % Check for correct data name
    success = [success newline 'Assigned: data'];
    if ~isfield(d.data, 'cycleTimes'); success = [success '>> FAILED']; end

    % Check for correct data size
    success = [success newline 'Imported: data'];
    if any(size(d.data.cycleTimes.data) ~= size(cycle_times)); success = [success '>> FAILED']; end

    % Check for correct data info
    success = [success newline 'Imported: info'];
    if ~strcmp(d.data.cycleTimes.info,str_cycles); success = [success '>> FAILED']; end


    % Check for correct data from get_cycleTimes
    success = [success newline 'Pulled: cycleTimes'];
    if any(size(get_cycleTimes(d)) ~= size(cycle_times)); success = [success '>> FAILED']; end
    
    % Check for correct strFldr:
    success = [success newline 'Assigned: strFldr'];
    if ~strcmp(d.strFldr, fpath); succes= [success '>> FAILED']; end
    
    % Check for correct arrFiles
    success = [success newline 'Assigned: arrFiles'];
    if ~isequal(d.arrFiles, fnames); success = [success '>> FAILED' ]; end
    
%     % NOTE: TEMPORARIlY THIS IS NOT WORKING BECAUSE WE ARE NOT USING ALL
%     THE FILES CURRENTLY. ALSO, THIS MAY NEVER WORK DEPENDING ON HOW WE
%     DECIDE TO OMIT CYCLES AND HOW WE DOCUMENT THAT. 

%     % Check that size of rawBehav is correct
%     success = [success newline 'Pulled: rawBehav']; 
%     if ~any(size(d.data.cycleTimes.data) == size(find(~cellfun('isempty',obj.rawBehav)))); end

    % Plot some random cycles that were specified
    % Identify the cycles that were specified
    if with_plots
        cycles = find(~cellfun('isempty', d.rawBehav));
        cyclesToPlot = round(linspace(1,length(cycles),5));
        F1 = figure();
        colors = {[1 0 0], [0 1 0], [0 0 1], [0 1 1], [1 0 1]};
        for i = 1:length(cyclesToPlot)
            plot(d.rawBehav{cycles(cyclesToPlot(i)),1}, 'color', colors{i})
            hold on
        end
    end

    % Get behavior in many different ways
    % TIME, RAW
    b1 = get_behavior(d, 'time', 'raw', 50 , 100, 11);
    if with_plots
        F2 = figure();
        for i = 1:length(cyclesToPlot)
            plot(b1(cycles(cyclesToPlot(i)),:), 'color', colors{i})
            hold on
        end
        title('Transformed Sample Cycles: Time, Raw')

        figure(F1);
        for i = 1:length(cyclesToPlot)
            idx1 = round((50/1000)*d.Fs);
            idx2 = idx1 + round((100/1000)*d.Fs);
            idxs = round(linspace(idx1,idx2, 11));
            plot(idxs,b1(cycles(cyclesToPlot(i)), :), 'Marker', 'o','MarkerSize', 10 , 'MarkerFaceColor', colors{i}, 'LineStyle', 'none', 'MarkerEdgeColor', [0 0 0])
        end
        title('Raw Sample Cycles with projected transform: Time, Raw')
    end
    
    % PHASE, Residual
    b2 = get_behavior(d, 'time', 'residual', 50 , 100, 11);
    
    % Check for correct size behavior residuals
    success = [success newline 'Pulled: behavior (time, residual)'];
    m = mean(b1,1, 'omitnan');
    if sum(sum(~isnan(b1))) ~= sum(sum(~isnan(b2)))
        success = [success '>> FAILED' ];    
    % Check that raw equals residual plus mean. 
    elseif ~isequaln(b1,(b2 + m))
        success = [success '>> FAILED' ];     
    end
    
    % Again, plot some random cycles that were specified
    % Identify the cycles that were specified
    if with_plots
        cycles = find(~cellfun('isempty', d.rawBehav));
        cyclesToPlot = round(linspace(1,length(cycles),5));
        F3 = figure();
        colors = {[1 0 0], [0 1 0], [0 0 1], [0 1 1], [1 0 1]};
        for i = 1:length(cyclesToPlot)
            plot(d.rawBehav{cycles(cyclesToPlot(i)),1}, 'color', colors{i})
            hold on
        end
    end
    
    if with_plots
        F4 = figure();
        for i = 1:length(cyclesToPlot)
            plot(b2(cycles(cyclesToPlot(i)),:), 'color', colors{i})
            hold on
        end
        title('Transformed Sample Cycles: Time, Residual')

        figure(F3);
        for i = 1:length(cyclesToPlot)
            idx1 = round((50/1000)*d.Fs);
            idx2 = idx1 + round((100/1000)*d.Fs);
            idxs = round(linspace(idx1,idx2, 11));
            plot(idxs,(b2(cycles(cyclesToPlot(i)), :) + m), 'Marker', 'o','MarkerSize', 10 , 'MarkerFaceColor', colors{i}, 'LineStyle', 'none', 'MarkerEdgeColor', [0 0 0])
        end
        title('Raw Sample Cycles with projected transform: Time, Residual')        
        
    end
    
    % TIME, PCA
    b1 = get_behavior(d, 'time', 'raw', 50, 100, 200);
    b2 = get_behavior(d, 'time', 'pca', 50 , 100, 200, 2);
    
    % Run PCA on b1 to compare with b2.
    m = mean(b1, 1, 'omitnan');
    [coeff, score, latent] = pca(b1);
    
    % Check that scores 1 and 2 match b2. 
    success = [success newline 'Pulled: behavior (time, PCA)'];
    if ~isequaln(score(:, 1:2), b2); success = [success '>> FAILED' ]; end
    
    % Check that the error on the full projections is low
    b1_Full_proj = score * coeff' + repmat(m, size(b1,1), 1);
    b1_proj_err = (b1 - b1_Full_proj)./b1;
    err = max(max(b1_proj_err));
    
    % Check that error between data and rpojections is less than .01%
    % RC20191202: NOTE: The error threshold here may need to be changed. 
    success = [success newline 'Matched: all PCs projection with raw'];
    if err >= 1e-4; success = [success '>> FAILED' ]; end
    
    
    % Plot the projections of the two PCs of interest with the raw data.    
    % Identify the cycles that were specified
    if with_plots
        cycles = find(~cellfun('isempty', d.rawBehav));
        cyclesToPlot = round(linspace(1,length(cycles),5));
        F5 = figure();
        colors = {[1 0 0], [0 1 0], [0 0 1], [0 1 1], [1 0 1]};
        for i = 1:length(cyclesToPlot)
            plot(d.rawBehav{cycles(cyclesToPlot(i)),1}, 'color', colors{i})
            hold on
        end
    end
    
    if with_plots
        F6 = figure();
        for i = 1:length(cyclesToPlot)
            plot(b2(cycles(cyclesToPlot(i)),:), 'color', colors{i})
            hold on
        end
        title('Transformed Sample Cycles: Time, PCA')

        figure(F5);
        coeff_t = coeff';
        b2_proj = b2 * coeff_t(1:2,:);
        m = mean(b1,1,'omitnan');
        colors2 = {};
        for i = 1:length(cyclesToPlot)
            idx1 = round((50/1000)*d.Fs);
            idx2 = idx1 + round((100/1000)*d.Fs);
            idxs = round(linspace(idx1,idx2, 200));
            colors2{i} = colors{i}*.7;
            plot(idxs,(b2_proj(cycles(cyclesToPlot(i)), :) + m),'LineWidth', 2 , 'LineStyle', ':', 'color', colors2{i})
        end
        title('Raw Sample Cycles with projected transform: Time, PCA')        
        
    end
    disp(success);
catch e
    e
    global_errs{end+1} = {'Instantiating mi_data_behavior with ID and verbose'};
    % Not possible to proceed without mi_data class
    
    disp([newline newline '===== ===== ===== ===== =====' newline 'GLOBAL ERRORS' newline]);
    for i=1:length(global_errs)
        disp(global_errs{i});
    end
    disp(['----- ----- ----- ----- -----' newline]);
    
    error('FATAL ERROR: Unable to construct mi_data_behavior objects');
end

%%  mi_analysis
try
    clear d
    clear a
    disp([newline '===== ===== ===== ===== =====']);
    disp(['RUNNING: mi_analysis()' newline newline]);

    d = mi_data_neural('test', 'verbose', 5);

    add_spikes(d, unit1, str_unit1, 30000, 'unit1');
    add_spikes(d, unit2, str_unit2, 30000, 'unit2');
    
    b = mi_data_pressure('test', 'verbose', 5);
    add_cycleTimes(b, cycle_times, str_cycles, 30000);
    
    % Construct mi_analysis object
    a = mi_analysis(d, b, {'unit1' , 'unit2'}, 'verbose', 5);
    
    
    % CHECK OBJECT FOR INSTANTIATION CONSISTENCY
    success = (['----- ----- ----- ----- -----' newline 'SUCCESSFUL:' newline]);
    
    % Check for correct objData
    success = [success newline 'Assigned: objData'];
    if ~isa(a.objData,'mi_data_neural'); success = [success '>> FAILED']; end
    
    % Check for correct varNames
    success = [success newline 'Assigned: varNames'];
    if ~isequal(a.varNames, {'unit1', 'unit2'}); success = [success '>> FAILED']; end
    
    % Check for verbose
    success = [success newline 'Assigned: verbose'];
    if a.verbose ~= 5; success = [success '>> FAILED']; end
    
    % Check for sim manager object
    success = [success newline 'Constructed: sim_manager'];
    if ~isa(a.sim_manager,'mi_ksg_sims'); success = [success '>> FAILED']; end
    
    % Check for integration with data object
    success = [success newline 'Matched: varNames to objData.data'];
    if ~isfield(a.objData.data, a.varNames{1}) || ~isfield(a.objData.data, a.varNames{2}); success = [success '>> FAILED']; end
    
    disp(success)
    
catch e
    e
    global_errs{end+1} = {'Instantiating mi_analysis'};
    % Not possible to proceed without mi_analysis class
    
    disp([newline newline '===== ===== ===== ===== =====' newline 'GLOBAL ERRORS' newline]);
    for i=1:length(global_errs)
        disp(global_errs{i});
    end
    disp(['----- ----- ----- ----- -----' newline]);
    
    error('FATAL ERROR: Unable to construct mi_analysis object');
end

%%  mi_analysis: calc_count_count
try
    clear d
    clear a
    disp([newline '===== ===== ===== ===== =====']);
    disp(['RUNNING: mi_analysis(): count_count' newline newline]);

    d = mi_data_neural('test', 'verbose', 5);

    add_spikes(d, unit1, str_unit1, 30000, 'unit1');
    add_spikes(d, unit2, str_unit2, 30000, 'unit2');
    
    b = mi_data_pressure('test', 'verbose', 5);
    add_cycleTimes(b, cycle_times, str_cycles, 30000);
    
    % Construct mi_analysis object
    a = calc_count_count(d, b, {'unit1' , 'unit2'}, 'verbose', 5);
    
    
    % CHECK OBJECT FOR INSTANTIATION CONSISTENCY
    success = (['----- ----- ----- ----- -----' newline 'SUCCESSFUL:' newline]);
    
    % Check for correct objData
    success = [success newline 'Assigned: objData'];
    if ~isa(a.objData,'mi_data_neural'); success = [success '>> FAILED']; end
    
    % Check for correct objBehav
    success = [success newline 'Assigned: objBehav'];
    if ~isa(a.objBehav,'mi_data_pressure'); success = [success '>> FAILED']; end
    
    % Check for correct varNames
    success = [success newline 'Assigned: varNames'];
    if ~isequal(a.varNames, {'unit1', 'unit2'}); success = [success '>> FAILED']; end
    
    % Check for verbose
    success = [success newline 'Assigned: verbose'];
    if a.verbose ~= 5; success = [success '>> FAILED']; end
    
    % Check for sim manager object
    success = [success newline 'Constructed: sim_manager'];
    if ~isa(a.sim_manager,'mi_ksg_sims'); success = [success '>> FAILED']; end
    
    % Check for integration with data object
    success = [success newline 'Matched: varNames to objData.data'];
    if ~isfield(a.objData.data, a.varNames{1}) || ~isfield(a.objData.data, a.varNames{2}); success = [success '>> FAILED']; end
    
    % Run buildMIs()
    a.buildMIs();
    

    success = [success newline 'Assigned: Subgroup IDs'];
    if sum(compVal) ~= 0; success = [success '>> FAILED']; end
    
    disp(success)
    
catch e
    e
    global_errs{end+1} = {'Instantiating mi_analysis'};
    % Not possible to proceed without mi_analysis class
    
    disp([newline newline '===== ===== ===== ===== =====' newline 'GLOBAL ERRORS' newline]);
    for i=1:length(global_errs)
        disp(global_errs{i});
    end
    disp(['----- ----- ----- ----- -----' newline]);
    
    error('FATAL ERROR: Unable to construct mi_analysis object');
end

%%  mi_analysis: calc_timing_count
try
    clear d
    clear a
    disp([newline '===== ===== ===== ===== =====']);
    disp(['RUNNING: mi_analysis(): timing_count' newline newline]);


    d = mi_data_neural('test', 'verbose', 5);

    add_spikes(d, unit1, str_unit1, 30000, 'unit1');
    add_spikes(d, unit2, str_unit2, 30000, 'unit2');
    
    b = mi_data_pressure('test', 'verbose', 5);
    add_cycleTimes(b, cycle_times, str_cycles, 30000);
    
    % Construct mi_analysis object
    a = calc_timing_count(d, b, {'unit1' , 'unit2'}, 'verbose', 5);
    
    
    % CHECK OBJECT FOR INSTANTIATION CONSISTENCY
    success = (['----- ----- ----- ----- -----' newline 'SUCCESSFUL:' newline]);
    
    % Check for correct objData
    success = [success newline 'Assigned: objData'];
    if ~isa(a.objData,'mi_data_neural'); success = [success '>> FAILED']; end
    
    % Check for correct objBehav
    success = [success newline 'Assigned: objBehav'];
    if ~isa(a.objBehav,'mi_data_pressure'); success = [success '>> FAILED']; end
    
    % Check for correct varNames
    success = [success newline 'Assigned: varNames'];
    if ~isequal(a.varNames, {'unit1', 'unit2'}); success = [success '>> FAILED']; end
    
    % Check for verbose
    success = [success newline 'Assigned: verbose'];
    if a.verbose ~= 5; success = [success '>> FAILED']; end
    
    % Check for correct timebase (specific to timing subclass)
    success = [success newline 'Assigned: timebase'];
    if ~isequal(a.timebase, 'time'); success = [success '>> FAILED']; end
    
    % Check for sim manager object
    success = [success newline 'Constructed: sim_manager'];
    if ~isa(a.sim_manager,'mi_ksg_sims'); success = [success '>> FAILED']; end
    
    % Check for integration with data object
    success = [success newline 'Matched: varNames to objData.data'];
    if ~isfield(a.objData.data, a.varNames{1}) || ~isfield(a.objData.data, a.varNames{2}); success = [success '>> FAILED']; end
    
    % Run buildMIs()
    a.buildMIs();
    

    
    disp(success)
    
catch e
    e
    global_errs{end+1} = {'Instantiating mi_analysis'};
    % Not possible to proceed without mi_analysis class
    
    disp([newline newline '===== ===== ===== ===== =====' newline 'GLOBAL ERRORS' newline]);
    for i=1:length(global_errs)
        disp(global_errs{i});
    end
    disp(['----- ----- ----- ----- -----' newline]);
    
    error('FATAL ERROR: Unable to construct mi_analysis object');
end

%%  mi_analysis: calc_count_behav
try
    clear d
    clear b
    clear a
    disp([newline '===== ===== ===== ===== =====']);
    disp(['RUNNING: mi_analysis(): count_behav' newline newline]);


    d = mi_data_neural('test', 'verbose', 5);

    add_spikes(d, unit1, str_unit1, 30000, 'unit1');
    
    b = mi_data_pressure('test', 'verbose', 5);
    add_cycleTimes(b, cycle_times, str_cycles, 30000);
    
    % Get behavior for pressure class
    set_data_files(b, fnames, fpath);
    
    build_behavior(b);
    
    % Construct mi_analysis object
    a = calc_count_behav(d, b, {'unit1'}, 'verbose', 5);
    
    
    % CHECK OBJECT FOR INSTANTIATION CONSISTENCY
    success = (['----- ----- ----- ----- -----' newline 'SUCCESSFUL:' newline]);
    
    % Check for correct objData
    success = [success newline 'Assigned: objData'];
    if ~isa(a.objData,'mi_data_neural'); success = [success '>> FAILED']; end
    
    % Check for correct objBehav
    success = [success newline 'Assigned: objBehav'];
    if ~isa(a.objBehav,'mi_data_pressure'); success = [success '>> FAILED']; end
    
    % Check for correct varNames
    success = [success newline 'Assigned: varNames'];
    if ~isequal(a.varNames, {'unit1'}); success = [success '>> FAILED']; end
    
    % Check for verbose
    success = [success newline 'Assigned: verbose'];
    if a.verbose ~= 5; success = [success '>> FAILED']; end
    
    % Check for correct timebase (specific to timing subclass)
    success = [success newline 'Assigned: timebase'];
    if ~isequal(a.b_timeBase, 'phase'); success = [success '>> FAILED']; end
    
    % Check for sim manager object
    success = [success newline 'Constructed: sim_manager'];
    if ~isa(a.sim_manager,'mi_ksg_sims'); success = [success '>> FAILED']; end
    
    % Check for integration with data object
    success = [success newline 'Matched: varNames to objData.data'];
    if ~isfield(a.objData.data, a.varNames{1}); success = [success '>> FAILED']; end
    
    % Run buildMIs()
    %a.buildMIs();
    

    
    disp(success)
    
catch e
    e
    global_errs{end+1} = {'Instantiating mi_analysis'};
    % Not possible to proceed without mi_analysis class
    
    disp([newline newline '===== ===== ===== ===== =====' newline 'GLOBAL ERRORS' newline]);
    for i=1:length(global_errs)
        disp(global_errs{i});
    end
    disp(['----- ----- ----- ----- -----' newline]);
    
    error('FATAL ERROR: Unable to construct mi_analysis object');
end
%%
diary off