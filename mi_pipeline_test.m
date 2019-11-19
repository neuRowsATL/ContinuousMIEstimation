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

% diary mi_pipeline_test_diary.txt
% diary on

%% RUN MI_DATA
load('TestData/20191018_bl21lb21_171218_spikes.mat');
unit1 = unit1*1000.;
unit2 = unit2*1000.;
unit3 = unit3*1000.;

str_unit1 = 'TestData/20191018_bl21lb21_171218_spikes.mat/unit1';
str_unit2 = 'TestData/20191018_bl21lb21_171218_spikes.mat/unit2';
str_unit3 = 'TestData/20191018_bl21lb21_171218_spikes.mat/unit3';

load('TestData/20191018_bl21lb21_171218_cycles.mat');
cycle_times = [cycle_times(1:end-1)' cycle_times(2:end)'];

str_cycles = 'TestData/20191018_bl21lb21_171218_cycles.mat';

%%
try
    disp([newline newline]);
    clear d;
    d = mi_data('test');
    
    
    % CHECK OBJECT FOR INSTANTIATION CONSISTENCY
    errs = ['----- ----- ----- ----- -----' newline 'ERRORS:' newline];
    disp([newline newline '===== ===== ===== ===== ===== ' newline 'SUCCESSFUL:' newline]);
    
    add_data(d, unit1, str_unit1, 30000);
    
    if all(size(d.data.noname.data) == size(unit1))
        disp('Imported: data');
    else
        errs = [errs newline 'Importing: data dims do not match for single dataset'];
    end
    
    if all(size(get_data(d)) == size(unit1))
        disp('Pulled: data');
    else
        errs = [errs newline 'Pulling: get_data() dims do not match for single dataset'];
    end
    
    disp(errs);
    clear d;
catch e
    e
    global_errs{end+1} = {'Intantiating mi_data with ID only'};
    disp([newline 'ERROR: Unable to instantiate mi_data with ID only']);
end

try
    disp([newline newline]);
    clear d;
    d = mi_data('test', 'verbose', 5);

    add_data(d, unit1, str_unit1, 30000, 'unit1');
    add_data(d, unit2, str_unit2, 30000, 'unit2');
    add_data(d, unit3, str_unit3, 30000, 'unit3');


    % CHECK OBJECT FOR INSTANTIATION CONSISTENCY
    errs = ['----- ----- ----- ----- -----' newline 'ERRORS:' newline];
    disp([newline newline '===== ===== ===== ===== ===== ' newline 'SUCCESSFUL:' newline]);
    if strcmp(d.ID, 'test')
        disp('Assigned: ID');
    else
        errs = [errs newline 'Assigning: ID'];
    end
    if d.Fs == 30000
        disp('Assigned: Fs');
    else
        errs = [errs newline 'Assigning: Fs'];
    end
    if d.verbose == 5
        disp('Assigned: verbose');
    else
        errs = [errs newline 'Assigning: verbose'];
    end
    if isfield(d.data, 'unit1') && isfield(d.data, 'unit2') && isfield(d.data, 'unit3')
        disp('assigned: data');
    else
        errs = [errs newline 'Assigning: data'];
    end

    if all(size(d.data.unit1.data) == size(unit1)) && ...
            all(size(d.data.unit2.data) == size(unit2)) && ...
            all(size(d.data.unit3.data) == size(unit3))
        disp('Imported: data');
    else
        errs = [errs newline 'Importing: data dims do not match'];
    end

    if strcmp(d.data.unit1.info,str_unit1) && ...
            strcmp(d.data.unit2.info,str_unit2) && ...
            strcmp(d.data.unit3.info,str_unit3)
        disp('Imported: info');
    else
        errs = [errs newline 'Importing: data info'];
    end

    if all(size(get_data(d,'unit1')) == size(unit1)) && ...
            all(size(get_data(d,'unit2')) == size(unit2)) && ...
            all(size(get_data(d,'unit3')) == size(unit3))
        disp('Pulled: data');
    else
        errs = [errs newline 'Pulling: get_data() dims do not match'];
    end
    
    disp(errs);
catch
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
clear d

d = mi_data_neural('test');

add_spikes(d, unit1, str_unit1, 30000);

get_spikes(d, 'format', 'raw');

get_count(d, cycle_times);
get_spikes(d, 'format', 'count', 'cycleTimes', cycle_times);

get_timing(d, cycle_times, 'timeBase', 'time');
get_spikes(d, 'format', 'timing', 'cycleTimes', cycle_times, 'timeBase', 'time');


%%

try
    disp([newline newline]);
    clear d;
    d = mi_data_behavior('test');
    clear d;
catch
    global_errs{end+1} = {'Intantiating mi_data_behavior with ID only'};
    disp([newline 'ERROR: Unable to instantiate mi_data_behavior with ID only']);
end

try
    disp([newline newline]);
    clear d
    d = mi_data_behavior('test', 'verbose', 5);

    add_cycleTimes(d, [unit1' unit1'], str_unit1, 30000);

    % CHECK OBJECT FOR INSTANTIATION CONSISTENCY
    errs = ['----- ----- ----- ----- -----' newline 'ERRORS:' newline];
    disp([newline newline '===== ===== ===== ===== ===== ' newline 'SUCCESSFUL:' newline]);
    if strcmp(d.ID, 'test')
        disp('Assigned: ID');
    else
        errs = [errs newline 'Assigning: ID'];
    end
    if d.Fs == 30000
        disp('Assigned: Fs');
    else
        errs = [errs newline 'Assigning: Fs'];
    end
    if d.verbose == 5
        disp('Assigned: verbose');
    else
        errs = [errs newline 'Assigning: verbose'];
    end
    if isfield(d.data, 'cycleTimes')
        disp('assigned: cycleTimes');
    else
        errs = [errs newline 'Assigning: cycleTimes'];
    end

    if all(size(d.data.cycleTimes.data) == [size(unit1,2) 2])
        disp('Imported: cycleTimes');
    else
        errs = [errs newline 'Importing: cycleTimes dims do not match'];
    end

    if strcmp(d.data.cycleTimes.info,str_unit1)
        disp('Imported: info');
    else
        errs = [errs newline 'Importing: cycleTimes info'];
    end

    if all(size(get_cycleTimes(d)) == [size(unit1,2) 2])
        disp('Pulled: cycleTimes');
    else
        errs = [errs newline 'Pulling: get_cycleTimes() dims do not match'];
    end
    
    disp(errs);
catch
    global_errs{end+1} = {'Instantiating mi_data_behavior with ID and verbose'};
    % Not possible to proceed without mi_data class
    
    disp([newline newline '===== ===== ===== ===== =====' newline 'GLOBAL ERRORS' newline]);
    for i=1:length(global_errs)
        disp(global_errs{i});
    end
    disp(['----- ----- ----- ----- -----' newline]);
    
    error('FATAL ERROR: Unable to construct mi_data_behavior objects');
end


%%
diary off