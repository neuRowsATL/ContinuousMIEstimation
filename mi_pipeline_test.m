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

try
    disp([newline newline]);
    d = mi_data('test');
    clear d;
catch
    global_errs{end+1} = {'Intantiating mi_data with ID only'};
    disp([newline 'ERROR: Unable to instantiate mi_data with ID only']);
end

try
    disp([newline newline]);
    d = mi_data('test', 'verbose', 5);

    str_unit1 = 'TestData/20191018_bl21lb21_171218_spikes.mat/unit1';
    str_unit2 = 'TestData/20191018_bl21lb21_171218_spikes.mat/unit2';
    str_unit3 = 'TestData/20191018_bl21lb21_171218_spikes.mat/unit3';

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




diary off