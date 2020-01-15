classdef mi_ksg_sims < handle
    % MI_KSG_sim_manager provides a structure within which mutual
    % information calculations can be run in parallel. MI_core objects are
    % added to a list and sim_manager then sets up a "master list" of data
    % and parameters to run in parallel. MI calculations are returned to
    % respective mi_core objects upon completion.
    
    properties
        verbose % set level of output for reference and debugging
        notes = {}
        
        mi_core_arr = {} % array of mi_core objects
        
        thread_count % number of available cores for parallel
        par_mode % flag to run in parallel mode
    end
    methods
        function obj = mi_ksg_sims(mode, verbose)
            if nargin == 1
                % initialize parallel mode and default to no output
                obj.par_mode = mode;
                obj.verbose = 0;
            elseif nargin == 2
                % initialize parallel mode and output level
                obj.par_mode = mode;
                obj.verbose = verbose;
            else
                % default to parallel mode and no output
                obj.par_mode = true;
                obj.verbose = 0;
            end
            
            if obj.verbose > 0; disp('Initializing MI_KSG_sim_manager...'); end
            
            % determine number of available cores
            c = parcluster('local');
            obj.thread_count = c.NumWorkers;
            
            if obj.verbose > 1; disp(obj); end
            
            % ===== ===== ===== ===== =====
            % Should a line be added to automatically determine whether it
            % will be more efficient to run the simulations in serial or in
            % parallel?
            % Also add a way to change the number of cores to use?
            
        end

        function add_sim(obj, obj_mi_core)
            % add a MI core object to list of simulations
            if obj.verbose > 0; disp('Adding MI_KSG_core to sims...'); end
            
            tmp_core_arr = obj.mi_core_arr; % create temp obj for convenience
            if isempty(tmp_core_arr) % if no core objs already associated
                key = num2str(dec2hex(round(rand(1)*100000)));
            else
                while 1 % make sure core objs have unique index
                    key = num2str(dec2hex(round(rand(1)*100000)));
                    if ~any(strcmp(tmp_core_arr(:,2), key))
                        break;
                    end
                end
            end
            obj.mi_core_arr = cat(1, tmp_core_arr, {obj_mi_core key}); % add core obj
            if obj.verbose > 0; disp('--> Done'); end
        end
        
        function remove_sim(obj, idx)            
            % remove a simulation from the list of sims
            obj.mi_core_arr(idx) = {};
        end
        
        function run_sims(obj)
            % run simulations in parallel or in serial
            
            if obj.verbose > 0; disp('Running MI_KSG_sim_manager simulations...'); end
            % retrive simulation data and parameters from each mi_core object
            % populate single cell matrix of all data/params
            if obj.verbose > 1; disp('>> Setting up simulations...'); end
            sim_set = cell(0,5);
            for i=1:size(obj.mi_core_arr,1)
                [tmp_core, tmp_key] = obj.mi_core_arr{i,:};
                tmp_set = get_core_dataset(tmp_core); % get MI data and parameters from core obj
                tmp_set(:,5) = {tmp_key}; % add core obj identifier for each data/param set
                sim_set = cat(1, sim_set, tmp_set); % add data/param set with identifiers to sim set
            end
                        
           % run MI calculations
           if obj.verbose > 1; disp('>> Running simulations...'); end
           sim_data = cell(size(sim_set,1),4); % pre-allocate memory
           if obj.par_mode > 0
               parfor i=1:size(sim_set,1) % run simulations in parallel
                   tmp_sim_set = sim_set(i,:); % needed for parfor
                   if length(tmp_sim_set{1}) < 5
                       MI = NaN;
                       if obj.verbose > 2; disp('MI = NaN'); end
                   elseif length(tmp_sim_set{1}) < tmp_sim_set{3}
                       MI = NaN;
                       if obj.verbose > 2; disp('MI = NaN'); end
                   else
                       MI = MIxnyn(tmp_sim_set{1}, tmp_sim_set{2}, tmp_sim_set{3}); % run MI calculation
                       if obj.verbose > 2; disp('MI estimated'); end
                   end
                   sim_data(i,:) = {MI/log(2) tmp_sim_set{3} tmp_sim_set{4} tmp_sim_set{5}}; % add results with params/index to data set
               end
           else
               for i=1:size(sim_set,1)
                   if obj.verbose > 2; disp(['  > Sim ' num2str(i)]); end
                   tmp_sim_set = sim_set(i,:); % for convenience
                   if length(tmp_sim_set{1}) < 5
                       MI = NaN;
                       if obj.verbose > 2; disp('MI = NaN'); end
                   elseif length(tmp_sim_set{1}) < tmp_sim_set{3}
                       MI = NaN;
                       if obj.verbose > 2; disp('MI = NaN'); end
                   else
                       MI = MIxnyn(tmp_sim_set{1}, tmp_sim_set{2}, tmp_sim_set{3}); % run MI calculation
                       if obj.verbose > 2; disp('MI estimated'); end
                   end
                   sim_data(i,:) = {MI/log(2) tmp_sim_set{3} tmp_sim_set{4} tmp_sim_set{5}}; % add results with params/index to data set
               end
           end
            
            % return MI calculations to respective mi_core objects
            if obj.verbose > 1; disp('>> Returning data to MI cores...'); end
            core_keys = unique(sim_data(:,4));
            for key_ix = 1:length(core_keys)
                data_ixs = find(strcmp(sim_data(:,4), core_keys{key_ix}) == 1); % find data entries that belong to core obj
                core_ix = find(strcmp([obj.mi_core_arr(:,2)], core_keys(key_ix)) == 1); % find core obj in list of mi_core
                set_core_data(obj.mi_core_arr{core_ix}, sim_data(data_ixs,1:3)); % send data back to respective mi_core objs
                if obj.mi_core_arr{core_ix}.opt_k == 0
                    find_k_value(obj.mi_core_arr{core_ix}); % optimize k-value if opt_k == 0
                end
            end
        end
        
    end
    
end
