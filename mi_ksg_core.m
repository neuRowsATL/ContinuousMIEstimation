classdef mi_ksg_core < handle
    % MI_KSG_core is used to set up sets of simulations to determine an
    % optimum k-value for a mutual information estimate and also calculates
    % mutual information and estimation error.
    
    properties
        verbose % for debugging purposes
        
        x % matrix of x data
        
        y % matrix of y data
        
        k_values % array of k-values
        
        mi_data % MI value, error estimate, data fraction, k-value
        
        opt_k % optimized k value; if -1, only runs MI calculation without any error estimate
        
        data_fracs = 10 % number of data fractions
        
        append % Taken from analysis object. 
               % Specify whether to re-run all analysis or to just run analysis for k-values that have not been previously included. 
        
        sim_obj % sim_manager object
    end
    
    methods
        function obj = mi_ksg_core(sim_obj, x, y, varargin)
            % This function generates the core object with the x and y
            % values, desired k values to run, and whether to estimate the error. 
            
            % Instantiate input parser
            p = inputParser;
            
            % Set up required inputs
            validate_sim_obj = @(x_var) assert(isa(x_var, 'mi_ksg_sims'), 'sim_obj must be a valid sim object');
            p.addRequired('sim_obj', validate_sim_obj);
            
            validate_x = @(x_var) assert(ismatrix(x_var), 'x must be a matrix');
            p.addRequired('x', validate_x);
            
            validate_y = @(x_var) assert(ismatrix(x_var), 'y must be a matrix');
            p.addRequired('y', validate_y);
            
            % Add optional
            
            % ks_arr
            default_ks_arr = 1:9;
            validate_ks_arr = @(x_var) assert(ismatrix(x_var), 'ks_arr must be a vector of integers');
            p.addParameter('ks_arr', default_ks_arr, validate_ks_arr);
            
            % opt_k
            valid_opt_k = [0 1 -1];
            default_opt_k = 1;
            validate_opt_k = @(x_var) assert(ismember(x_var, valid_opt_k), 'opt_k must be 0, 1 or -1');
            p.addParameter('opt_k', default_opt_k, validate_opt_k);
            
            % append
            default_append= true;
            validate_append = @(x_var) assert(islogical(x_var), 'append must be a logical value');
            p.addParameter('append', default_append, validate_append);
            
            % verbose
            default_verbose = 1; 
            validate_verbose = @(x_var) assert(isnumeric(x_var) && (rem(x_var,1) == 0), 'verbose must be an integer');
            p.addParameter('verbose', default_verbose, validate_verbose);
            
            % Parse the inputs
            p.KeepUnmatched = 1;
            p.parse(sim_obj, x, y, varargin{:});
            
            obj.x = p.Results.x;
            obj.y = p.Results.y;            
            obj.sim_obj = p.Results.sim_obj;
            obj.k_values = p.Results.ks_arr;
            obj.opt_k = p.Results.opt_k;
            obj.verbose = p.Results.verbose;
            
            add_sim(sim_obj, obj); % add this core obj to sim_manager list
        end
        
        % These methods are used to interface with other classes for data
        % analysis and visualization        
        function r = get_core_dataset(obj)
            
            
            % get cell array of data for MI calculation
            r = cell(0,4);
            if obj.opt_k < 0
                % only run MI calculation without error estimate
                for i=1:length(obj.k_values)
                    while 1
                        % generate unique key to track each simulation
                        key = num2str(dec2hex(round(rand(1)*100000)));
                        if ~any(strcmp(r(:,4), key))
                            break;
                        end
                    end
                    r = cat(1, r, {obj.x obj.y obj.k_values(i) key});
                end
            else
                % run MI calculation with error estimates
                 if obj.append
                    if isempty(obj.mi_data)
                        % Run estimates for all k values if none have been
                        % run yet
                        for i = 1:length(obj.k_values)
                            % create datasets for data fractions with unique key
                            % to track each simulation
                            r = cat(1, r, fractionate_data(obj, obj.k_values(i)));
                        end
                    else
                        mi_data = cell2mat(obj.mi_data);
                        k_finished = unique(mi_data(:,4));
                        for i = 1:length(obj.k_values)
                            if ismember(obj.k_values(i), k_finished)
                                continue
                            else
                                % create datasets for data fractions with unique key
                                % to track each simulation
                                r = cat(1, r, fractionate_data(obj, obj.k_values(i)));
                            end
                        end
                    end
                else
                    for i=1:length(obj.k_values)
                        % create datasets for data fractions with unique key
                        % to track each simulation
                        r = cat(1, r, fractionate_data(obj, obj.k_values(i)));
                    end
                end
            end 
        end
        
        function set_core_data(obj, dataset)
            % take sim_manager MI calculations and process results
            
            data_keys = unique([dataset(:,3)]); % extract simulation keys
            tmp_mi_data = cell(0,4);
            for key_ix = 1:length(data_keys) % iterate through each MI error estimation set
                tmp_match = strcmp([dataset(:,3)], data_keys(key_ix)); % find MI calculations that correspond to same data fractions
                count = sum(tmp_match); % determine number of data fractions
                data_ixs = find(tmp_match == 1); % identify which simulations to include
                
                mi = [dataset{data_ixs,1}];
                k = dataset{data_ixs(1),2};
                
                tmp_mi_data = cat(1, tmp_mi_data, {mean(mi) var(mi) count k}); % append MI with error estimation
            end
            
            if obj.append
                tmp_mi_data = cat(1, tmp_mi_data, obj.mi_data);
            end
            
            obj.mi_data = sortrows(tmp_mi_data,[4,3]);
        end
        
        function r = get_mi(obj, k, errThreshold)
            % get mutual information and error estimates
            data_ixs = cell2mat(obj.mi_data(:,4)) == k; % find MI calcs with k-value
            
            % calculate estimated error
            listSplitSizes = cell2mat(obj.mi_data(data_ixs,3));
            MIs = cell2mat(obj.mi_data(data_ixs,1));
            listVariances = cell2mat(obj.mi_data(data_ixs,2));
            listVariances = listVariances(2:end);
            
            k = listSplitSizes(2:end);
            variancePredicted = sum((k-1)./k.*listVariances)./sum((k-1));
            

            
% -------------THIS CALCULATES THE ERROR OF THE ERROR-----------------------------------            
%             N = size(obj.x,2);
%             Sml=variancePredicted*N;
%             varS = 2*Sml^2/sum((k-1)); %Estimated variance of the variance of the estimate of the mutual information at full N
%             stdvar = sqrt(varS/N^2); %the error bars on our estimate of the variance
% --------------------------------------------------------------------------------------

            % 2019107 BC
            % Adding hack to filter mutual information results that are
            % within three S.D. from 0
            if ((MIs(1) - errThreshold*(variancePredicted^0.5)) > 0 || errThreshold == 0) && (MIs(1) > 0)
                r.mi = MIs(1);
                r.err = variancePredicted^.5;
            else
                r.mi = 0;
                r.err = 0;
            end
            
            
%             % return MI value and error estimation
%             r.mi = MIs(1);
%             r.err = stdvar^0.5;
        end
        
        function r = find_k_value(obj)
            % determine best k-value to use
            data = cell2mat(obj.mi_data);
            
            k_vals = [obj.mi_data{:,4}];
            ks = unique(k_vals);
            
            weighted_dataFrac = zeros(length(ks),4);
            
            % Check for stability across at least four data fractions
            for ik = ks
                dataFracs = find(data(:,4) == ik);
                for iFrac = 1:length(dataFracs)
                    
                end
            end
            
            
            
%             % find k-value that is least sensitive to changing k-value
%             k_mi = zeros(1,length(ks));
%             for i=1:length(ks)
%                 k_ix = find([obj.mi_data{:,4}] == ks(i) & [obj.mi_data{:,3}] == 1);
%                 k_mi(i) = obj.mi_data{k_ix,1};
%             end
%             if length(ks) > 1
%                 for i=1:length(ks)
%                     if i == 1
%                         weighted_dataFrac(i,1) = (k_mi(i+1)-k_mi(i))/k_mi(i);
%                     elseif i == length(ks)
%                         weighted_dataFrac(i,1) = (k_mi(i) - k_mi(i-1))/k_mi(i);
%                     else
%                         weighted_dataFrac(i,1) = mean((k_mi([i-1 i+1])-k_mi(i))/k_mi(i));
%                     end
%                 end
% 
%                 % find k-value with stable data fractions
% 
%                 % first column
%                 % second column: percent difference from data frac = 1
%                 % third column: std of percent difference from data frac = 1
%                 % fourth column: std of data frac std
% 
%                 for i=1:length(ks)
%                     k_ixs = find(k_vals == ks(i));
%                     mi_vals = [obj.mi_data{k_ixs,1}];
%                     perc_mi_vals = (mi_vals - mi_vals(1))/mi_vals(1); % percent difference from 1 data fraction
% 
%                     weighted_dataFrac(i,2) = sign(mean(perc_mi_vals(2:end)))*mean(perc_mi_vals(2:end).^2)^0.5;
%                     weighted_dataFrac(i,3) = std(perc_mi_vals(2:end));
% 
%                     weighted_dataFrac(i,4) = std([obj.mi_data{k_ixs,2}]);
%                 end
% 
%                 % provide some quantification of confidence?
% 
%                 [~,ix] = min(sum(weighted_dataFrac.^2,2));
%                 obj.opt_k = ks(ix);
%             else
%                 obj.opt_k = obj.k_values;
%             end
        end
        
        function r = fractionate_data(obj, k)
            % return cell array of fractionated datasets with x-data,
            % y-data, k-value, and ix
            n = length(obj.x);
            r = cell(sum(1:obj.data_fracs),4);
            for frac_n = 1:obj.data_fracs
                % determine length of subsample
                a = randperm(n);
                l = round(linspace(0,n,frac_n+1));
                
                % generate unique key to track each simulation
                while 1
                    key = num2str(dec2hex(round(rand(1)*100000)));
                    if ~any(strcmp(r(:,4), key))
                        break;
                    end
                end
                    
                % select subsample of data and assign data and params to data cell array
                for j=1:frac_n 
                    xT = obj.x(a(l(j)+1:l(j+1)));
                    yT = obj.y(a(l(j)+1:l(j+1)));
                    r(sum(1:(frac_n-1))+j,:) = {xT yT k key};
                end
            end
        end
        
    end
end
