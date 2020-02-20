classdef mi_analysis < handle
    %  MI_KSG_data_analysis is a parent class used to set up a separate object for each pair of variables to feed into the MI calculations
    % 
    properties
        
        varNames % indicates which data(s) to pull
        
        objData % Reference to which data object to pull from
        
        objBehav % Reference to which behavior object to pull from (optional)
        
        % BC: This will be a list/cell array of objMIcore instances (may need to index)
        % BC: cell array with structure: {{objMICore} {coeff} {k-value} {coreID}}
        arrMIcore % Reference to MIcore object 
	
        sim_manager % Sim manager reference object
        
        append % Specify whether to re-run analysis or just for k-values that have not been previously included
        
        verbose % level of output for progress and troubleshooting/debugging
        
        notes %Indicates how much data has been omitted (optional)
        
    end

    methods
        function obj = mi_analysis(objData, objBehav, varNames, varargin)
            % This funtion inputs the data object reference and variable references
            
            % Instantiate input parser
            p = inputParser;
            
            % Set up required inputs
            validate_objData = @(x) assert(isa(x, 'mi_data'), 'objData must be a valid data object');
            p.addRequired('objData', validate_objData);
            
            validate_objBehav = @(x) assert(isa(x, 'mi_data'), 'objBehav must be a valid behavior object');
            p.addRequired('objBehav', validate_objBehav);
            
            validate_varNames = @(x) assert(iscell(x), 'varNames must be a cell array of strings');
            p.addRequired('varNames', validate_varNames);
            
            
            % Set up optional input
            
            % append
            default_append = true;
            validate_append = @(x) assert(islogical(x), 'append must be a logical value');
            p.addParameter('append', default_append, validate_append);
            
            % verbose
            default_verbose = 1;
            validate_verbose = @(x) assert(isnumeric(x) && rem(x,1) == 0, 'verbose must be an integer');
            p.addParameter('verbose', default_verbose, validate_verbose);
            
            % Parse the inputs
            % Set up InputParser to handle extra inputs from subclasses
            p.KeepUnmatched = 1;
            p.parse(objData, objBehav, varNames, varargin{:});
            
            obj.objData = p.Results.objData;
            obj.varNames = p.Results.varNames;
            obj.objBehav = p.Results.objBehav;
            obj.append = p.Results.append;
            obj.verbose = p.Results.verbose;
            
            % Temporarily set arrMIcore and instantiate a sim_manager
            % object
            obj.arrMIcore = {};
            obj.sim_manager = mi_ksg_sims(1,3);
            
            if obj.verbose > 0; disp([newline 'mi_analysis instantiated']); end
        end

        function buildMIs(obj, mi_data)

            % Set up empty array for obj.arrMIcore
            obj.arrMIcore = cell(size(mi_data,1),4);
            
            % Define groups to fill in arrMIcore
            xGroups = mi_data{1};
            yGroups = mi_data{2};
            coeffs = mi_data{3};
            
            v = obj.verbose;
            
            if v > 1; disp([newline 'Assigning ID to each subgroup...']); end
            
            for iGroup = 1:size(xGroups,1)
                x = xGroups{iGroup,1};
                y = yGroups{iGroup,1};

              
                while 1 % generate random key to keep track of which MI calculations belong together
                    key = num2str(dec2hex(round(rand(1)*100000)));
                    % break the while loop if the key has not already been
                    % assigned.
                    if iGroup == 1
                        break
                    elseif ~ismember({obj.arrMIcore{1:end,4}}, key)
                        break
                    end
                end
                
                
                % RC: Why do we set the k values in the core object and in
                % the arrMIcore?
                core1 = mi_ksg_core(obj.sim_manager, x, y, 'ks_arr', 1:9, 'opt_k',1, 'append', obj.append, 'verbose', obj.verbose);
                if v > 2; disp([newline '--> Core object instantiated for group: ' num2str(iGroup)]); end
                
                if v > 2; disp([newline '--> Group ' num2str(iGroup) ' has ' num2str(max(size(x))) ' data points']); end
                
	            obj.arrMIcore(iGroup,:) = {core1 coeffs{iGroup,1} 0 key};
                
                if v > 2; disp([newline '--> arrMIcore assigned']); end
	            % BC: The obj.findMIs function basically calls run_sims
            end
            
            mi_ksg_viz.audit_plots(obj);
            
	    % Sets up an MIcore object to calculate the MI values, and pushes the
	    % data from this object to the MIcore process. 
        
            if v > 0; disp(['COMPLETE: Added all data to arrMIcore']); end
        end
        
        function calcMIs(obj)
            v = obj.verbose;
            if v > 0; disp('Calculating mutual information...'); end
            run_sims(obj.sim_manager);
        end
    end
end
