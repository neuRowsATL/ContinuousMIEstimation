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
        
        verbose % level of output for progress and troubleshooting/debugging
        
        notes %Indicates how much data has been omitted (optional)
        
    end

    methods
        function obj = mi_analysis(objData, varNames, varargin)
            % This funtion inputs the data object reference and variable references
            
            % Instantiate input parser
            p = inputParser;
            
            % Set up required inputs
            validate_objData = @(x) assert(isobject(x), 'objData must be a valid data object');
            p.addRequired('objData', validate_objData);
            
            validate_varNames = @(x) assert(iscell(x), 'varNames must be a cell array of strings');
            p.addRequired('varNames', validate_varNames);
            
            % Set up optional inputs
            
            % objBehav
            % Default objBehav is empty struct for now? 
            default_objBehav = struct(); 
            validate_objBehav = @(x) assert(isobject(x) | isstruct(x), 'objBehav must be a valid behavior object');
            p.addOptional('objBehav', default_objBehav, validate_objBehav);
            
            % verbose
            default_verbose = 1;
            validate_verbose = @(x) assert(isnumeric(x) && rem(x,1) == 0, 'verbose must be an integer');
            p.addParameter('verbose', default_verbose, validate_verbose);
            
            % Parse the inputs
            p.parse(objData, varNames, varargin{:});
            
            obj.objData = p.Results.objData;
            obj.varNames = p.Results.varNames;
            obj.objBehav = p.Results.objBehav;
            obj.verbose = p.Results.verbose;
            
            % Temporarily set arrMIcore and instantiate a sim_manager
            % object
            obj.arrMIcore = {};
            obj.sim_manager = mi_ksg_sims(1,3);
            
            if obj.verbose > 0; disp([newline 'mi_analysis instantiated']); end
        end

        function buildMIs(obj, mi_data)
            % NOTE- we still need to change the default k-value for this. 
	        % BC: Move his for loop into the constructor for MI_KSG_data_analysis subclasses- DONE
            obj.arrMIcore = cell(size(mi_data,1),4);
            
            xGroups = mi_data{1};
            yGroups = mi_data{2};
            coeffs = mi_data{3};
            
            for iGroup = 1:size(xGroups,1)
                x = xGroups{iGroup,1};
                y = yGroups{iGroup,1};
		
		if obj.objData.reparamData == 1
                    for iDimension = 1:size(x,1)
                        x(iDimension,:) = obj.objData.reparameterizeData(x(iDimension,:));
                    end
                    for iDimension = 1:size(y,1)
                        y(iDimension,:) = obj.objData.reparameterizeData(y(iDimension,:));
                    end
                end
		
                % BC: Need to append new mi_core instance to the arrMICore object with associated information- DONE
                % RC-  Is it a problem that we name the core object the same thing each iteration? 
              
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
                core1 = mi_ksg_core(obj.sim_manager, x, y, 3:8, 0);
	            obj.arrMIcore(iGroup,:) = {core1 coeffs{iGroup,1} 0 key};
	            % BC: The obj.findMIs function basically calls run_sims
            end
	    % Sets up an MIcore object to calculate the MI values, and pushes the
	    % data from this object to the MIcore process. 
        end
        
        function calcMIs(obj)
            disp('Calculating mutual information...');
            run_sims(obj.sim_manager);
        end
    end
end
