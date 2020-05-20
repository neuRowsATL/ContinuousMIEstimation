<<<<<<< HEAD
classdef mi_analysis < handle
    %  MI_KSG_data_analysis is a parent class used to set up a separate object for each pair of variables to feed into the MI calculations
    % 
    properties
        
        varNames % indicates which data(s) to pull
        
        objData % Reference to which data object to pull from
        
        objBehav % Reference to which behavior object to pull from (optional)
        
        % BC: This will be a list/cell array of objMIcore instances (may need to index)
        % cell array with structure:{{objMICore} {coeff} {k-value} {MIestimate} {Error} {coreID}}
        arrMIcore % Reference to MIcore object 
	
        sim_manager % Sim manager reference object
        
        append % Specify whether to re-run analysis or just for k-values that have not been previously included
        
        verbose % level of output for progress and troubleshooting/debugging
        
        notes %Indicates how much data has been omitted (optional)
        
        reparam %Do you reparameterize
        
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
            
            default_reparam = 2;
            validate_reparam = @(x) assert(isnumeric(x) && rem(x,1) == 0, 'reparam must be an integer');
            p.addParameter('reparam', default_reparam, validate_reparam);
            
            % Parse the inputs
            % Set up InputParser to handle extra inputs from subclasses
            p.KeepUnmatched = 1;
            p.parse(objData, objBehav, varNames, varargin{:});
            
            obj.objData = p.Results.objData;
            obj.varNames = p.Results.varNames;
            obj.objBehav = p.Results.objBehav;
            obj.append = p.Results.append;
            obj.verbose = p.Results.verbose;
            obj.reparam = p.Results.reparam;
            
            % Temporarily set arrMIcore and instantiate a sim_manager
            % object
            obj.arrMIcore = {};
            obj.sim_manager = mi_ksg_sims(1,3);
            
            if obj.verbose > 0; disp([newline 'mi_analysis instantiated']); end
        end

        function buildMIs(obj, mi_data)

            % Set up empty array for obj.arrMIcore
            obj.arrMIcore = cell(size(mi_data,1),6);
            
            % Define groups to fill in arrMIcore
            xGroups = mi_data{1};
            yGroups = mi_data{2};
            coeffs = mi_data{3};
            
            v = obj.verbose;
            
            if v > 1; disp([newline 'Assigning ID to each subgroup...']); end
            
            for iGroup = 1:size(xGroups,1)
                x = xGroups{iGroup,1};
                y = yGroups{iGroup,1};

                if obj.reparam >= 1
                    if any(rem(x,1) ~= 0) % If continuous
                        for rep = 1:size(x,2) % Reparameterize x data
                            x(:,rep) = reparameterize_data(x(:,rep));
                        end
                    end
                    
                    if obj.reparam >= 2
                        if any(rem(y,1) ~= 0)
                            for rep = 1:size(y,2) % Reparameterize y data
                                y(:,rep) = reparameterize_data(y(:,rep));
                            end   
                        end
                    end
                end

                while 1 % generate random key to keep track of which MI calculations belong together
                    key = num2str(dec2hex(round(rand(1)*100000)));
                    % break the while loop if the key has not already been
                    % assigned.
                    if iGroup == 1
                        break
                    elseif ~ismember({obj.arrMIcore{1:end,6}}, key)
                        break
                    end
                end
                
                
                % RC: Why do we set the k values in the core object and in
                % the arrMIcore?
                core1 = mi_ksg_core(obj.sim_manager, x, y, 'ks_arr', 1:9, 'opt_k',1, 'append', obj.append, 'verbose', obj.verbose);
                if v > 2; disp([newline '--> Core object instantiated for group: ' num2str(iGroup)]); end
                
                if v > 2; disp([newline '--> Group ' num2str(iGroup) ' has ' num2str(max(size(x))) ' data points']); end
                
	            obj.arrMIcore(iGroup,:) = {core1 coeffs{iGroup,1} NaN NaN NaN key};
                
                if v > 2; disp([newline '--> arrMIcore assigned']); end
	            % BC: The obj.calcMIs function basically calls run_sims
            end
            
            % Audit plots
            % NOTE- Move these to functions in mi_ksg_viz eventually. 
            for iGroup = 1:size(xGroups)
%                 clear r_2 out_p
                coreObj = obj.arrMIcore{iGroup,1};
                
                if v > 4
                    % FOR NOW, NO AUDIT PLOTS FOR BEHAVIOR SUBCLASSES
                    if contains(class(obj), 'behav')
                        continue
                    else
                        
                        % Check for data type
                        % Histograms do not depend on data type
                        % First make histogram with x data
                        % RC 20191213: We should come back and set specific bin widths here.
                        % The only issue we may run into is 
                        x = coreObj.x;
%                         figure()
%                         histogram(x)
%                         hold on
%                         xlabel('X Value (binned)')
%                         ylabel('N Cycles')
%                         title('Histogram for X')

                        % Histogram for y
                        y = coreObj.y;
%                         figure()
%                         histogram(y)
%                         hold on
%                         xlabel('Y Value (binned)')
%                         ylabel('N Cycles')
%                         title('Histogram for Y')

                        % Also skip audit plots for data where both x and y are multi-dimensional
                        if all(size(x) > 1) & all(size(y) > 1)
                            continue
                        else
                            % Check for discrete data in both variables
                            if all(rem(x,1) == 0) & all(rem(y,1) == 0)
                                % For discrete data, plot a jittered histogram
                                
                                % Add noise for joint histogram
                                x_plot = x + 0.2*rand(size(x));
                                y_plot = y + 0.2*rand(size(y));
                                
                                % Make figure
%                                 figure()
%                                 hold on
%                                 plot(x_plot , y_plot, 'x')
%                                 hold on
%                                 xlabel('Discrete Value: X')
%                                 ylabel('Discrete Value: Y')
%                                 title('P(X,Y) Mixed Joint Distribution')                                
                            elseif all(rem(x,1) == 0) | all(rem(y,1) == 0)
                                % Add jitter only to the variable that is discrete, which for our data, will always be the second variable.
                                if all(rem(x,1) == 0)
                                    x_plot = x + 0.2*rand(size(x));
                                    x_L = 'Discrete Value: X';
                                else
                                    x_plot = x;
                                    x_L = 'Continuous Value: X';
                                end
                                if all(rem(y,1) == 0)
                                    y_plot = y + 0.2*rand(size(y));
                                    y_L = 'Discrete Value: Y';
                                else
                                    y_plot = y;
                                    y_L = 'Continuous Value: Y';
                                end
                                                                
                                % Make figure
                                figure()
                                hold on
                                plot(x_plot, y_plot, 'x')
                                hold off
                                xlabel(x_L)
                                ylabel(y_L)
                                title('P(X,Y) Mixed Joint Distribution')
                            else
                                % The assumption is that both distributions are continuous if neither of the above if statements are true.
                                
                                % Make figure
                                figure()
                                plot(x,y, 'x')
                                hold on
                                xlabel('Continuous Variable: X')
                                ylabel('Continuous Variable: Y')
                                title('P(X,Y) Continuous Joint Distribution')
                            end                            
                        end
                    end
                end
            end
                
	    % Sets up an MIcore object to calculate the MI values, and pushes the
	    % data from this object to the MIcore process. 
        
            if v > 0; disp(['COMPLETE: Added all data to arrMIcore']); end
        end
        
        function calcMIs(obj)
            v = obj.verbose;
            if v > 0; disp('Calculating mutual information...'); end
            run_sims(obj.sim_manager);
        end

        function getMIs(obj)
            for iCores = 1:size(obj.arrMIcore,1)
                core = obj.arrMIcore{iCores,1};

                % Set the warning message to empty
                lastwarn('')

                % Call core function to find k value and output MI and error
                r = core.get_mi(-1);

                % Grab last warning message
                w = lastwarn;

                if isempty(w)
                    obj.arrMIcore{iCores,3} = r.k;
                else
                    obj.arrMIcore{iCores,3} = w;
                end

                % Set temporary MI and error values
                obj.arrMIcore{iCores,4} = r.mi;
                obj.arrMIcore{iCores,5} = r.err;
                
            end           
            
        end
        
        function r = returnMIs(obj)
            % NOTE WE NEED TO EDIT THIS AS WE MAKE DESIGN DECISIONS!!
            r.mi = nansum(cell2mat(obj.arrMIcore(:,4)).*cell2mat(obj.arrMIcore(:,2)));
            r.err = nansum(cell2mat(obj.arrMIcore(:,5)).*cell2mat(obj.arrMIcore(:,2)));
        end
        
    end
end
=======
classdef mi_analysis < handle
    %  MI_KSG_data_analysis is a parent class used to set up a separate object for each pair of variables to feed into the MI calculations
    % 
    properties
        
        varNames % indicates which data(s) to pull
        
        objData % Reference to which data object to pull from
        
        objBehav % Reference to which behavior object to pull from (optional)
        
        % BC: This will be a list/cell array of objMIcore instances (may need to index)
        % cell array with structure:{{objMICore} {coeff} {k-value} {MIestimate} {Error} {coreID}}
        arrMIcore % Reference to MIcore object 
	
        sim_manager % Sim manager reference object
        
        append % Specify whether to re-run analysis or just for k-values that have not been previously included
        
        verbose % level of output for progress and troubleshooting/debugging
        
        notes %Indicates how much data has been omitted (optional)
        
        reparam %Do you reparameterize
        
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
            
            default_reparam = 2;
            validate_reparam = @(x) assert(isnumeric(x) && rem(x,1) == 0, 'reparam must be an integer');
            p.addParameter('reparam', default_reparam, validate_reparam);
            
            % Parse the inputs
            % Set up InputParser to handle extra inputs from subclasses
            p.KeepUnmatched = 1;
            p.parse(objData, objBehav, varNames, varargin{:});
            
            obj.objData = p.Results.objData;
            obj.varNames = p.Results.varNames;
            obj.objBehav = p.Results.objBehav;
            obj.append = p.Results.append;
            obj.verbose = p.Results.verbose;
            obj.reparam = p.Results.reparam;
            
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
            
            if v > 4; mi_ksg_viz.audit_plots(obj); end
            
	    % Sets up an MIcore object to calculate the MI values, and pushes the
	    % data from this object to the MIcore process. 
        
            if v > 0; disp(['COMPLETE: Added all data to arrMIcore']); end
        end
        
        function calcMIs(obj)
            v = obj.verbose;
            if v > 0; disp('Calculating mutual information...'); end
            run_sims(obj.sim_manager);
        end

        function getMIs(obj)
            for iCores = 1:size(obj.arrMIcore,1)
                core = obj.arrMIcore{iCores,1};

                % Set the warning message to empty
                lastwarn('')

                % Call core function to find k value and output MI and error
                r = core.get_mi(-1);

                % Grab last warning message
                w = lastwarn;

                if isempty(w)
                    obj.arrMIcore{iCores,3} = r.k;
                else
                    obj.arrMIcore{iCores,3} = w;
                end

                % Set temporary MI and error values
                obj.arrMIcore{iCores,4} = r.mi;
                obj.arrMIcore{iCores,5} = r.err;
                
            end           
            
        end
        
        function r = returnMIs(obj)
            % NOTE WE NEED TO EDIT THIS AS WE MAKE DESIGN DECISIONS!!
            r.mi = nansum(cell2mat(obj.arrMIcore(:,4)).*cell2mat(obj.arrMIcore(:,2)));
            r.err = nansum(cell2mat(obj.arrMIcore(:,5)).*cell2mat(obj.arrMIcore(:,2)));
        end
        
    end
end
>>>>>>> 012cb64efdd5ed3fe7d4fefe2e47a1e7139205f7
