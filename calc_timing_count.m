classdef calc_timing_count < mi_analysis
    %Each of these objects sets the stage to calculate the mutual
    %information between spike count and behavior and stores the results of
    %the calculation. 
    %   Detailed explanation goes here
    
    properties
        timebase
    end
    
    methods
        function obj = calc_timing_count(objData, objBehav, varNames, varargin)
            % Required arguments: objData, varNames
            % Check required inputs for validity using input parser
 
            % Set up input parser
            p = inputParser;
            
            % Set required inputs
            validate_objData = @(x) assert(isa(x, 'mi_data_neural'), 'objData must be a neural data subclass');
            p.addRequired('objData', validate_objData);
            
            validate_objBehav = @(x) assert(isa(x, 'mi_data_behavior'), 'objBehav must be a behavioral data subclass');
            p.addRequired('objBehav', validate_objBehav);
            
            validate_varNames = @(x) assert(iscell(x) && (length(x) == 2), 'varNames must be a cell of length 2');
            p.addRequired('varNames', validate_varNames);
            
            
            % Set parameters
            default_timebase = 'time';
            valid_timebases = {'time', 'phase'};
            validate_timebase = @(x) assert(ischar(x) && ismember(x, valid_timebases), 'timebase must be: time, phase');
            p.addParameter('timebase', default_timebase, validate_timebase); 
            
            % Prepare InputParser to parse only desired inputs
            p.KeepUnmatched = 1;
            p.parse(objData, objBehav, varNames, varargin{:});
                                   
            % Define validated inputs to parent constructor
            objData = p.Results.objData;
            objBehav = p.Results.objBehav;
            varNames = p.Results.varNames;
            
            % One more validation: Check that varNames references valid fields of objData
            for ivarNames = 1:length(varNames)
                assert(isfield(objData.data , varNames{ivarNames}), ['varName: ' varNames{ivarNames} 'is not a valid field of the neural data object']); 
            end
            
            % Call parent constructor
            obj@mi_analysis(objData, objBehav, varNames, varargin{:});
            
            % Define timebase property of subclass object
            obj.timebase = p.Results.timebase;

        end
        
        function buildMIs(obj)
            % Build the data and core objects necessary to run the sim manager for this analysis class            
            % First, segment neural data into breath cycles
            v = obj.verbose;
            
            % Find total spike count in a cycle for neuron 1 
            x_name  = obj.varNames{1};
            x = obj.objData.get_spikes('name', x_name , 'format', 'timing', 'cycleTimes', obj.objBehav.data.cycleTimes.data, 'timebase', obj.timebase);
           
            % Find different subgroups
            xCounts = obj.objData.get_spikes('name', x_name , 'format', 'count', 'cycleTimes', obj.objBehav.data.cycleTimes.data );
            xConds= unique(xCounts);

            % Next segment other neuron into cycles and find the count
            y_name = obj.varNames{2};
            y = obj.objData.get_spikes('name', y_name , 'format', 'count', 'cycleTimes', obj.objBehav.data.cycleTimes.data );


            % AS WRITTEN- we put each subgroup for the calculation into an array. 
            % NOTE currently as this code is written, we dont worry about data limitations. 
            xGroups = {};
            coeffs = {};
            yGroups = {};

            % Segment x and y data into roups based on x spike count
            % Set Group counter
            groupCount = 1;
            noteCount = 1;
            for iCond = 1:length(xConds)
                Cond = xConds(iCond);
                groupIdx = find(xCounts == Cond);
                if Cond == 0
                    num = length(groupIdx);
                    ratio = (num/length(xCounts))*100;
                    note = strcat('Omitting ', num2str(ratio), 'percent of cycles because zero spikes');
                    disp(note)
                    obj.notes{noteCount,1} = note;
                    noteCount = noteCount + 1;
                    % When there are zero spikes, the MI from timing is zero. This
                    % can't be accounted for in the calculation because
                    % there are no time values to send to MIxnyn.
                    % Therefore, we are setting the coeff for this group to
                    % zero. The percent will be accounted for in the rest
                    % of the Coeffs (the Coeffs will sum to 1 - n(zero)
                    continue
                elseif Cond > sum(xCounts == Cond)
                    num = length(groupIdx);
                    ratio = (num/length(xCounts))*100;
                    note = strcat('Omitting ', num2str(ratio), 'percent of cycles, where Cond = ' , num2str(Cond), 'because more spikes than data.');
                    disp(note)
                    obj.notes{noteCount,1} = note;
                    noteCount = noteCount + 1;
                    continue
                end
                ixGroup =  x(groupIdx,1:Cond);
                xGroups{groupCount,1} = ixGroup;
                coeffs{groupCount,1} = size(ixGroup,1)/length(xCounts);
                yGroups{groupCount,1} = y(groupIdx)';
                groupCount = groupCount + 1;
            end
	    
            buildMIs@mi_analysis(obj, {xGroups yGroups coeffs});          
            
        end
    end
end

