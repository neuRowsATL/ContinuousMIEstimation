classdef calc_timing_count < mi_analysis
    %Each of these objects sets the stage to calculate the mutual
    %information between spike count and behavior and stores the results of
    %the calculation. 
    %   Detailed explanation goes here
    
    properties
        n_timeBase
    end
    
    methods
        function obj = calc_timing_count(objData, objBehav, varNames, varargin)
            % Required arguments: objData, objBehav,  varNames
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
            default_n_timeBase = 'time';
            valid_n_timeBases = {'time', 'phase'};
            validate_n_timeBase = @(x) assert(ischar(x) && ismember(x, valid_n_timeBases), 'n_timeBase must be: time, phase');
            p.addParameter('n_timeBase', default_n_timeBase, validate_n_timeBase); 
            
            % Prepare InputParser to parse only desired inputs
            p.KeepUnmatched = 1;
            p.parse(objData, objBehav, varNames, varargin{:});
                                   
            % Define validated inputs to parent constructor
            objData = p.Results.objData;
            objBehav = p.Results.objBehav;
            varNames = p.Results.varNames;
            
            % One more validation: Check that varNames references valid fields of objData
            for ivarNames = 1:length(varNames)
                assert(isfield(objData.data , varNames{ivarNames}), ['varName: ' varNames{ivarNames} ' is not a valid field of the neural data object']); 
            end
            
            % Call parent constructor
            obj@mi_analysis(objData, objBehav, varNames, varargin{:});
            
            % Define timebase property of subclass object
            obj.n_timeBase = p.Results.n_timeBase;

        end
        
        function buildMIs(obj)
            % Build the data and core objects necessary to run the sim manager for this analysis class            
            % First, segment neural data into breath cycles
            v = obj.verbose;
            
            % Find spike timings in cycles for neuron 1 
            x_name  = obj.varNames{1};
            x = obj.objData.get_spikes('name', x_name , 'format', 'timing', 'cycleTimes', obj.objBehav.data.cycleTimes.data, 'timeBase', obj.n_timeBase);
            
%             % Audit Check
%             if sum(sum(~isnan(x))) ~= (sum(~isnan(obj.objData.data.(obj.varNames{1}).data)) - (sum(obj.objData.data.(obj.varNames{1}).data < obj.objBehav.data.cycleTimes.data(1,1) | obj.objData.data.(obj.varNames{1}).data > obj.objBehav.data.cycleTimes.data(end,2))))
%                 error('Error: N Spikes in x do not match that expected from objData.varNames{1}.');
%             end
           
            % Find different subgroups
            xCounts = obj.objData.get_spikes('name', x_name , 'format', 'count', 'cycleTimes', obj.objBehav.data.cycleTimes.data );
            xConds= unique(xCounts);
            
%             % Audit Check
%             if sum(xCounts) ~= (sum(~isnan(obj.objData.data.(obj.varNames{1}).data)) - (sum(obj.objData.data.(obj.varNames{1}).data < obj.objBehav.data.cycleTimes.data(1,1) | obj.objData.data.(obj.varNames{1}).data > obj.objBehav.data.cycleTimes.data(end,2))))
%                 error('Error: Spike Counts for x do not match that expected from objData.varNames{1}.');
%             end
%             if sum(xCounts) ~=  sum(sum(~isnan(x)))
%                 error('Error: Spike counts for x do not matach N Spikes in x.'); 
%             end


            % Next segment other neuron into cycles and find the count
            y_name = obj.varNames{2};
            y = obj.objData.get_spikes('name', y_name , 'format', 'count', 'cycleTimes', obj.objBehav.data.cycleTimes.data );

%             % Audit Check: number of spikes detected
%             if sum(y) ~= (sum(~isnan(obj.objData.data.(obj.varNames{2}).data)) - (sum(obj.objData.data.(obj.varNames{2}).data < obj.objBehav.data.cycleTimes.data(1,1) | obj.objData.data.(obj.varNames{2}).data > obj.objBehav.data.cycleTimes.data(end,2))))
%                 error('Error: N Spikes in x do not match that expected from objData.varNames{1}.');
%             end

            % AS WRITTEN- we put each subgroup for the calculation into an array. 
            % NOTE currently as this code is written, we dont worry about data limitations. 
            xGroups = {};
            coeffs = {};
            yGroups = {};

            % Segment x and y data into roups based on x spike count
            % Set Group counter
            groupCount = 1;
            noteCount = 1;
            omitCoeff = [];
            
            for iCond = 1:length(xConds)
                Cond = xConds(iCond);
                groupIdx = find(xCounts == Cond);
                if Cond == 0
                    % Find ratio and percent of data that will be omitted. 
                    num = length(groupIdx);
                    groupRatio = num/length(xCounts);
                    percent = groupRatio*100;
                    
                    % Document how much data is omitted.
                    note = strcat('Omitting ', num2str(percent), 'percent of cycles because zero spikes');
                    disp(note)
                    obj.notes{noteCount,1} = note;
                    
                    % Keep track of total omitted ratio
                    omitCoeff(noteCount) = groupRatio;
                    
                    % Increase note counter
                    noteCount = noteCount + 1;
                    % When there are zero spikes, the MI from timing is zero. This
                    % can't be accounted for in the calculation because
                    % there are no time values to send to MIxnyn.
                    % Therefore, we are setting the coeff for this group to
                    % zero. The percent will be accounted for in the rest
                    % of the Coeffs (the Coeffs will sum to 1 - n(zero)
                    continue
                elseif Cond > sum(xCounts == Cond)
                    % Find ratio and percent of data that will be omitted. 
                    num = length(groupIdx);
                    groupRatio = num/length(xCounts);
                    percent = groupRatio*100;
                    
                    % Document how much data is omitted. 
                    note = strcat('Omitting ', num2str(percent), 'percent of cycles, where Cond = ' , num2str(Cond), 'because more spikes than data.');
                    disp(note)
                    obj.notes{noteCount,1} = note;
                    
                    % Keep track of total omitted ratio
                    omitCoeff(noteCount) = groupRatio;
                    
                    % Increase note counter
                    noteCount = noteCount + 1;                    
                    continue
                end 
                % Define xdata for iCond
                ixGroup =  x(groupIdx,1:Cond);
            
                % Find coeff corresponding to iCond
                
                xGroups{groupCount,1} = ixGroup;
                

                
                coeffs{groupCount,1} = size(ixGroup,1)/length(xCounts);
                
                % Define y data for iCond
                yGroups{groupCount,1} = y(groupIdx)';
                
                % Increase group counter
                groupCount = groupCount + 1;
            end
            % Audit: Check that omit coeffs and group coeffs sum to 1 with a very small tolerance to account for matlab rounding error. 
            if ~ismembertol((sum(cell2mat(coeffs)) + sum(omitCoeff)), 1, 1e-12); error('Error: Sum of coeffs and omitted data ratios does not equal 1'); end
            
            % Audit: Is there still data left to analyze?
            if ismembertol(sum(omitCoeff),1, 1e-12); error('Error: All subgroups were omitted. Not enough data');end
            
            % Call parent class buildMIs()
            buildMIs@mi_analysis(obj, {xGroups yGroups coeffs});          
            
        end
    end
end

