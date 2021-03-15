classdef calc_timing_timing < mi_analysis
    %Each of these objects sets the stage to calculate the mutual
    %information between spike timing of neuron 1 and spike timing of neuron 2 and stores the results of
    %the calculation. 
    
    % For this type of calculation, it is likely that we will have to be creative
    % maybe we will use the average ISI in a breath cycle
    % maybe we will use only the first spike
    % We will have to see how much data we have. 
    
    properties
        n_timeBase
        discard_omittedData
            
    end
    
    methods
        function obj = calc_timing_timing(objData,objBehav, varNames, varargin)
        % Required Arguments: objData,objBehav,  varNames
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
            
            default_discard_omittedData = true;
            validate_discard_omittedData = @(x) assert(islogical(x), 'discard_omittedData must be a boolean value');
            p.addParameter('discard_omittedData', default_discard_omittedData, validate_discard_omittedData);
            
            
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
            obj.n_timeBase = p.Results.n_timeBase;
            obj.discard_omittedData = p.Results.discard_omittedData;

        end
        
        function buildMIs(obj)
            
            v = obj.verbose;
            
            % Find spike timings in cycles for neuron 1
            x_name = obj.varNames{1};

            cycles_interest = obj.objBehav.get_cycleTimes(obj.cycle_select);
            
            x = obj.objData.get_spikes('name', x_name , 'format', 'timing', 'cycleTimes', cycles_interest, 'timeBase', obj.n_timeBase);
% 
%             % Audit Check
%             if sum(sum(~isnan(x))) ~= (sum(~isnan(obj.objData.data.(obj.varNames{1}).data)) - (sum(obj.objData.data.(obj.varNames{1}).data < cycles_interest(1,1) | obj.objData.data.(obj.varNames{1}).data > cycles_interest(end,2))))
%                 error('Error: N Spikes in x do not match that expected from objData.varNames{1}.');
%             end

            % Find different subgroups for y
            xCounts = obj.objData.get_spikes('name', x_name , 'format', 'count', 'cycleTimes', cycles_interest );
            xConds= unique(xCounts);

%             % Audit Check
%             if sum(xCounts) ~= (sum(~isnan(obj.objData.data.(obj.varNames{1}).data)) - (sum(obj.objData.data.(obj.varNames{1}).data < cycles_interest(1,1) | obj.objData.data.(obj.varNames{1}).data > cycles_interest(end,2))))
%                 error('Error: Spike Counts for x do not match that expected from objData.varNames{1}.');
%             end
%             if sum(xCounts) ~=  sum(sum(~isnan(x)))
%                 error('Error: Spike counts for x do not matach N Spikes in x.'); 
%             end

            % Find spike timings in cycles for neuron 2
            y_name = obj.varNames{2};
            y = obj.objData.get_spikes('name', y_name , 'format', 'timing', 'cycleTimes', cycles_interest, 'timeBase', obj.n_timeBase);

%             % Audit Check
%             if sum(sum(~isnan(y))) ~= (sum(~isnan(obj.objData.data.(obj.varNames{2}).data)) - (sum(obj.objData.data.(obj.varNames{2}).data < cycles_interest(1,1) | obj.objData.data.(obj.varNames{2}).data > cycles_interest(end,2))))
%                 error('Error: N Spikes in y do not match that expected from objData.varNames{2}.');
%             end

            % Find different subgroups for y
            yCounts = obj.objData.get_spikes('name', y_name , 'format', 'count', 'cycleTimes', cycles_interest );
            yConds= unique(yCounts);

%             % Audit Check
%             if sum(yCounts) ~= (sum(~isnan(obj.objData.data.(obj.varNames{2}).data)) - (sum(obj.objData.data.(obj.varNames{2}).data < cycles_interest(1,1) | obj.objData.data.(obj.varNames{2}).data > cycles_interest(end,2))))
%                 error('Error: Spike Counts for y do not match that expected from objData.varNames{1}.');
%             end
%             if sum(yCounts) ~=  sum(sum(~isnan(y)))
%                 error('Error: Spike counts for y do not matach N Spikes in y.'); 
%             end

            % AS WRITTEN- we put each subgroup for the calculation into an array. 
            % NOTE currently as this code is written, we dont worry about data limitations.
            xGroups = {};
            yGroups = {};
            coeffs = {};
            
            % Set Group Counter
            noteCount = 1;
            groupCounter = 1;
%             omitCoeff = [];
            
            for ixCond = 1:length(xConds)
                xCond = xConds(ixCond);
                xgroupIdx = find(xCounts == xCond);
%                 if xCond == 0
%                     % Find ratio and percent of data that will be omitted.
%                     num = length(xgroupIdx);
%                     groupRatio = (num/length(xCounts));
%                     percent = groupRatio*100;
% 
%                     % Document how much data is omitted.
%                     note = strcat('Omitting ', num2str(percent), ' percent of cycles because zero spikes in x.');
%                     disp(note)
%                     obj.notes{noteCount,1} = note;
% 
%                     % Keep track of total omitted ratio
%                     omitCoeff(noteCount) =  groupRatio;
% 
%                     % Increase not counter.
%                     noteCount = noteCount + 1;
%                     continue
%                 end
                for iyCond = 1:length(yConds)
                    yCond = yConds(iyCond);
                    ygroupIdx = find(yCounts == yCond);
                    xygroupIdx = intersect(xgroupIdx,ygroupIdx);
%                     if yCond == 0
%                         num = length(xygroupIdx);
%                         groupRatio = (num/length(yCounts));
%                         percent = groupRatio * 100;
%                         note = strcat('Omitting ', num2str(percent),'where xCond = ', num2str(xCond),'and yCond = ',num2str(yCond), ' percent of cycles because zero spikes in y.');
%                         disp(note)
%                         obj.notes{noteCount,1} = note;
% 
%                         % Keep track of total omitted ratio
%                         omitCoeff(noteCount) = groupRatio;
% 
%                         % Increase group count
%                         noteCount = noteCount + 1;
%                         continue
%                         
%                     elseif xCond > length(xygroupIdx)
%                         num = length(xygroupIdx);
%                         groupRatio = (num/length(xCounts));
%                         percent = groupRatio * 100;
% 
%                         % Document how much data is omitted.
%                         note = strcat('Omitting ', num2str(percent), ' percent of cycles, where xCond = ', num2str(xCond),'and yCond = ', num2str(yCond), 'because more spikes in x than data.');
%                         disp(note)
%                         obj.notes{noteCount,1} = note;
% 
%                         % Keep track of total omitted ratio
%                         omitCoeff(noteCount) = groupRatio;
% 
%                         % Increase note counter
%                         noteCount = noteCount + 1;
%                         continue
%                         
%                     elseif yCond > length(xygroupIdx)
%                         num = length(xygroupIdx);
%                         groupRatio = (num/length(yCounts));
%                         percent = groupRatio * 100;
%                         
%                         % Document how much data is omitted. 
%                         note = strcat('Omitting ', num2str(percent), ' percent of cycles, where xCond = ',num2str(xCond), ' and yCond = ', num2str(yCond), 'because more spikes than data.');
%                         disp(note)
%                         obj.notes{noteCount,1} = note;
% 
%                         % Keep track of total omitted ratio
%                         omitCoeff(noteCount) = groupRatio;
% 
%                         % Increase note counter
%                         noteCount = noteCount + 1;
%                         continue   
%                     end
                    % Check for data in subgroup
                    coeff = size(xygroupIdx,2)/length(xCounts);
                    if coeff == 0
                        continue
                    else
                        xGroup = x(xygroupIdx,1:xCond);
                        xGroups{groupCounter,1} = xGroup;
                        yGroup = y(xygroupIdx,1:yCond);
                        yGroups{groupCounter,1} = yGroup;

                        coeffs{groupCounter,1} = coeff;
                        groupCounter = groupCounter + 1;
                    end
                end
            end

            % Audit: Check that omit coeffs and group coeffs sum to 1 with a very small tolerance to account for matlab rounding error. 
            if ~ismembertol(sum(cell2mat(coeffs)), 1, 1e-12); error('Error: Sum of coeffs and omitted data ratios does not equal 1'); end

%            % Audit: Is there still data left to analyze?
%            if ismembertol(sum(omitCoeff),1, 1e-12); error('Error: All subgroups were omitted. Not enough data');end
            
            buildMIs@mi_analysis(obj, {xGroups yGroups coeffs});     
            
        end
    end
end

