classdef calc_timing_timing_behav < mi_analysis
    
    %Each of these objects sets the stage to calculate the mutual
    %information between spike timing of neuron 1 and spike timing of neuron 2 and stores the results of
    %the calculation. 
    
    % For this type of calculation, it is likely that we will have to be creative
    % maybe we will use the average ISI in a breath cycle
    % maybe we will use only the first spike
    % We will have to see how much data we have. 
    
    properties
        n_timeBase
        b_timeBase
        feature
        start
        dur
        nSamp
        nPC
    end
    
    methods
       function obj = calc_timing_timing_behav(objData,objBehav, varNames, varargin)
        % Set up input parser
            p = inputParser;
            
            % Set required inputs
            validate_objData = @(x) assert(isa(x, 'mi_data_neural'), 'objData must be a neural data subclass');
            p.addRequired('objData', validate_objData);
            
            validate_objBehav = @(x) assert(isa(x, 'mi_data_behavior'), 'objBehav must be a behavioral data subclass');
            p.addRequired('objBehav', validate_objBehav);
            
            validate_varNames = @(x) assert(iscell(x) && (length(x) == 2), 'varNames must be a cell of length 1');
            p.addRequired('varNames', validate_varNames);
                        
            % Set parameters
            default_n_timeBase = 'time';
            valid_n_timeBases = {'time', 'phase'};
            validate_n_timeBase = @(x) assert(ischar(x) && ismember(x, valid_n_timeBases), 'n_timeBase must be: time, phase');
            p.addParameter('n_timeBase', default_n_timeBase, validate_n_timeBase); 

                                    
            % Set parameters
            default_b_timeBase = 'phase';
            valid_b_timeBases = {'time', 'phase'};
            validate_b_timeBase = @(x) assert(ischar(x) && ismember(x, valid_timeBases), 'b_timeBase must be: time, phase');
            p.addParameter('b_timeBase', default_b_timeBase, validate_b_timeBase); 

            default_feature = 'residual';
            valid_feature = {'raw', 'pca', 'residual'};
            validate_feature = @(x) assert(ischar(x) && ismember(x, valid_feature), 'feature  must be: raw, pca, residual');
            p.addParameter('feature', default_feature, validate_feature); 

            default_start = pi/2;
            validate_start = @(x) assert((0 < x < 2*pi) || isinteger(x), 'start must be in units of radians or ms, and must match b_timeBase');
            p.addParameter('start', default_start, validate_start); 

            default_dur = pi;
            validate_dur = @(x) assert((0 < x < 2*pi) || isinteger(x), 'dur must be in units of radians or ms, and must match b_timeBase');
            p.addParameter('dur', default_dur, validate_dur); 

            default_nSamp = 11;
            validate_nSamp = @(x) assert(isinteger(x), 'nSamp must be an integer');
            p.addParameter('nSamp', default_nSamp, validate_nSamp);

            default_nPC = 3;
            validate_nPC = @(x) assert(isinteger(x), 'nPC must be an integer');
            p.addParameter('nPC', default_nPC, validate_nPC);

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
            
            % Define other properties of subclass object
            obj.n_timeBase = p.Results.n_timeBase;
            obj.b_timeBase = p.Results.b_timeBase;
            obj.feature = p.Results.feature;
            obj.start = p.Results.start;
            obj.dur = p.Results.dur;
            obj.nSamp = p.Results.nSamp;
            obj.nPC = p.Results.nPC;
        end
        
        function buildMIs(obj)
            % So I propose that we use this method to prep the
            % count_behavior data for the MI core and go ahead and run MI
            % core from here. Then we can use the output of MI core to fill
            % in the MI, kvalue, and errors.
            

            % First, segment neural data into breath cycles
            neuron = obj.vars(1);
            n1 = obj.objData.getTiming(neuron);
           
            % Find different subgroups for neuron 1
            n1Counts = obj.objData.getCount(neuron);
            n1Conds = unique(n1Counts);

            % Segment neuron 2 into breath cycles
            % Segment neuron 2 into breath cycles
            neuron = obj.vars(2);
            n2 = obj.objData.getTiming(neuron);

            % Find different subgroups for neuron 2
            n2Counts = obj.objData.getCount(neuron);
            n2Conds = unique(n2Counts);
            
            % Segment behavioral data into cycles

            if nargin < 6
                y = obj.objData.processBehavior();
            elseif nargin == 6
                y = obj.objData.processBehavior();
            end
            
            %Both neurons collectively will make up the x group. We will
            %concatonate each condition. 
            xGroups = {};
            yGroups = {};
            coeffs = {};
            groupCounter = 1;
            noteCount = 1;
            for in1Cond = 1:length(n1Conds)
                n1Cond = n1Conds(in1Cond);
                n1groupIdx = find(n1Counts == n1Cond);
                if n1Cond == 0
                    num = sum(n1Counts == n1Cond);
                    ratio = (num/length(n1Counts))*100;
                    note = strcat('Omitting ', num2str(ratio), ' percent of cycles because zero spikes in x.');
                    disp(note)
                    obj.notes{noteCount,1} = note;
                    noteCount = noteCount + 1;
                    continue
                end
                for in2Cond = 1:length(n2Conds)
                    n2Cond = n2Conds(in2Cond);
                    n2groupIdx = find(n2Counts == n2Cond);
                    xgroupIdx = intersect(n1groupIdx,n2groupIdx);
                    if n2Cond == 0
                        num = sum(n2Counts == n2Cond);
                        ratio = (num/length(n2Counts))*100;
                        note = strcat('Omitting ', num2str(ratio), ' percent of cycles because zero spikes in y.');
                        disp(note)
                        obj.notes{noteCount,1} = note;
                        noteCount = noteCount + 1;
                        continue
                    elseif n1Cond + n2Cond > length(xgroupIdx)
                        num = sum(n2Counts == n2Cond);
                        ratio = (num/length(n2Counts))*100;
                        note = strcat('Omitting ', num2str(ratio), ' percent of cycles,','where n1Cond = ',num2str(n1Cond), ' and n2Cond = ', num2str(n2Cond), 'because more spikes than data.');
                        disp(note)
                        obj.notes{noteCount,1} = note;
                        noteCount = noteCount + 1;
                        continue
                    end
%                     elseif n1Cond > length(xgroupIdx)
%                         num = sum(n1Counts == n1Cond);
%                         ratio = (num/length(n1Counts))*100;
%                         note = strcat('Omitting ', num2str(ratio), ' percent of cycles, where n1Cond = ', num2str(n1Cond), 'because more spikes than data.');
%                         disp(note)
%                         obj.notes{noteCount,1} = note;
%                         noteCount = noteCount + 1;
%                         continue 
%                     elseif n2Cond > length(xgroupIdx)
%                         num = sum(n2Counts == n2Cond);
%                         ratio = (num/length(n2Counts))*100;
%                         note = strcat('Omitting ', num2str(ratio), ' percent of cycles, where n2Cond = ', num2str(n2Cond), 'because more spikes than data.');
%                         disp(note)
%                         obj.notes{noteCount,1} = note;
%                         noteCount = noteCount + 1;
%                         continue   
                    n1Group = n1(xgroupIdx,1:n1Cond);
                    n2Group = n2(xgroupIdx,1:n2Cond);
                    xGroup = [n1Group,n2Group];

                    xGroups{groupCounter,1} = xGroup;
                    yGroup = y(xgroupIdx,1:end);
                    if length(xGroup) ~= length(yGroup)
                        keyboard
                    end
                    yGroups{groupCounter,1} = yGroup;
                    coeffs{groupCounter,1} = length(xgroupIdx)/length(n1Counts);
                    groupCounter = groupCounter + 1;
                end
                
            end
            buildMIs@mi_analysis(obj, {xGroups yGroups coeffs}); 
            
        end
    end
end


