classdef calc_timing_count_behav < mi_analysis
    
    %Each of these objects sets the stage to calculate the mutual
    %information between spike timing of neuron 1 and spike timing of neuron 2 and stores the results of
%the calculation.
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
        function obj = calc_timing_count_behav(objData, objBehav, varNames, varargin)
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
            
        % Build the data and core objects necessary to run the sim manager for this analysis class
            
            % Find spike timings in cycles for neuron 1 
            n1_name  = obj.varNames{1};
            n1 = obj.objData.get_spikes('name', n1_name , 'format', 'timing', 'cycleTimes', obj.objBehav.data.cycleTimes.data, 'timeBase', obj.n_timeBase);
            
            % Audit Check
            if sum(sum(~isnan(n1))) ~= (sum(~isnan(obj.objData.data.(obj.varNames{1}).data)) - (sum(obj.objData.data.(obj.varNames{1}).data < obj.objBehav.data.cycleTimes.data(1,1) | obj.objData.data.(obj.varNames{1}).data > obj.objBehav.data.cycleTimes.data(end,2))))
                error('Error: N Spikes in n1 do not match that expected from objData.varNames{1}.');
            end
           
            % Find different subgroups
            n1Counts = obj.objData.get_spikes('name', n1_name , 'format', 'count', 'cycleTimes', obj.objBehav.data.cycleTimes.data );
            n1Conds= unique(n1Counts);
            
            % Audit Check
            if sum(n1Counts) ~= (sum(~isnan(obj.objData.data.(obj.varNames{1}).data)) - (sum(obj.objData.data.(obj.varNames{1}).data < obj.objBehav.data.cycleTimes.data(1,1) | obj.objData.data.(obj.varNames{1}).data > obj.objBehav.data.cycleTimes.data(end,2))))
                error('Error: Spike Counts for n1 do not match that expected from objData.varNames{1}.');
            end
            if sum(n1Counts) ~=  sum(sum(~isnan(n1)))
                error('Error: Spike counts for n1 do not matach N Spikes in n1.'); 
            end
            
            % First, segment neural data into breath cycles
            n2_name = obj.varNames{2}
            n2Counts =obj.objData.get_spikes('name', n2_name , 'format', 'count', 'cycleTimes', obj.objBehav.data.cycleTimes.data, 'timeBase', obj.n_timeBase);

            % Audit Check: n2
            if sum(n2Counts) ~= (sum(~isnan(obj.objData.data.(obj.varNames{2}).data)) - (sum(obj.objData.data.(obj.varNames{2}).data < obj.objBehav.data.cycleTimes.data(1,1) | obj.objData.data.(obj.varNames{2}).data > obj.objBehav.data.cycleTimes.data(end,2))))
                error('Error: Spike Counts for n2 do not match that expected from objData.varNames{1}.');
            end


            % Get the behavioral data for analysis
            y = get_behavior(obj.objBehav, obj.b_timeBase, obj.feature, obj.start, obj.dur, obj.nSamp,'nPC', obj.nPC );
                        
            %Both neurons collectively will make up the x group. We will
            %concatonate each condition.
            
            xGroups = {};
            yGroups = {};
            coeffs = {};
            
            iGroup = 1;
            noteCount = 1;
            omitCoeff = [];
            
            for in1Cond = 1:length(n1Conds)
                n1Cond = n1Conds(in1Cond);
                n1groupIdx = find(n1Counts == n1Cond);
                if n1Cond == 0
                    num = length(n1groupIdx);
                    groupRatio = (num/length(n1Counts));
                    percent = groupRatio * 100;
                    
                    note = strcat('Omitting ', num2str(percent), 'percent of cycles because zero spikes');
                    disp(note)
                    obj.notes{noteCount,1} = note;

                    % Keep track of total omitted ratio
                    omitCoeff(noteCount) = groupRatio;

                    % Increase note counter. 
                    noteCount = noteCount + 1;
                    
                    % When there are zero spikes, the MI from timing is zero. This
                    % can't be accounted for in the calculation because
                    % there are no time values to send to MIxnyn.
                    % Therefore, we are setting the coeff for this group to
                    % zero. The percent will be accounted for in the rest
                    % of the Coeffs (the Coeffs will sum to 1 - n(zero)
                    continue
                elseif n1Cond > sum(n1Counts == n1Cond)
                    num = length(n1groupIdx);
                    groupRatio = (num/length(n1Counts));
                    percent = groupRatio * 100;

                    % Document how much data is omitted
                    note = strcat('Omitting ', num2str(percent), 'percent of cycles, where n1Cond = ' , num2str(n1Cond), 'because more spikes than data.');                    
                    disp(note)
                    obj.notes{noteCount,1} = note;
                    
                    % Keep track of total omitted ratio.
                    omitCoeff(noteCount) = groupRatio;

                    % Increase note counter
                    noteCount = noteCount + 1;
                    continue
                end
                % Define x data for iCond
                n1Group = n1(n1groupIdx,1:n1Cond);
                n2Group = n2Counts(n1groupIdx)';
                xGroup = [n1Group,n2Group];
                xGroups{iGroup,1} = xGroup;

                % Define y data for iCond
                yGroup = y(n1groupIdx,1:end);
                yGroups{iGroup,1} = yGroup;

                % Find coeff corresponding to iCond
                coeffs{iGroup,1} = length(n1groupIdx)/length(n1Counts);
                % Increase group counter
                iGroup = iGroup + 1;

                
            end
            % Audit: Check that omit coeffs and group coeffs sum to 1 with a very small tolerance to account for matlab rounding error. 
            if ~ismembertol((sum(cell2mat(coeffs)) + sum(omitCoeff)), 1, 1e-12); error('Error: Sum of coeffs and omitted data ratios does not equal 1'); end
            
            % Call parent class buildMIs()
            buildMIs@mi_analysis(obj, {xGroups yGroups coeffs});          
            end

            % From here, each entry in xGroups and yGroups will feed into
            % the MI calculator. 
            % Figure out how each subgroup is going to feed into the 
            % MI_sim_manager and set up the data for that (maybe via
            % different lists). 
            
            
        end
    end

