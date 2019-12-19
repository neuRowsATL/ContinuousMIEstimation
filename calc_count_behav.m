classdef calc_count_behav < mi_analysis
    %Each of these objects sets the stage to calculate the mutual
    %information between spike count and behavior and stores the results of
    %the calculation. 
    %   Detailed explanation goes here
    
    properties
        b_timeBase
        feature
        start
        dur
        nSamp
        nPC
    end
    
    methods


       function obj = calc_count_behav(objData, objBehav, varNames, varargin)
       % Required arguments: objData, varNames
       % Optional, string-style inputs: b_timebase, feature, start, dur, nSamp
            % Check required inputs for validity using input parser
 
            % Set up input parser
            p = inputParser;
            
            % Set required inputs
            validate_objData = @(x) assert(isa(x, 'mi_data_neural'), 'objData must be a neural data subclass');
            p.addRequired('objData', validate_objData);
            
            validate_objBehav = @(x) assert(isa(x, 'mi_data_behavior'), 'objBehav must be a behavioral data subclass');
            p.addRequired('objBehav', validate_objBehav);
            
            validate_varNames = @(x) assert(iscell(x) && (length(x) == 1), 'varNames must be a cell of length 1');
            p.addRequired('varNames', validate_varNames);
                        
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
            
            % Define timebase property of subclass object
            obj.b_timeBase = p.Results.b_timeBase;
            obj.feature = p.Results.feature;
            obj.start = p.Results.start;
            obj.dur = p.Results.dur;
            obj.nSamp = p.Results.nSamp;
            obj.nPC = p.Results.nPC;

        end

        function buildMIs(obj)

        % Build the data and core objects necessary to run the sim manager for this analysis class.

            v = obj.verbose;

            % Find the total spike count in a cycle for neuron 1
            x_name = obj.varNames{1};
            x = obj.objData.get_spikes('name', x_name, 'format', 'count', 'cycleTimes', obj.objBehav.data.cycleTimes.data);

            % Audit Check
            if sum(x) ~= (sum(~isnan(obj.objData.data.(obj.varNames{1}).data)) - (sum(obj.objData.data.(obj.varNames{1}).data < obj.objBehav.data.cycleTimes.data(1,1) | obj.objData.data.(obj.varNames{1}).data > obj.objBehav.data.cycleTimes.data(end,2))))
                error('Error: Spike Counts for x do not match that expected from objData.varNames{1}.');
            end

            % Set the groups that will serve as the x variable.
            xGroups{1,1} = x;

            % Get the behavioral data for analysis
            y = get_behavior(obj.objBehav, obj.b_timeBase, obj.feature, obj.start, obj.dur, obj.nSamp,'nPC', obj.nPC );

            % Set the groups that will serve as the y variable.
            yGroups{1,1} = y;
            
            % For data with only one group, the coeff is 1. 
            
            coeffs = {1};

            buildMIs@mi_analysis(obj, {xGroups yGroups coeffs});          
            
            
        end
    end
end

