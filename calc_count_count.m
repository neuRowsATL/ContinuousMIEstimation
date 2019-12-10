classdef calc_count_count < mi_analysis
    %Each of these objects sets the stage to calculate the mutual
    %information between spike count and behavior and stores the results of
    %the calculation. 
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods

       function obj = calc_count_count(objData,varNames, varargin)
            % Required arguments: objData, varNames
            % Check required inputs for validity using input parser
            
            % Set up input parser
            p = inputParser;
            
            % Set required inputs
            validate_objData = @(x) assert(isa(x, 'mi_data_neural'), 'objData must be a neural data subclass');
            p.addRequired('objData', validate_objData);
            
            validate_varNames = @(x) assert(iscell(x) && (length(x) == 2), 'varNames must be a cell of length 2');
            p.addRequired('varNames', validate_varNames);
            
            p.parse(objData, varNames);
            
            objData = p.Results.objData;
            varNames = p.Results.varNames;
            
            % Check that varNames references valid fields of objData
            for ivarNames = 1:length(varNames)
                assert(isfield(objData.data , varNames{ivarNames}), ['varName: ' varNames{ivarNames} 'is not a valid field of the neural data object']); 
            end
            
            % Call parent constructor
            obj@mi_analysis(objData, varNames, varargin{:});
        end
        
        function  buildMIs(obj)
            % Build the data and core objects necessary to run the sim manager for this analysis class. 


            % Find total spike count in a cycle for neuron 1 
            x_name  = obj.varNames{1};
            x = obj.objData.get_spikes('name', x_name , 'format', 'count', 'cycleTimes', obj.objData.data.cycleTimes.data );

            % Set groups that will serve as x variable
            % RC20191210: Nemenman code automatically adds noise. Removed
            % the extra noise
            xGroups{1,1} = x;

            % Next find spike count for neuron 2
            y_name = obj.varNames{2};
            y = obj.objData.get_spikes('name', y_name, 'format', 'count', 'cycleTimes', obj.objData.data.cycleTimes.data );

            % Set groups that will serve as y variable
            % RC20191210: Nemenman code automatically adds noise. Removed
            % the extra noise
            yGroups{1,1} = y;

            % Set coefficients for groups
            coeffs = {1};
            
            % Call parent buildMIs()
            buildMIs@mi_analysis(obj, {xGroups yGroups coeffs});
        end
    end
end

