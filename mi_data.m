classdef mi_data < handle
    %  MI_KSG_data is used to set up a data object with all of the data
    %  for a given recording session
    % WE MAY WANT TO CHANGE PROPERTY NAMES TO GENERIC VARIABLES
    properties
        neurons % cell array
                % 1 x N array of spike timing vectors- spike times in MS. 
                % where N is the number of neurons recorded during this session. 
        n_timebase % either 'phase' or 'time'
%         behavior % N x n vector of the continuous pressure where N is the total cycles and n is the maximum 
                 % number of samples per cycle
        obj_behav % behavioral data object
%         cycleTimes % {1 x 2} array with 1: 1 x N vector of onset times in seconds of each breath cycle
                    % and 2: 1 x N vector of the peaks corresponding to the breathcycles
                    % where N is the total number of breath cycles 
        b_timebase % either 'phase' or 'time' 
                    % DEFAULT: 'time'
        b_length % An integer to indicate number of points to keep for each behavioral cycle
                    % DEFAULT: 11
        b_windowOfInterest % Either a phase window or a time window in milliseconds regardless of timebase
                    % This value determines what window of the data we want
                    % to use for our calculation. 
                    % DEFAULT: pi or 100ms
        b_startPhase % either a radian angle or a time in ms indicating to relativel time/phase
                   % to document the behavior relative to the cycle onset
                   % time. This must be in the same units as b_timebase
                   % DEFAULT: .8pi or 50ms
        b_dataTransform % either 'none, 'pca', or 'residual' - THIS CAN BE ADDED TO       
        bFs % sample frequency of pressure wave
        nFs % sample frequency of neural data 
        
        reparamData % 0 if we don't want to reparameterize our data to have a Gaussian distribution and 1 if we do.
        
        dataInfo % either a reference to the generator object or to a data file/name
        
        verbose % level of output for progress and troubleshooting/debugging
    end

    methods
        function obj = mi_data(nFs,pFs)
%            % This function documents the sample frequencies
%            % Note that I need to add the proper functions once Bryce sends me his code  
%           
%            % Initiate input parser
%            p = inputParser;
%            
%            % Set up required inputs
%            p.addRequired('nFs');
%            p.addRequired('bFs');
%            p.addRequired('data_ref');
%            
%            % Set up defaults and optional input for obj.n_timebase
%            default_n_timebase = 'time';
%            validate_n_timebase = @(x) assert(ismember(x,{'time','phase'}),'n_timebase must be either phase, time');
%            
%            % Default and optional input for obj.reparamData
%            default_reparamData = false;
%            validate_reparamData = @(x) isboolean(x);
%            
%            % Default and optional input for obj.b_timebase
%            default_b_timebase = 'phase';
%            validate_b_timebase = @(x) assert(ismember(x,{'phase','time'}), 'timebase must be either phase or time');
%            
%            % Default and optional input for obj.b_Length
%            default_b_Length = 11;
%            validate_b_Length = @(x) assert(isinteger(x),'length must be an integer value');
%             
%            % Default and optional input for obj.b_dataTransform
%            default_b_dataTransform = 'residual'; 
%            validate_b_dataTransform = @(x) assert(ismember(x,{'none','residual','pca'}), 'dataTransform must be none, residual, or pca');
% 
%             
%            % Set up the optional input: timebase
%            p.addParameter('n_timebase',default_n_timebase,validate_n_timebase);
%            % Set up optional input: reparamData
%            p.addParameter('reparamData', default_reparamData, validate_reparamData);
%            % Set up optional input: b_timebase
%            p.addParameter('b_timebase',default_b_timebase,validate_b_timebase);
%            % Set up optional input: b_Length
%            p.addParameter('b_Length',default_b_Length, validate_b_Length);
%            % Set up optional input: b_dataTransform
%            p.addParameter('b_dataTransform',default_b_dataTransform, validate_b_dataTransform);
%            
%            % Parse the inputs
%            p.parse(nFs, bFs, data_ref, varargin{:});
%            
%            obj.bFs = p.Results.bFs;
%            obj.nFs = p.Results.nFs;
%            obj.data_ref = p.Results.data_ref;
%            obj.n_timebase = p.Results.n_timebase;
%            obj.reparamData = p.Results.reparamData;
%            obj.b_timebase = p.Results.b_timebase;
%            obj.b_Length = p.Results.b_Length;
%            obj.b_dataTransform = p.Results.b_dataTransform;
%            
%            % Set other properties to empty arrays temporarily 
%            obj.neurons = {};
%            obj.behavior = {};
%            obj.cycleTimes = {};
%            obj.b_windowOfInterest = {};
%            obj.b_startPhase = {};
        end
        
        function add_data()
            
        end
        
        function get_data()

        end
       
        function r = reparam_data(obj,data)
            r = reparameterize_data(data);
        end
   
    end
end

    
