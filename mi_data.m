classdef mi_data < handle
    %  MI_KSG_data is used to set up a data object with all of the data
    %  for a given recording session
    % WE MAY WANT TO CHANGE PROPERTY NAMES TO GENERIC VARIABLES
    properties
        ID % unique identifier for data object
        
        data % array of data
        dataInfo % either a reference to the generator object or to a data file/name
        
        Fs % sampling rate in Hz
        
        verbose % level of output for progress and troubleshooting/debugging
    end

    methods
        function obj = mi_data(ID)
           % This function documents the sample frequencies
           % Note that I need to add the proper functions once Bryce sends me his code  
          
           % Initiate input parser
           p = inputParser;
           
           % Set up required inputs
           p.addRequired('ID');
           
           % Default and optional input for obj.b_Length
           default_verbose = 1;
           validate_verbose = @(x) assert(isinteger(x),'verbose must be an integer');
           % Set up optional input: verbose
           p.addParameter('verbose',default_verbose, validate_verbose);

           % Parse the inputs
           p.parse(ID, varargin{:});

           obj.data = []; % array of data
           obj.dataInfo = {}; % cell array of filename/generators
           obj.Fs = -1; % initially set to invalid value
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

    
