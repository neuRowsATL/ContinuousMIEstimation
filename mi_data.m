classdef mi_data < handle
    %  MI_KSG_data is used to set up a data object with all of the data
    %  for a given recording session
    % WE MAY WANT TO CHANGE PROPERTY NAMES TO GENERIC VARIABLES
    properties
        ID % unique identifier for data object
        
        data % array of data
        
        Fs % sampling rate in Hz
        
        verbose % level of output for progress and troubleshooting/debugging
    end

    methods
        function obj = mi_data(ID, varargin)
            % This function documents the sample frequencies
            % Note that I need to add the proper functions once Bryce sends me his code  

            % Initiate input parser
            p = inputParser;

            % Set up required inputs
            validate_ID = @(x) assert(ischar(x), 'ID must be a string/char');
            p.addRequired('ID');

            % Default and optional input for obj.b_Length
            default_verbose = 1;
            validate_verbose = @(x) assert(isnumeric(x) && rem(x,1) == 0,'verbose must be an integer');
            % Set up optional input: verbose
            p.addParameter('verbose',default_verbose, validate_verbose);

            % Parse the inputs
            p.parse(ID, varargin{:});

            obj.ID = p.Results.ID;
            obj.verbose = p.Results.verbose;

            obj.data = struct(); % struct of data
            obj.Fs = -1; % initially set to invalid value            
            
            if obj.verbose > 0; disp([newline 'mi_data instantiated: ', obj.ID]); end
        end
        
        function add_data(obj, data, dataInfo, Fs, varargin)
            % Add data from file

            v = obj.verbose;
            
            % Initiate input parser
            p = inputParser;
            
            % Required data file name
            validate_data = @(x) assert(ismatrix(x),'data must be a matrix');
            p.addRequired('data',validate_data);

            % Required data file name
            validate_dataInfo = @(x) assert(ischar(x),'dataInfo must be a string file name');
            p.addRequired('dataInfo',validate_dataInfo);

            % Required sampling frequency
            validate_Fs = @(x) assert(isnumeric(x) && rem(x,1) == 0 && x > 0,'Fs must be specified as sampling rate in Hz > 0');
            p.addRequired('Fs', validate_Fs);

            % Optional data name
            default_name = '';
            validate_name = @(x) assert(ischar(x) && isempty(str2num(x(1))), 'Name must be a string that does not start with a number');
            p.addOptional('name', default_name, validate_name);
            
            % Parse the inputs
            p.parse(data, dataInfo, Fs, varargin{:});
            
            if v>1; disp([newline '--> Adding data...']); end
            
            obj.Fs = p.Results.Fs;
            name = p.Results.name;
            
            if v>2; disp([newline 'Fs: ' num2str(obj.Fs) newline 'name: ' name]); end
            
            if isempty(name) % check to see if name is specified
                if v>3; disp('--> --> name is empty'); end
                if isempty(fields(obj.data)) % if mi_data has no data specified...
                    if v>3; disp('--> --> obj.data is empty'); end
                    
                    % Set obj.data and obj.dataInfo with default field  
                    obj.data.noname = struct();
                    obj.data.noname.data = p.Results.data;
                    obj.data.noname.info = p.Results.dataInfo;
                    
                    warning('No name specified, default field implemented.');
                    if v>2; disp([newline 'data: ' regexprep(num2str(size(obj.data.noname.data)), '\s*', ' x ') newline 'info: ' obj.data.noname.info]); end
                else
                    error('Name argument required if using multiple data matrices');
                end
            else
                if isfield(obj.data, 'noname') % check that the existing data has valid name
                    error('Existing data was not assigned a name. Must specify name for multiple data matrices.');
                else
                    if isfield(obj.data, name) % check that the specified name is valid
                        error('Name already exists in obj.data. Use a different name.');
                    end
                    
                    % assign data and dataInfo
                    obj.data.(name) = struct();
                    obj.data.(name).data = p.Results.data;
                    obj.data.(name).info = p.Results.dataInfo;
                    if v>2; disp([newline 'data: ' regexprep(num2str(size(obj.data.(name).data)), '\s*', ' x ') newline 'dataInfo: ' obj.data.(name).info]); end
                end
            end
            if v>0; disp(['COMPLETE: Added data to mi_data(' obj.ID ')']); end
            if v>1
                disp(['mi_data(' obj.ID ') with fields:']);
                fields(obj.data)
            end
            if v>2
                fs = fields(obj.data);
                for i=1:length(fs)
                    disp(obj.data.(fs{i}))
                end
            end
        end
        
        function r = get_data(obj, varargin)
            % Returns data in raw or processed form for analysis
            v = obj.verbose; 
            
            p = inputParser;
            
            % Optional name argument
            default_name = 'noname';
            validate_name = @(x) assert(ischar(x), 'Name must be a string/char');
            p.addOptional('name', default_name, validate_name);
            
            p.parse(varargin{:});
            
            name = p.Results.name;
            r = obj.data.(name).data; % if no name was specified when adding data, "noname" is default field
            if v>0; disp(['COMPLETE: Data returned']); end
        end
    end
end  
