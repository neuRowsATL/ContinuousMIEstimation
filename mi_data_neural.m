classdef mi_data_neural < mi_data
    methods
        function obj = mi_data_neural(ID, varargin)
            % Required arguments: ID
            obj@mi_data(ID,varargin{:});
        end
        
        function add_spikes(obj, data, dataInfo, Fs, varargin)
            add_data(obj, data, dataInfo, Fs, varargin{:});
        end
        
        function add_cycleTimes(obj, data, dataInfo, Fs, varargin)
            add_data(obj, data, dataInfo, Fs, 'cycleTimes', varargin{:});
        end
        
        function r = get_spikes(obj, varargin)
            v = obj.verbose;
            p = inputParser;
            
            % Parameter format
            default_format = 'raw';
            valid_formats = {'raw', 'count', 'timing'};
            validate_format = @(x) assert(ischar(x) && ismember(x, valid_formats), 'format must be: raw, count, timing');
            p.addParameter('format', default_format, validate_format);
            
            % Parameter cycleTimes
            default_cycleTimes = [];
            validate_cycleTimes = @(x) assert(isnumeric(x) && size(x,2) == 2, 'cycleTimes must be N x 2 numeric');
            p.addParameter('cycleTimes', default_cycleTimes, validate_cycleTimes);
            
            % Parameter timebase
            default_timebase = 'time';
            valid_timebases = {'time', 'phase'};
            validate_timebase = @(x) assert(ischar(x) && ismember(x,valid_timebases), 'timebase must be: time, phase');
            p.addParameter('timeBase', default_timebase, validate_timebase);
            
            % Parameter name
            default_name = 'noname';
            validate_name = @(x) assert(ischar(x), 'name must be string/char');
            p.addParameter('name', default_name, validate_name);
            
            p.parse(varargin{:});
            format = p.Results.format;
            cycleTimes = p.Results.cycleTimes;
            timeBase = p.Results.timeBase;
            name = p.Results.name;
            
            
            if v>1; disp([newline '--> Getting spike data...']); end
            
            switch format
                case 'raw'
                    if v>1; disp('--> Raw data'); end
                    r = get_data(obj, name);
                case 'count'
                    r = get_count(obj, cycleTimes, name);
                    if v>1; disp('--> Spike count data'); end
                case 'timing'
                    r = get_timing(obj, cycleTimes, 'timeBase', timeBase, 'name', name);
                    if v>1; disp('--> Spike timing data'); end
            end
            if v>0; disp('COMPLETE: Spike data retrieved!'); end
        end
        
        function r = get_count(obj, cycleTimes, varargin)
            v = obj.verbose;
            
            p = inputParser;

            % Required cycleTimes
            validate_cycleTimes = @(x) assert(isnumeric(x) && size(x,2) == 2, 'cycleTimes must be N x 2 numeric');
            p.addRequired('cycleTimes', validate_cycleTimes);

            % Parameter name
            default_name = 'noname';
            validate_name = @(x) assert(ischar(x), 'name must be string/char');
            p.addOptional('name', default_name, validate_name);
            
            p.parse(cycleTimes, varargin{:});
            cycle_ts = p.Results.cycleTimes;
            name = p.Results.name;

            spike_ts = obj.data.(name).data;

            
            if v>1; disp([newline '--> Calculating spike count...']); end
            
            if v>2
                disp(['cycleTimes: ' regexprep(num2str(size(cycle_ts)), '\s*', ' x ')]);
                disp(['name: ' name]);
            end
            
            % Find the number of spikes in each cycle
            % We include data that comes after the onset of the first cycle
            % and before the onset of the last cycle
            cycle_spike_counts = zeros(1,size(cycle_ts,1)-1);
            
            % AVG run time: 0.4599 s
            for cycle_ix = 1:size(cycle_ts,1)
                cycle_spike_counts(cycle_ix) = sum((spike_ts > cycle_ts(cycle_ix,1)) & (spike_ts < cycle_ts(cycle_ix,2)));
            end
            
            
            r = cycle_spike_counts;
            
            if v>0; disp('COMPLETE: Spike count calculated!'); end
        end
        
        function r = get_timing(obj, cycleTimes, varargin)
            % Makes matrices of pressure data  need to decide what units to use and spiking data based on what we have      
            % NOTE- currently as the code is written, we omit any neural or pressure data that occurs	 
            % before the onset of the first cycle or after the onset of the last cycle
            %  - additionally, we are segmenting spikes based on the cycle times rather than trying
            % to keep bursts together and using negative spike times 
            v = obj.verbose; 

            p = inputParser;

            % Required cycleTimes
            validate_cycleTimes = @(x) assert(isnumeric(x) && size(x,2) == 2, 'cycleTimes must be N x 2 numeric');
            p.addRequired('cycleTimes', validate_cycleTimes);

            % Parameter timebase
            default_timebase = 'time';
            valid_timebases = {'time', 'phase'};
            validate_timebase = @(x) assert(ischar(x) && ismember(x,valid_timebases), 'timebase must be: time, phase');
            p.addParameter('timeBase', default_timebase, validate_timebase);

            % Parameter name
            default_name = 'noname';
            validate_name = @(x) assert(ischar(x), 'name must be string/char');
            p.addParameter('name', default_name, validate_name);

            p.parse(cycleTimes, varargin{:});
            cycle_ts = p.Results.cycleTimes;
            name = p.Results.name;
            timeBase = p.Results.timeBase;

            spike_ts = obj.data.(name).data;

            
            if v>1; disp([newline '--> Calculating spike timing...']); end
            
            if v>2
                disp(['cycleTimes: ' regexprep(num2str(size(cycle_ts)), '\s*', ' x ')]);
                disp(['timeBase: ' timeBase]);
                disp(['name: ' name]);
            end 
            
            % Find the number of spikes in each cycle
            if v>1; disp('--> Getting spike count...'); end
            cycle_spike_counts = obj.get_count(cycle_ts, name);

            % Calculate relative spike times for each breathing cycle
            % if verbose > 1; disp('-> Calculating relative spike times by cycle'); end
            if v>1; disp('--> Calculating timing based on count...'); end
            cycle_spike_ts = nan(size(cycle_ts,1), max(cycle_spike_counts));

            
            % AVG run time: 0.9658 s
            for cycle_ix = 1:(size(cycle_ts,1))
               cycle_spikes_ix = find((spike_ts > cycle_ts(cycle_ix,1)) & (spike_ts < cycle_ts(cycle_ix,2)));
               if ~isempty(cycle_spikes_ix)
                   cycle_spike_ts(cycle_ix,1:length(cycle_spikes_ix)) = spike_ts(cycle_spikes_ix) - cycle_ts(cycle_ix,1);
               end
            end
            
            
            switch timeBase
               case 'phase'
                   if v>1; disp('--> Transforming timing to phase...'); end
                   % Convert spike times to phase values in radians. 
                   cycle_durs = diff(cycle_ts');
                   % Calculate the phase conversion for each cycle
                   phase_factor = (2*pi)./cycle_durs;
                   % Propogate the phase conversion to a matrix
                   phase_factor_matrix = repmat(phase_factor,size(cycle_spike_ts,2),1);
                   % Multiply each spike time by the phase factor for the
                   % respective cycle
                   cycle_spike_phase = cycle_spike_ts.*phase_factor_matrix';
                   % Output the variable. 
                   r = cycle_spike_phase;

               case 'time'
                   r = cycle_spike_ts;
            end
            
            if v>0; disp('COMPLETE: Spike timing calculated'); end
        end
    end
end