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
           % This function documents the sample frequencies
           % Note that I need to add the proper functions once Bryce sends me his code  
          
           % Initiate input parser
           p = inputParser;
           
           % Set up required inputs
           p.addRequired('nFs');
           p.addRequired('bFs');
           p.addRequired('data_ref');
           
           % Set up defaults and optional input for obj.n_timebase
           default_n_timebase = 'time_ms';
           validate_n_timebase = @(x) assert(ismember(x,{'time_s','time_ms','phase'}),'n_timebase must be either phase, time_s, or time_ms');
           
           % Default and optional input for obj.reparamData
           default_reparamData = false;
           validate_reparamData = @(x) isboolean(x);
           
           % Default and optional input for obj.b_timebase
           default_b_timebase = 'phase';
           validate_b_timebase = @(x) assert(ismember(x,{'phase','time'}), 'timebase must be either phase or time');
           
           % Default and optional input for obj.b_Length
           default_b_Length = 11;
           validate_b_Length = @(x) assert(isinteger(x),'length must be an integer value');
            
           % Default and optional input for obj.b_dataTransform
           default_b_dataTransform = 'residual'; 
           validate_b_dataTransform = @(x) assert(ismember(x,{'none','residual','pca'}), 'dataTransform must be none, residual, or pca');

            
           % Set up the optional input: timebase
           p.addParameter('n_timebase',default_n_timebase,validate_n_timebase);
           % Set up optional input: reparamData
           p.addParameter('reparamData', default_reparamData, validate_reparamData);
           % Set up optional input: b_timebase
           p.addParameter('b_timebase',default_b_timebase,validate_b_timebase);
           % Set up optional input: b_Length
           p.addParameter('b_Length',default_b_Length, validate_b_Length);
           % Set up optional input: b_dataTransform
           p.addParameter('b_dataTransform',default_b_dataTransform, validate_b_dataTransform);
           
           % Parse the inputs
           p.parse(nFs, bFs, data_ref, varargin{:});
           
           obj.bFs = p.Results.bFs;
           obj.nFs = p.Results.nFs;
           obj.data_ref = p.Results.data_ref;
           obj.n_timebase = p.Results.n_timebase;
           obj.reparamData = p.Results.reparamData;
           obj.b_timebase = p.Results.b_timebase;
           obj.b_Length = p.Results.b_Length;
           obj.b_dataTransform = p.Results.b_dataTransform;
           
           % Set other properties to empty arrays temporarily 
           obj.neurons = {};
           obj.behavior = {};
           obj.cycleTimes = {};
           obj.b_windowOfInterest = {};
           obj.b_startPhase = {};
       end
       
       function add_spikes(obj, spike_times, varargin)                 
           % BC-20190123: Can add an index to neurons for explicit tracking
           % across classes?
           % RC: 20190517: Not sure what this would look like...

           % Set default optional parameter value
           default_n_timebase = 'time';
           validate_n_timebase = @(x) assert(ismember(x,{'time','phase'}),'timebase must be either phase or time');
           
           % Set up inputParser
           p = inputParser;
           % SpikeTimes are required if they haven't been specified
           if isempty(obj.neurons)
               p.addRequired('spike_times');
               p.parse(spike_times);
               
           elseif ~isempty(obj.neurons)
               p.addOptional('spike_times',[]);
               p.parse(spike_times);
           end

           % Add a set of spike times to object unless the spikes weren't
           % specified
           if ~isempty(spike_times)
               obj.neurons{end+1} = spike_times;
           end
           
           % Determine n_timebase has been set
           if isempty(obj.n_timebase)
               % If n_timebase has not been defined yet, set it to the
               % default or specified value
               p.addParameter('timebase',default_n_timebase,validate_n_timebase)
           elseif ~isempty(obj.n_timebase)
               % If n_timebase has been defined, the default is the
               % pre-assigned value, and it should only change if timebase
               % is an input. 
               p.addParameter('timebase',obj.n_timebase,validate_n_timebase);
           end
           p.parse(spike_times);
           obj.n_timebase = p.Results.timebase;

       end
       
       function add_behavior(obj, behavior)
           obj.obj_behav = behavior;
           obj.set_behavior();
       end
       
%        function set_behavior(obj, behavior, cycleFreq, cutoffFreq, filterFreq, varargin)
        function set_behavior(obj, varargin)
            p = inputParser;
            
            % BC 20190821: add behavior object and overwrite existing
%             p.addRequired('behavior');
            
            
%             if isempty(obj.behavior)
%                % Behavior input is required if it has not been defined yet.
%                p.addRequired('behavior');
%             elseif ~isempty(obj.behavior)
%                % If the behavior is already specified for the object, then
%                % it defaults to the pre-set value, and the input is not
%                % necessary.
%                p.addOptional('behavior',obj.behavior);
%             end

            % Convert behavior to cycles if its not already in cycles
            % --> BC 20190820: this is taken care of by pressure cycle
            % class
%             if size(behavior,1) == 1
% 
%                 % Set optional inputs relevant for obtaining cycle times
%                 % --> BC 20190710: add an object parameter that is a struct of arguments used to analyze behavior cycles
%                 p.addOptional('cycleFreq',2);
%                 p.addOptional('cutoffFreq',10);
%                 p.addOptional('filterFreq', 100)
%                 parse(p,behavior);            
% % --------------Consider making this into a separate function-----------
%                 % Add cycle times to object. 
%                 obj.cycleTimes = obj.make_cycleTimes(behavior, p.Results.cycleFreq, p.Results.cutoffFreq);
% 
%                % Find cycle onset times
%                cycle_times = obj.cycleTimes{1,1};
%                % Convert from onset times in seconds to samples
%                cycle_samples = ceil(cycle_times .* obj.bFs);
%                % Find the length of each cycle
%                cycle_lengths = diff(cycle_samples);
% 
%                % Find maximum cycle length
%                maxLength = max(cycle_lengths);
% 
%                % Find total number of cycles
%                nCycles = length(cycle_lengths);
% 
%                % Make an NaN matrix to hold cycle data. 
%                behaviorCycles = nan(nCycles,maxLength);
% 
%                % Filter Pressure Waves
%                filterData = obj.filterBehavior(behavior, p.Results.cycleFreq,p.Results.filterFreq);
% 
%                 % Assign pressure waves to matrix rows
%                 for iCycle = 1:nCycles
%                     behaviorCycles(iCycle,1:cycle_lengths(iCycle)) = filterData(cycle_samples(iCycle):(cycle_samples(iCycle+1)-1));
%                 end
% %-----------------------------------------------------------------------
%                 % Store behavior
%                 obj.behavior = behaviorCycles;
%             end
%             obj.obj_behav = behavior;

            % Set behavioral property defaults
           
            
            default_b_timebase = 'phase';
            validate_b_timebase = @(x) assert(ismember(x,{'phase','time'}), 'timebase must be either phase or time');
            p.addParameter('timebase', default_b_timebase, validate_b_timebase);
            p.KeepUnmatched = true; % allows for arguments that are not yet added to parser
            p.parse(varargin{:});
            
            if isempty(obj.b_timebase)
                timebase = p.Results.timebase;
            else
                if sum(strcmp(varargin, 'timebase')) == 1
                    timebase = p.Results.timebase;
                else
                    timebase = obj.b_timebase;
                end
            end
            
            default_b_length = 11;
            validate_b_length = @(x) assert(isnumeric(x) && x >= 0 && floor(x) == x,'length must be an integer value');
            p.addParameter('length', default_b_length, validate_b_length);
            
            default_b_dataTransform = 'residual'; 
            validate_b_dataTransform = @(x) assert(ismember(x,{'none','residual','pca'}), 'dataTransform must be none, residual, or pca');
            p.addParameter('dataTransform', default_b_dataTransform, validate_b_dataTransform);

%             % Set behavior timebase value
%             if isempty(obj.b_timebase)
%                 % set n_timebase to default or inputed value
%                 p.addParameter('timebase',default_b_timebase, validate_b_timebase);
%             elseif ~isempty(obj.b_timebase)
%                 % overwrite current n_timebase setting only if the value is
%                 % inputted
%                 p.addParameter('timebase',obj.b_timebase, validate_b_timebase);
%             end
%             p.parse(varargin{:});
% %             p.parse(behavior,varargin{:});
%             % Set behavior timebase property
%             obj.b_timebase = p.Results.timebase;

%             % Set behavior length 
%             if isempty(obj.b_length)
%                 % set b_length to default or inputed value
%                 p.addParameter('Length',default_b_length, validate_b_length);
%             elseif ~isempty(obj.b_length)
%                 % overwrite current b_length setting only if the value is
%                 % inputted
%                 p.addParameter('Length',obj.b_length, validate_b_length);
%             end


%             % Set behavior Length property
% %             p.parse(behavior,varargin{:});
%             p.parse(varargin{:});
%             obj.b_length = p.Results.Length;

            % Set behavior startPhase 
            % Adjust startPhase default and validation depending on timebase property
%             if strcmp(obj.b_timebase, 'phase')
            if strcmp(timebase, 'phase')
                default_b_startPhase = .8*pi;
                validate_b_startPhase = @(x) assert(0 <= x && x < 2*pi,'startPhase units must be in radians to match timebase');
            elseif strcmp(timebase, 'time')
                default_b_startPhase = 50;
                validate_b_startPhase = @(x) assert(isnumeric(x) && x >= 0 && floor(x) == x,'startPhase units must be a positive integer in milliseconds to match timebase');
                % Can't use isinteger() because varargin automatically
                % casts argument as float, will deal with at end of
                % function during value assignments
            end
            p.addParameter('startPhase', default_b_startPhase, validate_b_startPhase);
            
%             % Set behavior startPhase
%             if isempty(obj.b_startPhase)
%                 % set n_startPhase to default or inputed value
%                 p.addParameter('startPhase',default_b_startPhase, validate_b_startPhase);
%             elseif ~isempty(obj.b_startPhase)
%                 % overwrite current b_startPhase setting only if the value is
%                 % inputted
%                 p.addParameter('startPhase',obj.b_startPhase, validate_b_startPhase);
%             end
            % Set behavior startPhase property
%             p.parse(behavior,varargin{:});
%             p.parse(varargin{:});
%             obj.b_startPhase = p.Results.startPhase;

            % Set behavior windowOfInterest

            % Adjust startPhase default and validation depending on timebase property
            if strcmp(timebase, 'phase')
                default_b_windowOfInterest = pi;
                validate_b_windowOfInterest = @(x) assert(isnumeric(x),'startPhase must be numerical'); % check for valid time window at end of function
            elseif strcmp(timebase, 'time')
                default_b_windowOfInterest = 100;
                validate_b_windowOfInterest = @(x) assert(isnumeric(x) && x >= 0,'windowOfInterest units must be a positive integer in milliseconds to match timebase');
            end
            p.addParameter('windowOfInterest', default_b_windowOfInterest, validate_b_windowOfInterest);
%             % Set behavior windowOfInterest
%             if isempty(obj.b_windowOfInterest)
%                 % set b_windowOfInterest to default or inputed value
%                 p.addParameter('windowOfInterest',default_b_windowOfInterest, validate_b_windowOfInterest);
%             elseif ~isempty(obj.b_windowOfInterest)
%                 % overwrite current b_windowOfInterest setting only if the value is
%                 % inputted
%                 p.addParameter('windowOfInterest',obj.b_windowOfInterest, validate_b_windowOfInterest);
%             end
%             % Set behavior startPhase property
% %             p.parse(behavior,varargin{:});
%             p.parse(varargin{:});
%             obj.b_windowOfInterest = p.Results.windowOfInterest;

%             % Set behavior dataTransform
%             if isempty(obj.b_dataTransform)
%                 % set n_dataTransform to default or inputed value
%                 p.addParameter('dataTransform',default_b_dataTransform, validate_b_dataTransform);
%             elseif ~isempty(obj.b_dataTransform)
%                 % overwrite current b_dataTransform setting only if the value is
%                 % inputted
%                 p.addParameter('dataTransform',obj.b_dataTransform, validate_b_dataTransform);
%             end
%             % Set behavior dataTransform property
% %             p.parse(behavior,varargin{:});
%             p.parse(varargin{:});
%             obj.b_dataTransform = p.Results.dataTransform;

            p.parse(varargin{:});

            % NEED TO CHECK FOR VALIDITIY OF PRESSURE WAVE WITHIN SAME
            % CYCLE OF SPIKES
%             if startPhase + windowOfInterest > 2*pi
%                 warning('Invalid limits: windowOfInterest plus startPhase must be less than 2*pi');
%             end
            
            obj.b_timebase = p.Results.timebase;
            obj.b_length = p.Results.length;            
            obj.b_startPhase = double(p.Results.startPhase);
            obj.b_windowOfInterest = p.Results.windowOfInterest;
            obj.b_dataTransform = p.Results.dataTransform;
            
       end
       
%        function [Times] = make_cycleTimes(obj, behavior, cycleFreq, cutoffFreq)
%         function make_cycleTimes(obj)
%            % This function takes in raw behavioral data, applies a low pass filter
%            % and identifies the onset of cycle times based on the
%            % negative peaks. 
%            % NOTE- this function does not deal with cycles that need to be
%            % omitted. 
%            % Set sample Frequency
%            Fs = obj.bFs;
%            
%            % Convert cycle frequency to Samples
%            cycleLengthSeconds = 1/cycleFreq;
%            cycleLengthSamples = cycleLengthSeconds*Fs;           
% 
%            % Convert cutoff freq to width of gaussian in samples
%            cutoffSeconds = 1/cutoffFreq;
%            cutoffSamples = cutoffSeconds*Fs;
% 
%            % Find alpha value for input into gausswin Function
%            alpha = (cycleLengthSamples - 1)/(2*cutoffSamples);
%            
%            % Generate Gaussian Window
%            g = gausswin(cycleLengthSamples,alpha);
% 
%            % Convolve the gaussian window with the behavior data
%            behaviorSmoothed = conv(behavior,g,'same');
% 
%            % Find the negative peaks of the pressure cycles to determine the onset
%            % times
% 
%            % We use the negative pressure vector to find negative peaks
%            behaviorForPeaks = -1*behaviorSmoothed;
%            [pks, locs] = findpeaks(behaviorForPeaks,Fs, 'MinPeakDistance', cycleLengthSeconds/1.3);
%            
%            % BC 20190515: struct may come with additional overheadx
%            % RC 20190520: I don't remember how we wanted to change this...
%            Times = {locs,pks};
%            obj.cycleTimes = Times;

%             obj.cycleTimes = obj.obj_behav.cycleTimes; % reference times out of pressure data object

%        end
       
       function [filterData] = filterBehavior(obj, behavior, cycleFreq, filterFreq)
           % This function prepares the raw behavioral data for analysis
           % Convert cycle freq to length of gaussian in samples
           cycleLengthSeconds = 1/cycleFreq;
           cycleLengthSamples = cycleLengthSeconds * obj.bFs;
           % Convert filter freq to width of gaussian in samples
           filterWinSeconds = 1/filterFreq;
           filterWinSamples = filterWinSeconds * obj.bFs;
           
           % Find alpha value for input to gaussian window function.
           alpha = (cycleLengthSamples - 1)/(2*filterWinSamples);
           
           % Generate the gaussian window for filter
           g = gausswin(cycleLengthSamples, alpha);
           
           filterData = conv(behavior,g,'same');

       end

       function r = getTiming(obj, dataNum)
           % Makes matrices of pressure data  need to decide what units to use and spiking data based on what we have      
           % NOTE- currently as the code is written, we omit any neural or pressure data that occurs	 
           % before the onset of the first cycle or after the onset of the last cycle
           %  - additionally, we are segmenting spikes based on the cycle times rather than trying
           % to keep bursts together and using negative spike times 
           %
           % INPUT
           % dataNum : positive integer neuron number to specify neuron of interest
           %
             % Return an m x n matrix of data with m cycles and n sample points (if pressure) 
             % or n maximum spikes (if neuron data) 
           %
           % if verbose > 1; disp([newline 'Running: dataByCycles' newline]); end

           spike_ts = obj.neurons{dataNum{1}};
%            cycle_ts = obj.cycleTimes{1,1};
           cycle_ts = obj.obj_behav.cycleTimes'*1000.; % BC20190820

           % Find the number of spikes in each cycle
           cycle_spike_counts = obj.getCount(dataNum);

           % Calculate relative spike times for each breathing cycle
           % if verbose > 1; disp('-> Calculating relative spike times by cycle'); end
           cycle_spike_ts = nan(size(cycle_ts,1)-1, max(cycle_spike_counts));
           for cycle_ix = 1:(size(cycle_ts,1)-1)
               cycle_spikes_ix = find((spike_ts > cycle_ts(cycle_ix)) & (spike_ts < cycle_ts(cycle_ix+1)));
               if ~isempty(cycle_spikes_ix)
                   cycle_spike_ts(cycle_ix,1:length(cycle_spikes_ix)) = spike_ts(cycle_spikes_ix)-cycle_ts(cycle_ix);
                   disp('goodbye');
               end
               disp('goodbye2');
           end
           switch(obj.n_timebase)
               case('phase')
                   % Convert spike times to phase values in radians. 
                   cycle_lengths = diff(cycle_ts);
                   % Find the dimensions of the cycle_spike_ts matrix
                   dimension_cycle_spike_ts = size(cycle_spike_ts);
                   % Calculate the phase conversion for each cycle
                   phase_factor = (2*pi)./cycle_lengths;
                   % Propogate the phase conversion to a matrix
                   phase_factor_matrix = repmat(phase_factor,1,dimension_cycle_spike_ts(2));
                   % Multiply each spike time by the phase factor for the
                   % respective cycle
                   cycle_spike_phase = cycle_spike_ts.*phase_factor_matrix;
                   % Output the variable. 
                   r = cycle_spike_phase;
                   
               case('time')
                   r = cycle_spike_ts;
           end
           
           disp('hello');
       
       end
       
       function r = getCount(obj, dataNum)
           spike_ts = obj.neurons{dataNum{1}};
%            cycle_ts = obj.cycleTimes{1,1};
            cycle_ts = obj.obj_behav.cycleTimes'*1000.; % BC 20190820

           % Find the number of spikes in each cycle
           % We include data that comes after the onset of the first cycle
           % and before the onset of the last cycle
           cycle_spike_counts = zeros(1,size(cycle_ts,1)-1);
           for cycle_ix = 1:(size(cycle_ts,1)-1)
               cycle_spikes_ix = find((spike_ts > cycle_ts(cycle_ix)) & (spike_ts < cycle_ts(cycle_ix+1)));
               if ~isempty(cycle_spikes_ix)
                 cycle_spike_counts(cycle_ix) = length(cycle_spikes_ix);
               end
           end
           r = cycle_spike_counts;
       end
       
       
        function r = processBehavior(obj)
            
            behavData = obj.obj_behav.cycleData;
            
%             if obj.verbose > 2
%                 figure()
%                 h_avgBehav = plot(nanmean(behavData), 'k-', 'LineWidth', 3);
%                 hold on;
%             end
            
            
            switch(obj.b_dataTransform)
            case('pca')
                switch(obj.b_timebase)
                case('phase')
                    % Find the lengths of the cycles in samples
%                     cycleLengths_samples = sum(~isnan(obj.behavior),2);
                    cycleLengths_samples = sum(~isnan(behavData),2);

                    % Find the sample associated with the start phase for each
                    % cycle
                    start_samples = ceil(obj.b_startPhase.*(cycleLengths_samples./(2*pi)));

                    % Find the number of samples that encompases the window of
                    % interest for the cycles. 
                    windowOfInterest_samples = ceil(obj.b_windowOfInterest.*(cycleLengths_samples./(2*pi)));

                    % Find the stop sample of the window of interest for each
                    % cycle
                    stop_samples = start_samples + windowOfInterest_samples;

                   % Set up empty matrix to store pressure data.
                   nCycles = length(cycleLengths_samples);
                   cycle_behavior = nan(nCycles, 1000);

                   for cycle_ix = 1:nCycles
                       % Document all of the data points for the window of
                       % interest
%                        cycle_data = obj.behavior(cycle_ix,start_samples(cycle_ix):stop_samples(cycle_ix));
                       cycle_data = behavData(cycle_ix,start_samples(cycle_ix):stop_samples(cycle_ix));
                       % We resample the cyclic data to make the data have
                       % uniform length --> I am arbitrarily picking 1000
                       % data points. 
                       resampled_cycle_data = resample(cycle_data,1000,length(cycle_data));
                       cycle_behavior(cycle_ix, 1:1000) = resampled_cycle_data; 
                   end

                case('time')
                % Find the lengths of the cycles in samples
%                     cycleLengths_samples = sum(~isnan(obj.behavior),2);
                    cycleLengths_samples = sum(~isnan(behavData),2);

                    % Find the sample associated with the start time for each
                    % cycle
                    % Convert startPhase from ms to seconds
                    startTime = obj.b_startPhase/1000;
                    start_samples = ceil(startTime*obj.bFs);

                    % Find the number of samples that encompases the window of
                    % interest for the cycles. 
                    % Convert windowOFInterest from ms to seconds
                    windowOfInterest_seconds = obj.b_windowOfInterest./1000;
                    windowOfInterest_samples = ceil(windowOfInterest_seconds*obj.bFs);

                    % Find the stop sample of the window of interest for each
                    % cycle
                    stop_samples = start_samples + windowOfInterest_samples;

                   % Set up empty matrix to store pressure data.
                   nCycles = length(cycleLengths_samples);
                   cycle_behavior = nan(nCycles, windowOfInterest_samples);

                   for cycle_ix = 1:nCycles
                       % Document all of the data points for the window of
                       % interest
%                        cycle_data = obj.behavior(cycle_ix,start_samples:stop_samples); 
                       cycle_data = behavData(cycle_ix,start_samples:stop_samples); 
                       cycle_behavior(cycle_ix, 1:windowOfInterest_samples+1) = cycle_data;
                   end
                end
                % Find the PCs of the window of interest cycle data
                [~,score,~] = pca(cycle_behavior);
                cycle_dataPCs = score(:,1:obj.b_length);
                r = cycle_dataPCs;
                
                otherwise % 'none' or 'residual'
                    switch(obj.b_timebase)
                        case('phase')
                            % Find the lengths of the cycles in samples
        %                     cycleLengths_samples = sum(~isnan(obj.behavior),2);
                            cycleLengths_samples = sum(~isnan(behavData),2);

                            % Find the sample associated with the start phase for each
                            % cycle
                            start_samples = ceil(obj.b_startPhase.*(cycleLengths_samples./(2*pi)));

                            % Find the number of samples that encompases the window of
                            % interest for the cycles. 
                            windowOfInterest_samples = ceil(obj.b_windowOfInterest.*(cycleLengths_samples./(2*pi)));

                            % Find the stop sample of the window of interest for each
                            % cycle
                            stop_samples = start_samples + windowOfInterest_samples;

                           % Set up empty matrix to store pressure data.
                           nCycles = length(cycleLengths_samples);
                           cycle_behavior = nan(nCycles, obj.b_length);
                           
                           for cycle_ix = 1:nCycles
                               % Document all of the data points for the window of
                               % interest
        %                        cycle_data = obj.behavior(cycle_ix,start_samples(cycle_ix):stop_samples(cycle_ix));
                               cycle_data = behavData(cycle_ix,start_samples(cycle_ix):stop_samples(cycle_ix));

%                                if obj.verbose > 2; plot(start_samples(cycle_ix):stop_samples(cycle_ix), cycle_data); end
                               
                               % Resample to get only the desired number of points
                               resampled_cycle_data = resample(cycle_data,obj.b_length,length(cycle_data));
                               cycle_behavior(cycle_ix, 1:obj.b_length) = resampled_cycle_data;   
                           end

                        case('time')
                        % Find the lengths of the cycles in samples
        %                     cycleLengths_samples = sum(~isnan(obj.behavior),2);
                            cycleLengths_samples = sum(~isnan(behavData),2);

                            % Find the sample associated with the start time for each
                            % cycle
                            % Convert startPhase from ms to seconds
                            startTime = obj.b_startPhase/1000;
                            start_samples = ceil(startTime*obj.bFs);

                            % Find the number of samples that encompases the window of
                            % interest for the cycles. 
                            % Convert windowOFInterest from ms to seconds
                            windowOfInterest_seconds = obj.b_windowOfInterest/1000;
                            windowOfInterest_samples = ceil(windowOfInterest_seconds*obj.bFs);

                            % Find the stop sample of the window of interest for each
                            % cycle
                            stop_samples = start_samples + windowOfInterest_samples;

                           % Set up empty matrix to store pressure data.
                           nCycles = length(cycleLengths_samples);
                           cycle_behavior = nan(nCycles, obj.b_length);

                           figure();
                           
                           for cycle_ix = 1:nCycles
                               % Document all of the data points for the window of
                               % interest
        %                        cycle_data = obj.behavior(cycle_ix,start_samples:stop_samples);
                               cycle_data = behavData(cycle_ix,start_samples:stop_samples);
                               % Resample to get only the desired number of points
                               
                               plot(behavData(cycle_ix,:));
                               hold on;
                               plot(start_samples:stop_samples, cycle_data);
                               
                               resampled_cycle_data = resample(cycle_data,obj.b_length,length(cycle_data));
                               cycle_behavior(cycle_ix, 1:obj.b_length) = resampled_cycle_data;   
                           end
                    end

                    % Transform the behavioral data
                    switch(obj.b_dataTransform)
                        case('none')
                            r = cycle_behavior;
                        case('residual')
                            avg_cycle = mean(cycle_behavior,1);
                            avg_cycleMatrix = repmat(avg_cycle, length(cycle_behavior),1);
                            cycle_residuals = cycle_behavior - avg_cycleMatrix;
                            r = cycle_residuals;

                    end
            end

%             if obj.verbose > 2; uistack(h_avgBehav, 'top'); end

        end

        function r = reparameterizeData(obj,data)
            r = reparameterize_data(data);
        end
   
    end
end

    
