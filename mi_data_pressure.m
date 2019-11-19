classdef mi_data_pressure < mi_data_behavior
    %    
    % Class defined to load and process data from Intan RHD files for
    % air pressure assuming that the air pressure is recorded in the
    % board_adc channel
    %
   properties
       
   end
   methods
       function add_cycles()
           
       end
       
        function set_source(obj, obj_files)
%             %% takes a struct array of data files
%             obj.raw_source = obj_files; 
        end
        
        function calc_cycles(obj, method)
%             %% take raw pressure data, segment into cycles, save MAT file, calc onset times
%             
%             % Need to add Parser input to validate method
%             obj.cycleMethod = method;
%             
%             cycle_wavs = [];
%             for i=1:length(obj.raw_source)
%                 if obj.verbose > 0
%                     disp('===== ===== ===== ===== =====');
%                     disp(['Processing file ' num2str(i) ' of ' num2str(length(obj.raw_source))]);
%                     disp(['File: ' obj.raw_source(i).folder '\' obj.raw_source(i).name]);
%                 end
%                 [pressure_ts pressure_wav] = read_Intan_RHD2000_nongui_adc([obj.raw_source(i).folder '\' obj.raw_source(i).name], obj.verbose);
%                 
%                 % use RC's low pass filter
%                 % Convolve the gaussian window with the behavior data
%                 behav_smooth = smoothdata(pressure_wav, 'gaussian', 1024); % provides best local smoothing by trial and error
% 
%                 cycle_ixs = [];
%                 if strcmp(obj.cycleMethod, 'threshold');
%                     cycle_thresh = (max(behav_smooth) - min(behav_smooth))*0.25 + min(behav_smooth); % use 25% of waveform amplitude to trigger cycle onsets
%                     cycle_ixs = find(diff(sign(behav_smooth-cycle_thresh)) == 2);
%                     cycle_ixs = cycle_ixs(2:end-1); % automatically drop first and last due to impartial cycles
% 
%                     % scrub cycles for too short or too long
%                     cycle_len = diff(cycle_ixs);
%                     
%                     bad_cycle = find(isoutlier(cycle_len)==1);
%                     if obj.verbose > 1; disp(['--> Dropping ' num2str(length(bad_cycle)) ' cycles!']); end
%                     
%                     if ~isempty(bad_cycle)
%                         cycle_ixs(bad_cycle) = [];
%                         cycle_len(bad_cycle) = [];
%                     end
%                     
%                     % plot resulting behavior cycles to visualize analysis
%                     if obj.verbose > 2
%                         figure();
%                         plot(pressure_ts, pressure_wav, 'k-', 'LineWidth', 2);
%                         hold on;
%                         plot(pressure_ts, behav_smooth, 'g-');
%                         plot(pressure_ts([1 end]), [cycle_thresh cycle_thresh], 'r-');
%                         scatter(pressure_ts(cycle_ixs), behav_smooth(cycle_ixs), 80, 'r', 'filled');
%                     end
%                 end
%             
%                 % populate data structures
%                 tmp_cycle_wavs = nan(length(cycle_ixs), max(cycle_len));   
%                 
%                 % matrix for current data file
%                 for j=1:length(cycle_ixs)
%                     if j < length(cycle_ixs)
%                         wav_len = min(cycle_ixs(j+1) - cycle_ixs(j), size(tmp_cycle_wavs,2)); % pull waveform duration matching length of dataset
%                     else
%                         wav_len = min(size(pressure_wav,2) - cycle_ixs(j), size(tmp_cycle_wavs,2));
%                     end
%                    tmp_cycle_wavs(j,1:wav_len) = pressure_wav(cycle_ixs(j):cycle_ixs(j)+wav_len-1);
%                 end
%                 
%                 % add data to full dataset
%                 obj.cycleTimes = horzcat(obj.cycleTimes, pressure_ts(cycle_ixs));
%             
%                 if size(cycle_wavs,2) > size(tmp_cycle_wavs,2)
%                     nCols = size(cycle_wavs,2)-size(tmp_cycle_wavs,2);
%                     tmp_cycle_wavs(:,end+1:end+nCols) = nan(size(tmp_cycle_wavs,1), nCols);
%                 elseif size(cycle_wavs,2) < size(tmp_cycle_wavs,2)
%                     nCols = size(tmp_cycle_wavs,2)-size(cycle_wavs,2);
%                     cycle_wavs(:,end+1:end+nCols) = nan(size(cycle_wavs,1), nCols);
%                 end
%                 cycle_wavs(end+1:end+size(tmp_cycle_wavs,1),:) = tmp_cycle_wavs;
%                 
%                 disp('');
%             
%             end
%             
%             % save data in object
%             obj. cycleData = cycle_wavs;
        end
        
%        function set_behavior(obj, behavior, cycleFreq, cutoffFreq, filterFreq, varargin)
%         function set_behavior(obj, varargin)
%             p = inputParser;
%             
%             % BC 20190821: add behavior object and overwrite existing
% %             p.addRequired('behavior');
%             
%             
% %             if isempty(obj.behavior)
% %                % Behavior input is required if it has not been defined yet.
% %                p.addRequired('behavior');
% %             elseif ~isempty(obj.behavior)
% %                % If the behavior is already specified for the object, then
% %                % it defaults to the pre-set value, and the input is not
% %                % necessary.
% %                p.addOptional('behavior',obj.behavior);
% %             end
% 
%             % Convert behavior to cycles if its not already in cycles
%             % --> BC 20190820: this is taken care of by pressure cycle
%             % class
% %             if size(behavior,1) == 1
% % 
% %                 % Set optional inputs relevant for obtaining cycle times
% %                 % --> BC 20190710: add an object parameter that is a struct of arguments used to analyze behavior cycles
% %                 p.addOptional('cycleFreq',2);
% %                 p.addOptional('cutoffFreq',10);
% %                 p.addOptional('filterFreq', 100)
% %                 parse(p,behavior);            
% % % --------------Consider making this into a separate function-----------
% %                 % Add cycle times to object. 
% %                 obj.cycleTimes = obj.make_cycleTimes(behavior, p.Results.cycleFreq, p.Results.cutoffFreq);
% % 
% %                % Find cycle onset times
% %                cycle_times = obj.cycleTimes{1,1};
% %                % Convert from onset times in seconds to samples
% %                cycle_samples = ceil(cycle_times .* obj.bFs);
% %                % Find the length of each cycle
% %                cycle_lengths = diff(cycle_samples);
% % 
% %                % Find maximum cycle length
% %                maxLength = max(cycle_lengths);
% % 
% %                % Find total number of cycles
% %                nCycles = length(cycle_lengths);
% % 
% %                % Make an NaN matrix to hold cycle data. 
% %                behaviorCycles = nan(nCycles,maxLength);
% % 
% %                % Filter Pressure Waves
% %                filterData = obj.filterBehavior(behavior, p.Results.cycleFreq,p.Results.filterFreq);
% % 
% %                 % Assign pressure waves to matrix rows
% %                 for iCycle = 1:nCycles
% %                     behaviorCycles(iCycle,1:cycle_lengths(iCycle)) = filterData(cycle_samples(iCycle):(cycle_samples(iCycle+1)-1));
% %                 end
% % %-----------------------------------------------------------------------
% %                 % Store behavior
% %                 obj.behavior = behaviorCycles;
% %             end
% %             obj.obj_behav = behavior;
% 
%             % Set behavioral property defaults
%            
%             
%             default_b_timebase = 'phase';
%             validate_b_timebase = @(x) assert(ismember(x,{'phase','time'}), 'timebase must be either phase or time');
%             p.addParameter('timebase', default_b_timebase, validate_b_timebase);
%             p.KeepUnmatched = true; % allows for arguments that are not yet added to parser
%             p.parse(varargin{:});
%             
%             if isempty(obj.b_timebase)
%                 timebase = p.Results.timebase;
%             else
%                 if sum(strcmp(varargin, 'timebase')) == 1
%                     timebase = p.Results.timebase;
%                 else
%                     timebase = obj.b_timebase;
%                 end
%             end
%             
%             default_b_length = 11;
%             validate_b_length = @(x) assert(isnumeric(x) && x >= 0 && floor(x) == x,'length must be an integer value');
%             p.addParameter('length', default_b_length, validate_b_length);
%             
%             default_b_dataTransform = 'residual'; 
%             validate_b_dataTransform = @(x) assert(ismember(x,{'none','residual','pca'}), 'dataTransform must be none, residual, or pca');
%             p.addParameter('dataTransform', default_b_dataTransform, validate_b_dataTransform);
% 
% %             % Set behavior timebase value
% %             if isempty(obj.b_timebase)
% %                 % set n_timebase to default or inputed value
% %                 p.addParameter('timebase',default_b_timebase, validate_b_timebase);
% %             elseif ~isempty(obj.b_timebase)
% %                 % overwrite current n_timebase setting only if the value is
% %                 % inputted
% %                 p.addParameter('timebase',obj.b_timebase, validate_b_timebase);
% %             end
% %             p.parse(varargin{:});
% % %             p.parse(behavior,varargin{:});
% %             % Set behavior timebase property
% %             obj.b_timebase = p.Results.timebase;
% 
% %             % Set behavior length 
% %             if isempty(obj.b_length)
% %                 % set b_length to default or inputed value
% %                 p.addParameter('Length',default_b_length, validate_b_length);
% %             elseif ~isempty(obj.b_length)
% %                 % overwrite current b_length setting only if the value is
% %                 % inputted
% %                 p.addParameter('Length',obj.b_length, validate_b_length);
% %             end
% 
% 
% %             % Set behavior Length property
% % %             p.parse(behavior,varargin{:});
% %             p.parse(varargin{:});
% %             obj.b_length = p.Results.Length;
% 
%             % Set behavior startPhase 
%             % Adjust startPhase default and validation depending on timebase property
% %             if strcmp(obj.b_timebase, 'phase')
%             if strcmp(timebase, 'phase')
%                 default_b_startPhase = .8*pi;
%                 validate_b_startPhase = @(x) assert(0 <= x && x < 2*pi,'startPhase units must be in radians to match timebase');
%             elseif strcmp(timebase, 'time')
%                 default_b_startPhase = 50;
%                 validate_b_startPhase = @(x) assert(isnumeric(x) && x >= 0 && floor(x) == x,'startPhase units must be a positive integer in milliseconds to match timebase');
%                 % Can't use isinteger() because varargin automatically
%                 % casts argument as float, will deal with at end of
%                 % function during value assignments
%             end
%             p.addParameter('startPhase', default_b_startPhase, validate_b_startPhase);
%             
% %             % Set behavior startPhase
% %             if isempty(obj.b_startPhase)
% %                 % set n_startPhase to default or inputed value
% %                 p.addParameter('startPhase',default_b_startPhase, validate_b_startPhase);
% %             elseif ~isempty(obj.b_startPhase)
% %                 % overwrite current b_startPhase setting only if the value is
% %                 % inputted
% %                 p.addParameter('startPhase',obj.b_startPhase, validate_b_startPhase);
% %             end
%             % Set behavior startPhase property
% %             p.parse(behavior,varargin{:});
% %             p.parse(varargin{:});
% %             obj.b_startPhase = p.Results.startPhase;
% 
%             % Set behavior windowOfInterest
% 
%             % Adjust startPhase default and validation depending on timebase property
%             if strcmp(timebase, 'phase')
%                 default_b_windowOfInterest = pi;
%                 validate_b_windowOfInterest = @(x) assert(isnumeric(x),'startPhase must be numerical'); % check for valid time window at end of function
%             elseif strcmp(timebase, 'time')
%                 default_b_windowOfInterest = 100;
%                 validate_b_windowOfInterest = @(x) assert(isnumeric(x) && x >= 0,'windowOfInterest units must be a positive integer in milliseconds to match timebase');
%             end
%             p.addParameter('windowOfInterest', default_b_windowOfInterest, validate_b_windowOfInterest);
% %             % Set behavior windowOfInterest
% %             if isempty(obj.b_windowOfInterest)
% %                 % set b_windowOfInterest to default or inputed value
% %                 p.addParameter('windowOfInterest',default_b_windowOfInterest, validate_b_windowOfInterest);
% %             elseif ~isempty(obj.b_windowOfInterest)
% %                 % overwrite current b_windowOfInterest setting only if the value is
% %                 % inputted
% %                 p.addParameter('windowOfInterest',obj.b_windowOfInterest, validate_b_windowOfInterest);
% %             end
% %             % Set behavior startPhase property
% % %             p.parse(behavior,varargin{:});
% %             p.parse(varargin{:});
% %             obj.b_windowOfInterest = p.Results.windowOfInterest;
% 
% %             % Set behavior dataTransform
% %             if isempty(obj.b_dataTransform)
% %                 % set n_dataTransform to default or inputed value
% %                 p.addParameter('dataTransform',default_b_dataTransform, validate_b_dataTransform);
% %             elseif ~isempty(obj.b_dataTransform)
% %                 % overwrite current b_dataTransform setting only if the value is
% %                 % inputted
% %                 p.addParameter('dataTransform',obj.b_dataTransform, validate_b_dataTransform);
% %             end
% %             % Set behavior dataTransform property
% % %             p.parse(behavior,varargin{:});
% %             p.parse(varargin{:});
% %             obj.b_dataTransform = p.Results.dataTransform;
% 
%             p.parse(varargin{:});
% 
%             % NEED TO CHECK FOR VALIDITIY OF PRESSURE WAVE WITHIN SAME
%             % CYCLE OF SPIKES
% %             if startPhase + windowOfInterest > 2*pi
% %                 warning('Invalid limits: windowOfInterest plus startPhase must be less than 2*pi');
% %             end
%             
%             obj.b_timebase = p.Results.timebase;
%             obj.b_length = p.Results.length;            
%             obj.b_startPhase = double(p.Results.startPhase);
%             obj.b_windowOfInterest = p.Results.windowOfInterest;
%             obj.b_dataTransform = p.Results.dataTransform;
%             
%        end
       
       
%        function [filterData] = filterBehavior(obj, behavior, cycleFreq, filterFreq)
%            % This function prepares the raw behavioral data for analysis
%            % Convert cycle freq to length of gaussian in samples
%            cycleLengthSeconds = 1/cycleFreq;
%            cycleLengthSamples = cycleLengthSeconds * obj.bFs;
%            % Convert filter freq to width of gaussian in samples
%            filterWinSeconds = 1/filterFreq;
%            filterWinSamples = filterWinSeconds * obj.bFs;
%            
%            % Find alpha value for input to gaussian window function.
%            alpha = (cycleLengthSamples - 1)/(2*filterWinSamples);
%            
%            % Generate the gaussian window for filter
%            g = gausswin(cycleLengthSamples, alpha);
%            
%            filterData = conv(behavior,g,'same');
% 
%        end        
        
        function r = get_feature(obj, format)
%             %% Implemented to return different features of behavior cycles
%             % i.e., raw data, PCA, residual, area
%            
%             switch(format)
%                 case('raw')
%                     % need to implement from RC's code
%                 case('phase')
%                     % need to implement from RC's code                    
%             end
            
        end       
   end    
end