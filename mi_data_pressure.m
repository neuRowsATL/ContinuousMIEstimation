classdef mi_data_pressure < mi_data_behavior
    %    
    % Class defined to load and process data from Intan RHD files for
    % air pressure assuming that the air pressure is recorded in the
    % board_adc channel
    %
   properties
       
   end
   methods
       function obj = mi_data_pressure(ID, varargin)
            % Required arguments: ID
            obj@mi_data_behavior(ID,varargin{:}); 
       end
       
       function set_data_files(obj, arrDataFiles, varargin)
           % Input needs to be cell array of file names
           % Optional folder argument defaults to current working directory
           v = obj.verbose;
           
           p = inputParser;
           
           % Required input: arrDataFiles
           validate_arrDataFiles = @(x) assert(iscell(x) && ischar(x{1}), 'arrDataFiles needs to be cell array of string file names');
           p.addRequired('arrDataFiles', validate_arrDataFiles);
           
           % Optional input: strFldrPath
           default_strFldrPath = pwd;
           validate_strFldrPath = @(x) assert(true, 'strFldrPath needs to be a string path');
           p.addOptional('strFldrPath', default_strFldrPath, validate_strFldrPath);
           
           p.parse(arrDataFiles, varargin{:});

           if v>1; disp([newline '--> Setting data file info...']); end
           
           if isfolder(p.Results.strFldrPath)
               if v>2; disp(['--> --> Folder exists: ' p.Results.strFldrPath]); end
               filepaths = fullfile(p.Results.strFldrPath, p.Results.arrDataFiles);
               for i=1:length(filepaths)
                   if ~isfile(filepaths(i)); error(['File does not exist: ' filepaths(i)]); end
                   if v>2 && isfile(filepaths(i)); disp(['File: ' filepaths(i)]); end
               end
                obj.arrFiles = p.Results.arrDataFiles;
                obj.strFldr = p.Results.strFldrPath;
           else
               error(['strFldrPath does not exist: ' p.Results.strFldrPath]);
           end
            if v>0; disp('COMPLETE: Data files set!'); end
       end
       
       function build_behavior(obj)
            % Process behavior pulls data to build raw waveform matrix
            v = obj.verbose;
            if v>1; disp([newline '--> Building behavioral data...']); end
            
            % Find cycle onset times
            cycle_times = obj.data.cycleTimes.data;

            % Make a cell array to hold cycle data
            nCycles = size(cycle_times,1);
            behaviorCycles = cell(nCycles,1);
            
            
            behavOffset = 0;
            % Iterate and load data file info
%             for i=1:length(obj.arrFiles)
            for i=1:3
                if obj.verbose > 1
                    disp('===== ===== ===== ===== =====');
                    disp(['Processing file ' num2str(i) ' of ' num2str(length(obj.arrFiles))]);
                    disp(['File: ' obj.strFldr '\' obj.arrFiles{i}]);
                end
                [pressure_ts, pressure_wav] = read_Intan_RHD2000_nongui_adc(fullfile(obj.strFldr, obj.arrFiles{i}), obj.verbose);

                
                if v>2; disp(['Start time: ' num2str(pressure_ts(1)*1000.) ' ms']); end
                
                % Filter Pressure Waves
%                 filterData = obj.filterBehavior(pressure_wav, obj.Fs, filterFreq); % This will change once we update filterBehavior func
                if v>3; disp('--> --> Filtering data...');
                filterData = pressure_wav;

                
                % Find cycles that occur within the limits of the current
                % data file
                cycleIxs = cycle_times(:,1) >= pressure_ts(1)*1000. & cycle_times(:,2) <= pressure_ts(end)*1000.;
                validCycles = cycle_times(cycleIxs,:);
                tStart = pressure_ts(1)*1000.; % beginning of data file in ms
                
                if v>3; disp(['# Cycles: ' num2str(sum(cycleIxs))]); end
                
                % Assign pressure waves to matrix rows
                % Consider alternative ways to save speed here?
                for iCycle = 1:size(validCycles,1)
                   cycleStart = ceil((validCycles(iCycle,1)-tStart)*obj.Fs/1000.);
                   cycleEnd = floor((validCycles(iCycle,2)-tStart)*obj.Fs/1000.);
                   behaviorCycles{iCycle+behavOffset} = filterData(cycleStart:cycleEnd);
                end
                
                behavOffset = behavOffset + size(validCycles,1);
            end
            
            obj.rawBehav = behaviorCycles;
            if v>0; disp('COMPLETE: Behavioral data loaded!'); end
       end
       
       
       function [filterData] = filterBehavior(obj, behavior, cycleFreq, filterFreq)
           
            if v>1; disp([newline '--> Filtering behavioral data...']); end 
           
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

            if v>0; disp('COMPLETE: Behavioral data filetered!'); end
       end        
        
       function r = get_behavior(obj, timeBase, feature, start, dur, nSamp, varargin)
            %% Implemented to return different features of behavior cycles after processing raw waveform matrix data
            % i.e., raw data, PCA, residual, area

            p = inputParser;
            
            % Required: timeBase
            % 'time' or 'phase' to indicate how to calculate time
            valid_timeBase = {'time' 'phase'};
            validate_timeBase = @(x) assert(ischar(x) && ismember(x, valid_timeBase), 'timeBase must be a char/string: time, phase');
            p.addRequired('timeBase', validate_timeBase);
            
            % Required: feature
            % 'raw', 'pca', 'residual'
            valid_feature = {'raw' 'pca' 'residual'};
            validate_feature = @(x) assert(ischar(x) && ismember(x, valid_feature), 'feature must be a char/string: raw, pca, residual');
            p.addRequired('feature', validate_feature);
            
            % Required: start
            % in ms or rad/phase based on timeBase
            validate_start = @(x) assert(isnumeric(x), 'start must be numeric corresponding to timeBase');
            p.addRequired('start', validate_start);
            
            % Required: dur
            % in ms or rad/phase based on timeBase
            validate_dur = @(x) assert(isnumeric(x), 'dur must be numeric indicating the duration of the vector corresponding to timeBase');
            p.addRequired('dur', validate_dur);
            
            % Required: nSamp
            % in samples/data points
            validate_nSamp = @(x) assert(isnumeric(x) && mod(x,1) == 0, 'nSamp must be an integer indicating number of data points in sample');
            p.addRequired('nSamp', validate_nSamp);
            
            % Optional: nPC
            % integer to indicate number of principal components to return
            default_nPC = 3;
            validate_nPC = @(x) assert(isnumeric(x) && mod(x,1) == 0, 'nPC must be an integer indicating number of principal components to return');
            p.addOptional('nPC', default_nPC, validate_nPC);
            
            p.parse(timeBase, feature, start, dur, nSamp, varargin{:});
            
            timeBase = p.Results.timeBase;
            feature = p.Results.feature;
            start = p.Results.start;
            dur = p.Results.dur;
            nSamp = p.Results.dur;
            nPC = p.Results.nPC;
            
            
            if v>1; disp([newline '--> Getting formatted behavioral data...']); end
            
            switch(timeBase)
                case('time')
                    if v>2; disp([newline '--> --> Formatting behavior: time']); end
                    % Find the lengths of the cycles in samples
                    cycleLengths_samples = cell2mat(cellfun(@(x) length(x), obj.rawBehav, 'UniformOutput', false));

                    % Find the sample associated with the start time for each
                    % cycle in ms
                    start_samples = ceil(start*obj.Fs/1000.);

                    % Find the number of samples that encompases the window of
                    % interest for the cycles. 
                    % Convert windowOFInterest from ms to seconds
                    dur_samples = ceil(dur*obj.Fs/1000.);

                    % Find the stop sample of the window of interest for each
                    % cycle
                    stop_samples = start_samples + dur_samples;

                    if v>3
                        disp(['Cycle Start Time: ' num2str(start) ' ms']);
                        disp(['Cycle Sample Duration: ' num2str(dur) ' ms']);
                        disp(['Sampled Data Points: ' num2str(nSamp)]);
                    end
                    
                    % Set up empty matrix to store pressure data.
                    nCycles = length(cycleLengths_samples);
                    cycle_behavior = zeros(nCycles, nSamp);

                    tic
                    for cycle_ix = 1:nCycles
                       % Document all of the data points for the window of
                       % interest
                       dat = obj.rawBehav{cycle_ix};
                       if length(dat) >= stop_samples
                           cycle_data = dat(start_samples:stop_samples);
                           % Resample to get only the desired number of points
                           newSamples = round(linspace(1,length(cycle_data),nSamp));

                           resampled_cycle_data = cycle_data(newSamples);
                           cycle_behavior(cycle_ix,:) = resampled_cycle_data;   
                       else
                           cycle_behavior(cycle_ix,:) = nan(1,nSamp);
                       end
                    end
                    toc
                    
                    if v>2; warning(['WARNING: Ommitted ' sum(isnan(cycle_behavior(:,1))) ' empty cycles']); end

                case('phase')
                    if v>2; disp([newline '--> --> Formatting behavior: phase']); end
                    % Find the lengths of the cycles in samples
                    cycleLengths_samples = cell2mat(cellfun(@(x) length(x), obj.rawBehav, 'UniformOutput', false));

                    % Find the sample associated with the start phase for each
                    % cycle
                    start_samples = ceil(cycleLengths_samples.*(start/(2*pi)));

                    % Find the number of samples that encompases the window of
                    % interest for the cycles.
                    windowOfInterest_samples = ceil((cycleLengths_samples).*(dur./(2*pi)));

                    % Find the stop sample of the window of interest for each
                    % cycle
                    stop_samples = start_samples + windowOfInterest_samples;

                    if v>3
                        disp(['Cycle Start Phase: ' num2str(start) ' rad']);
                        disp(['Cycle Sample Duration: ' num2str(dur) ' rad']);
                        disp(['Sampled Data Points: ' num2str(nSamp)]);
                    end
                    
                    % Set up empty matrix to store pressure data.
                    nCycles = length(cycleLengths_samples);
                    cycle_behavior = zeros(nCycles, nSamp);

                    tic
                    for cycle_ix = 1:nCycles
                        % Document all of the data points for the window of
                        % interest
                        if stop_samples(cycle_ix) > start_samples(cycle_ix)
                            dat = obj.rawBehav{cycle_ix};
                            cycle_data = dat(start_samples(cycle_ix):stop_samples(cycle_ix));

                            newSamples = round(linspace(1,length(cycle_data),nSamp));
                            % Resample to get only the desired number of points
                            resampled_cycle_data = cycle_data(newSamples);
                            cycle_behavior(cycle_ix, 1:nSamp) = resampled_cycle_data;   
                        else
                            cycle_behavior(cycle_ix,:) = nan(1,nSamp);
                        
                        end
                    end              
                    toc
                    
                    if v>2; warning(['WARNING: Ommitted ' sum(isnan(cycle_behavior(:,1))) ' empty cycles']); end
            end
            
            switch(feature)
                case('raw')
                    if v>2; disp([newline '--> --> Processing behavior: raw']); end
                    r = cycle_behavior;
                case('pca')
                    if v>2; disp([newline '--> --> Processing behavior: PCA (' num2str(nPC) ') PCs']); end
                    [~,score,~] = pca(cycle_behavior);
                    r = score(:,1:nPC);
                case('resid')
                    if v>2; disp([newline '--> --> Processing behavior: residual']); end
                    r = cycle_behavior - mean(cycle_behavior,1);     
            end
            
            if v>0; disp('COMPLETE: Behavioral data retrieved!'); end
        end       
   end    
end

%            obj.b_startPhase = {};