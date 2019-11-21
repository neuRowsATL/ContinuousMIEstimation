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
           
           p = inputParser;
           
           % Required input: arrDataFiles
           validate_arrDataFiles = @(x) assert(iscell(x) && ischar(x{1}), 'arrDataFiles needs to be cell array of string file names');
           p.addRequired('arrDataFiles', validate_arrDataFiles);
           
           % Optional input: strFldrPath
           default_strFldrPath = pwd;
           validate_strFldrPath = @(x) assert(true, 'strFldrPath needs to be a string path');
           p.addOptional('strFldrPath', default_strFldrPath, validate_strFldrPath);
           
           p.parse(arrDataFiles, varargin{:});

           if isfolder(p.Results.strFldrPath)
               filepaths = fullfile(p.Results.strFldrPath, p.Results.arrDataFiles);
               for i=1:length(filepaths)
                   if ~isfile(filepaths(i)); error(['File does not exist: ' filepaths(i)]); end
               end
                obj.arrFiles = p.Results.arrDataFiles;
                obj.strFldr = p.Results.strFldrPath;
           else
               error(['strFldrPath does not exist: ' p.Results.strFldrPath]);
           end
       end
       
       function build_behavior(obj)
            % Process behavior pulls data to build raw waveform matrix
           
            % Find cycle onset times
            cycle_times = obj.data.cycleTimes.data;

            % Make an NaN cell array to hold cycle data. 
            nCycles = size(cycle_times,1);
            behaviorCycles = cell(nCycles,1);
            
            
            behavOffset = 0;
            % Iterate and load data file info
%             for i=1:length(obj.arrFiles)
            for i=1:3
                if obj.verbose > 2
                    disp('===== ===== ===== ===== =====');
                    disp(['Processing file ' num2str(i) ' of ' num2str(length(obj.arrFiles))]);
                    disp(['File: ' obj.strFldr '\' obj.arrFiles{i}]);
                end
                [pressure_ts, pressure_wav] = read_Intan_RHD2000_nongui_adc(fullfile(obj.strFldr, obj.arrFiles{i}), obj.verbose);

                
                % Filter Pressure Waves
%                 filterData = obj.filterBehavior(pressure_wav, obj.Fs, filterFreq); % This will change once we update filterBehavior func
                filterData = pressure_wav;

                
                % Find cycles that occur within the limits of the current
                % data file
                cycleIxs = cycle_times(:,1) >= pressure_ts(1)*1000. & cycle_times(:,2) <= pressure_ts(end)*1000.;
                validCycles = cycle_times(cycleIxs,:);
                tStart = pressure_ts(1)*1000.; % beginning of data file in ms
                
                
                % Assign pressure waves to matrix rows
                % Consider alternative ways to save speed here
                for iCycle = 1:size(validCycles,1)
                   cycleStart = ceil((validCycles(iCycle,1)-tStart)*obj.Fs/1000.);
                   cycleEnd = floor((validCycles(iCycle,2)-tStart)*obj.Fs/1000.);
                   behaviorCycles{iCycle+behavOffset} = filterData(cycleStart:cycleEnd);
                end
                
                behavOffset = behavOffset + size(validCycles,1);
            end
            
            obj.rawBehav = behaviorCycles;
       end
       
       
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
        
       function r = get_behavior(obj, timeBase, feature, start, dur, nSamp)
%         function r = get_behavior(obj, format, feature, start, dur, nSamp, nPC)
            %% Implemented to return different features of behavior cycles after processing raw waveform matrix data
            % i.e., raw data, PCA, residual, area

            switch(timeBase)
                case('time')
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

                case('phase')
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
            end
            
            switch(feature)
                case('raw')
                    r = cycle_behavior;
                case('pca')
                    [~,score,~] = pca(cycle_behavior);
                    r = score(:,1:nPC);
                case('resid')
                    r = cycle_behavior - mean(cycle_behavior,1);     
            end
        end       
   end    
end

%            obj.b_startPhase = {};