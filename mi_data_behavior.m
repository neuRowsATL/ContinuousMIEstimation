classdef mi_data_behavior < handle
    properties
        raw_source % array or struct of raw files and data information
        
        cycleMethod % method used to calculate cycle times
        cycleTimes % cycle onset times
        cycleData % MAT object of raw data for each cycle
        
        verbose % boolean flag used to track progress and errors
    end
    
    methods
        function obj = mi_data_behavior(varargin)
            p = inputParser;
            addParameter(p, 'verbose', 1);
            parse(p, varargin{:});
            obj.verbose = p.Results.verbose;
        end
        
        function set_source()
            %% This function is implemented to set the raw source variable based on data collection and storage methods
            warning('Not implemented error: set_source()');
        end
        
        function calc_cycles()
            %% Implemented to take raw behavior recording and identify onset times
            warning('Not implemented error: calc_cycles()');
        end
        
        function get_feature()
            %% Implemented to return different features of behavior cycles
            % i.e., raw data, PCA, residual, area
            warning('Not implemented error: get_feature()');
        end
        
        %         function r = processBehavior(obj)
%             
%             behavData = obj.obj_behav.cycleData;
%             
% %             if obj.verbose > 2
% %                 figure()
% %                 h_avgBehav = plot(nanmean(behavData), 'k-', 'LineWidth', 3);
% %                 hold on;
% %             end
%             
%             
%             switch(obj.b_dataTransform)
%             case('pca')
%                 switch(obj.b_timebase)
%                 case('phase')
%                     % Find the lengths of the cycles in samples
% %                     cycleLengths_samples = sum(~isnan(obj.behavior),2);
%                     cycleLengths_samples = sum(~isnan(behavData),2);
% 
%                     % Find the sample associated with the start phase for each
%                     % cycle
%                     start_samples = ceil(obj.b_startPhase.*(cycleLengths_samples./(2*pi)));
% 
%                     % Find the number of samples that encompases the window of
%                     % interest for the cycles. 
%                     windowOfInterest_samples = ceil(obj.b_windowOfInterest.*(cycleLengths_samples./(2*pi)));
% 
%                     % Find the stop sample of the window of interest for each
%                     % cycle
%                     stop_samples = start_samples + windowOfInterest_samples;
% 
%                    % Set up empty matrix to store pressure data.
%                    nCycles = length(cycleLengths_samples);
%                    cycle_behavior = nan(nCycles, 1000);
% 
%                    for cycle_ix = 1:nCycles
%                        % Document all of the data points for the window of
%                        % interest
% %                        cycle_data = obj.behavior(cycle_ix,start_samples(cycle_ix):stop_samples(cycle_ix));
%                        cycle_data = behavData(cycle_ix,start_samples(cycle_ix):stop_samples(cycle_ix));
%                        % We resample the cyclic data to make the data have
%                        % uniform length --> I am arbitrarily picking 1000
%                        % data points. 
%                        resampled_cycle_data = resample(cycle_data,1000,length(cycle_data));
%                        cycle_behavior(cycle_ix, 1:1000) = resampled_cycle_data; 
%                    end
% 
%                 case('time')
%                 % Find the lengths of the cycles in samples
% %                     cycleLengths_samples = sum(~isnan(obj.behavior),2);
%                     cycleLengths_samples = sum(~isnan(behavData),2);
% 
%                     % Find the sample associated with the start time for each
%                     % cycle
%                     % Convert startPhase from ms to seconds
%                     startTime = obj.b_startPhase/1000;
%                     start_samples = ceil(startTime*obj.bFs);
% 
%                     % Find the number of samples that encompases the window of
%                     % interest for the cycles. 
%                     % Convert windowOFInterest from ms to seconds
%                     windowOfInterest_seconds = obj.b_windowOfInterest./1000;
%                     windowOfInterest_samples = ceil(windowOfInterest_seconds*obj.bFs);
% 
%                     % Find the stop sample of the window of interest for each
%                     % cycle
%                     stop_samples = start_samples + windowOfInterest_samples;
% 
%                    % Set up empty matrix to store pressure data.
%                    nCycles = length(cycleLengths_samples);
%                    cycle_behavior = nan(nCycles, windowOfInterest_samples);
% 
%                    for cycle_ix = 1:nCycles
%                        % Document all of the data points for the window of
%                        % interest
% %                        cycle_data = obj.behavior(cycle_ix,start_samples:stop_samples); 
%                        cycle_data = behavData(cycle_ix,start_samples:stop_samples); 
%                        cycle_behavior(cycle_ix, 1:windowOfInterest_samples+1) = cycle_data;
%                    end
%                 end
%                 % Find the PCs of the window of interest cycle data
%                 [~,score,~] = pca(cycle_behavior);
%                 cycle_dataPCs = score(:,1:obj.b_length);
%                 r = cycle_dataPCs;
%                 
%                 otherwise % 'none' or 'residual'
%                     switch(obj.b_timebase)
%                         case('phase')
%                             % Find the lengths of the cycles in samples
%         %                     cycleLengths_samples = sum(~isnan(obj.behavior),2);
%                             cycleLengths_samples = sum(~isnan(behavData),2);
% 
%                             % Find the sample associated with the start phase for each
%                             % cycle
%                             start_samples = ceil(obj.b_startPhase.*(cycleLengths_samples./(2*pi)));
% 
%                             % Find the number of samples that encompases the window of
%                             % interest for the cycles. 
%                             windowOfInterest_samples = ceil(obj.b_windowOfInterest.*(cycleLengths_samples./(2*pi)));
% 
%                             % Find the stop sample of the window of interest for each
%                             % cycle
%                             stop_samples = start_samples + windowOfInterest_samples;
% 
%                            % Set up empty matrix to store pressure data.
%                            nCycles = length(cycleLengths_samples);
%                            cycle_behavior = nan(nCycles, obj.b_length);
%                            
%                            for cycle_ix = 1:nCycles
%                                % Document all of the data points for the window of
%                                % interest
%         %                        cycle_data = obj.behavior(cycle_ix,start_samples(cycle_ix):stop_samples(cycle_ix));
%                                cycle_data = behavData(cycle_ix,start_samples(cycle_ix):stop_samples(cycle_ix));
% 
% %                                if obj.verbose > 2; plot(start_samples(cycle_ix):stop_samples(cycle_ix), cycle_data); end
%                                
%                                % Resample to get only the desired number of points
%                                resampled_cycle_data = resample(cycle_data,obj.b_length,length(cycle_data));
%                                cycle_behavior(cycle_ix, 1:obj.b_length) = resampled_cycle_data;   
%                            end
% 
%                         case('time')
%                         % Find the lengths of the cycles in samples
%         %                     cycleLengths_samples = sum(~isnan(obj.behavior),2);
%                             cycleLengths_samples = sum(~isnan(behavData),2);
% 
%                             % Find the sample associated with the start time for each
%                             % cycle
%                             % Convert startPhase from ms to seconds
%                             startTime = obj.b_startPhase/1000;
%                             start_samples = ceil(startTime*obj.bFs);
% 
%                             % Find the number of samples that encompases the window of
%                             % interest for the cycles. 
%                             % Convert windowOFInterest from ms to seconds
%                             windowOfInterest_seconds = obj.b_windowOfInterest/1000;
%                             windowOfInterest_samples = ceil(windowOfInterest_seconds*obj.bFs);
% 
%                             % Find the stop sample of the window of interest for each
%                             % cycle
%                             stop_samples = start_samples + windowOfInterest_samples;
% 
%                            % Set up empty matrix to store pressure data.
%                            nCycles = length(cycleLengths_samples);
%                            cycle_behavior = nan(nCycles, obj.b_length);
% 
%                            figure();
%                            
%                            for cycle_ix = 1:nCycles
%                                % Document all of the data points for the window of
%                                % interest
%         %                        cycle_data = obj.behavior(cycle_ix,start_samples:stop_samples);
%                                cycle_data = behavData(cycle_ix,start_samples:stop_samples);
%                                % Resample to get only the desired number of points
%                                
%                                plot(behavData(cycle_ix,:));
%                                hold on;
%                                plot(start_samples:stop_samples, cycle_data);
%                                
%                                resampled_cycle_data = resample(cycle_data,obj.b_length,length(cycle_data));
%                                cycle_behavior(cycle_ix, 1:obj.b_length) = resampled_cycle_data;   
%                            end
%                     end
% 
%                     % Transform the behavioral data
%                     switch(obj.b_dataTransform)
%                         case('none')
%                             r = cycle_behavior;
%                         case('residual')
%                             avg_cycle = mean(cycle_behavior,1);
%                             avg_cycleMatrix = repmat(avg_cycle, length(cycle_behavior),1);
%                             cycle_residuals = cycle_behavior - avg_cycleMatrix;
%                             r = cycle_residuals;
% 
%                     end
%             end
% 
% %             if obj.verbose > 2; uistack(h_avgBehav, 'top'); end
% 
%         end
        
    end
end