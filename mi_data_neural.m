classdef mi_data_neural < mi_data
    
    methods
        function add_spikes()
            
        end
        
        %        function r = getTiming(obj, dataNum)
%            % Makes matrices of pressure data  need to decide what units to use and spiking data based on what we have      
%            % NOTE- currently as the code is written, we omit any neural or pressure data that occurs	 
%            % before the onset of the first cycle or after the onset of the last cycle
%            %  - additionally, we are segmenting spikes based on the cycle times rather than trying
%            % to keep bursts together and using negative spike times 
%            %
%            % INPUT
%            % dataNum : positive integer neuron number to specify neuron of interest
%            %
%              % Return an m x n matrix of data with m cycles and n sample points (if pressure) 
%              % or n maximum spikes (if neuron data) 
%            %
%            % if verbose > 1; disp([newline 'Running: dataByCycles' newline]); end
% 
%            spike_ts = obj.neurons{dataNum{1}};
% %            cycle_ts = obj.cycleTimes{1,1};
%            cycle_ts = obj.obj_behav.cycleTimes'*1000.; % BC20190820
% 
%            % Find the number of spikes in each cycle
%            cycle_spike_counts = obj.getCount(dataNum);
% 
%            % Calculate relative spike times for each breathing cycle
%            % if verbose > 1; disp('-> Calculating relative spike times by cycle'); end
%            cycle_spike_ts = nan(size(cycle_ts,1)-1, max(cycle_spike_counts));
%            for cycle_ix = 1:(size(cycle_ts,1)-1)
%                cycle_spikes_ix = find((spike_ts > cycle_ts(cycle_ix)) & (spike_ts < cycle_ts(cycle_ix+1)));
%                if ~isempty(cycle_spikes_ix)
%                    cycle_spike_ts(cycle_ix,1:length(cycle_spikes_ix)) = spike_ts(cycle_spikes_ix)-cycle_ts(cycle_ix);
%                    disp('goodbye');
%                end
%                disp('goodbye2');
%            end
%            switch(obj.n_timebase)
%                case('phase')
%                    % Convert spike times to phase values in radians. 
%                    cycle_lengths = diff(cycle_ts);
%                    % Find the dimensions of the cycle_spike_ts matrix
%                    dimension_cycle_spike_ts = size(cycle_spike_ts);
%                    % Calculate the phase conversion for each cycle
%                    phase_factor = (2*pi)./cycle_lengths;
%                    % Propogate the phase conversion to a matrix
%                    phase_factor_matrix = repmat(phase_factor,1,dimension_cycle_spike_ts(2));
%                    % Multiply each spike time by the phase factor for the
%                    % respective cycle
%                    cycle_spike_phase = cycle_spike_ts.*phase_factor_matrix;
%                    % Output the variable. 
%                    r = cycle_spike_phase;
%                    
%                case('time')
%                    r = cycle_spike_ts;
%            end
%            
%            disp('hello');
%        
%        end
%        
%        function r = getCount(obj, dataNum)
%            spike_ts = obj.neurons{dataNum{1}};
% %            cycle_ts = obj.cycleTimes{1,1};
%             cycle_ts = obj.obj_behav.cycleTimes'*1000.; % BC 20190820
% 
%            % Find the number of spikes in each cycle
%            % We include data that comes after the onset of the first cycle
%            % and before the onset of the last cycle
%            cycle_spike_counts = zeros(1,size(cycle_ts,1)-1);
%            for cycle_ix = 1:(size(cycle_ts,1)-1)
%                cycle_spikes_ix = find((spike_ts > cycle_ts(cycle_ix)) & (spike_ts < cycle_ts(cycle_ix+1)));
%                if ~isempty(cycle_spikes_ix)
%                  cycle_spike_counts(cycle_ix) = length(cycle_spikes_ix);
%                end
%            end
%            r = cycle_spike_counts;
%        end
end