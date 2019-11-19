classdef mi_data_behavior < mi_data
    properties
        rawBehav % struct of raw behavior data and information
        
        arrFiles % list of files with raw behavioral data
    end
    
    methods
        function obj = mi_data_behavior(ID, varargin)
            % Required arguments: ID
            obj@mi_data(ID,varargin{:});
        end
        
        function add_cycleTimes(obj, data, dataInfo, Fs, varargin)
            add_data(obj, data, dataInfo, Fs, 'cycleTimes', varargin{:});
        end
        
        function build_behav(obj)
            % Pull waveform data from data files according to
            % data/cycleTimes
            warning('NOT IMPLEMENTED ERROR');
        end
        
        function r = get_behav(obj)
            % Return behavior data as matrix of raw waveform/PCA/ICA/etc.
            warning('NOT IMPLEMENTED ERROR');
        end
        
        function r = get_cycleTimes(obj)
            % Return behavior cycle times as matrix of on/off pairs
            r = get_data(obj, 'cycleTimes');
        end
    end
end