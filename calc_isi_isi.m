classdef calc_isi_isi < mi_analysis
    %Each of these objects sets the stage to calculate the mutual
    %information between spike count and behavior and stores the results of
    %the calculation. 
    %   Detailed explanation goes here
    
    properties
        isi_cutoff % ms
        isi_offset % number of ISIs to offset
        noise % ms; std of noise to add to isi
        shuffle % T/F to run code on mixed up spike times
    end
    
    methods
% NEED TO CONSIDER: where do we want to store ISI specific properties to
% modify how the data is processed.  we can pull get spikes in raw format
% here and then process the data.

        function obj = calc_isi_isi(objData, objBehav, varNames, varargin) 
            % isi_offset, isi_cutoff, noise
            % Required arguments: objData, varNames
            % Check required inputs for validity using input parser
            
            % Set up input parser
            p = inputParser;
            
            % Set required inputs
            validate_objData = @(x) assert(isa(x, 'mi_data_neural'), 'objData must be a neural data subclass');
            p.addRequired('objData', validate_objData);
            
            validate_objBehav = @(x) assert(isa(x, 'mi_data_behavior'), 'objBehav must be a behavioral data subclass');
            p.addRequired('objBehav', validate_objBehav);
            
            validate_varNames = @(x) assert(iscell(x) && (length(x) == 1), 'varNames must be a cell of length 1');
            p.addRequired('varNames', validate_varNames);
            
            % Set parameters
            default_isi_cutoff = 200;
            validate_isi_cutoff = @(x) assert(isinteger(x), 'isi_cutoff must be an integer in ms units');
            p.addParameter('isi_cutoff', default_isi_cutoff, validate_isi_cutoff);
            
            default_isi_offset = 1;
            validate_isi_offset = @(x) assert(isinteger(x), 'isi_offset must be an integer');
            p.addParameter('isi_offset', default_isi_offset, validate_isi_offset);
            
            default_noise = 0;
            validate_noise = @(x) assert(isnumeric(x), 'noise must be a numeric quantity in ms');
            p.addParameter('noise', default_noise, validate_noise);
            
            default_shuffle = false;
            validate_shuffle = @(x) assert(isboolean(x), 'shuffle must be a boolean value');
            p.addParameter('shuffle', default_shuffle, validate_shuffle);
            
            
            % Prepare InputParser to parse only desired inputs
            p.KeepUnmatched = 1;
            p.parse(objData, objBehav, varNames, varargin{:});
            
            % Define validated inputs to parent constructor
            objData = p.Results.objData;
            objBehav = p.Results.objBehav;
            varNames = p.Results.varNames;

            
            % One more validation: Check that varNames references valid fields of objData
            for ivarNames = 1:length(varNames)
                assert(isfield(objData.data , varNames{ivarNames}), ['varName: ' varNames{ivarNames} 'is not a valid field of the neural data object']); 
            end
            
            % Call parent constructor
            obj@mi_analysis(objData, objBehav, varNames, varargin{:});
            
            % Define other properties of subclass object
            obj.isi_cutoff = p.Results.isi_cutoff;
            obj.isi_offset = p.Results.isi_offset;
            obj.noise = p.Results.noise;
            obj.shuffle = p.Results.shuffle;
            
        end
        
        function buildMIs(obj)
            % So I propose that we use this method to prep the
            % count_behavior data for the MI core and go ahead and run MI
            % core from here. Then we can use the output of MI core to fill
            % in the MI, kvalue, and errors.
            
            % First, get spike times from neuron
            spikeTimes = obj.objData.get_spikes('format', 'raw', 'name', obj.varNames{1,1});
            
            % BC-20190129: DO WE NEED TO REWRITE THIS CODE FOR INCREASED FLEXIBILITY...
            % WE CAN REWRITE IT SO THAT IT CAN TAKE THE MI BETWEEN ANY TWO SERIES OF ISI... ?
            % For example:
            % - interspike intervals within same spike train but with varying time delays
            % --> ISI_n | ISI_n+1
            % --> ISI_n | ISI_n+2
            % --> ISI_n | ISI_n+3
            
            % Find ISIs from spike times
            ISIs = diff(spikeTimes);
            if obj.shuffle
                ISIs = ISIs(randperm(length(ISIs)));
            end
            
            % I assume we want to keep all ISIs that we are comparing to
            % within a burst/cycle...
            consec_isi = [ISIs(1:end-obj.isi_offset) ISIs(obj.isi_offset + 1:end)];
            consec_sum = sum(consec_isi, 2);
            select_isi = consec_sum <= obj.isi_cutoff;
            consec_isi = consec_isi(select_isi,1:end);
            
            % should this be normrnd(0,noise)?
            jitter = obj.noise*randn(size(consec_isi));
            mi_isis = consec_isi + jitter;

            % Make a vector of the first ISIs
            x = mi_isis(:,1);
            
            avg_x = mean(x);
            std_x = std(x);
            
            % RC 20200407: This should be covered by the reparam flag in
            % the parent class now.
%            xGroups{1,1} = (x-avg_x)/std_x;
             xGroups{1,1} = x;
            

            % Make a vector of the second ISIs
            y = mi_isis(:,2);
            
            avg_y = mean(y);
            std_y = std(y);
            
            % RC 20200407: This should be covered by the reparam flag in
            % the parent class now. 
%            yGroups{1,1} = (y-avg_y)/std_y;
             yGroups{1,1} = y;
            
            coeffs = {1};
            
            buildMIs@mi_analysis(obj, {xGroups yGroups coeffs});
        end
    end
end

