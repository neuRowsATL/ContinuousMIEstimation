classdef mi_ksg_viz < handle
    properties
        
    end
    methods
        function plot_generic()
            % make generic plot of two variables
        end
         function r_plot = plot_data_fraction(obj, obj_core, k, f, ax)
            % make data fraction plot
            if nargin == 3
                fig = figure();
            elseif nargin > 3
                fig = figure(f);
            end
            
            
            bool_ixs = cell2mat(obj_core.mi_data(:,4)) == k;
            xs = cell2mat(obj_core.mi_data(bool_ixs,3));
            ys = cell2mat(obj_core.mi_data(bool_ixs,1));
            err = cell2mat(obj_core.mi_data(bool_ixs,2));
            r_plot = errorbar(ax, xs, ys, err, '-b', 'Marker', '.', 'MarkerSize', 15);
            
            %xlabel('Data Fraction (1/N)');
            %ylabel('Mutual Information');
            title({[ 'k = ' num2str(k)]});
            
            xlim([min(xs)*0.8 max(xs)*1.1]);
            
        end
         function r_plot = plot_k_dependence(obj, obj_core, f, ax)
            % make k-dependence plot
            if nargin == 2
                fig = figure();
            elseif nargin > 2
                fig = figure(f);
            end
            
            
            ks = obj_core.k_values;
            ys = zeros(1, length(ks));
            err = zeros(1, length(ks));
            for k_ix=1:length(ks)
                dat = get_mi(obj_core, ks(k_ix));
                ys(k_ix) = dat.mi;
                err(k_ix) = dat.err;
            end
            r_plot = errorbar(ks, ys, err, '-b', 'Marker', '.', 'Markersize', 15);
            
            %xlabel('k-value');
            %ylabel('Mutual Information');
            title({'Kraskov-Stoegbauer-Grassberger' 'k-dependence'});
            
            xlim([min(ks)*0.8 max(ks)*1.1]);
        end
        function plot_mi_estimates()
            % make plot of mutual information vs. parameter
        end
        function r_plot = make_ksg_graph(obj, obj_core, f)
            % make plot of KSG data points for k-NN estimator visualization
            if nargin == 2
                fig = figure();
            elseif nargin > 2
                fig = figure(f);
            end
            
            ax = subplot(1,1,1);
            r_plot = scatter(ax, obj_core.x, obj_core.y, 15, 'b', 'filled');
            xlabel('x');
            ylabel('y');
            title({'Kraskov-Stoegbauer-Grassberger' 'Nearest-Neighbors Plot'});
        end
        function make_mi_scale_graph()
            % make plot of values with scaled color in z-plane
        end
        function plot_error_estimate()
            % make linear regression plot to estimate error
        end
    end
end
