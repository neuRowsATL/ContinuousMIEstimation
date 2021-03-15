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
            var = cell2mat(obj_core.mi_data(bool_ixs,2));
            err = var.^.5;
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
                dat = get_mi(obj_core, -1,'k', ks(k_ix));
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
    
    methods(Static)
        
        % Audit Plots from mi_analysis class
        function audit_plots(mi_analysis)
            
	    % Iterating through the core objects
            for iGroup = 1:size(mi_analysis.arrMIcore,1)
                coreObj = mi_analysis.arrMIcore{iGroup,1};
                
                % FOR NOW, NO AUDIT PLOTS FOR BEHAVIOR SUBCLASSES
                if contains(class(mi_analysis), 'behav')
                    continue
                else
                    
                    % Check for data type
                    % Histograms do not depend on data type
                    % First make histogram with x data
                    % RC 20191213: We should come back and set specific bin widths here.
                    % The only issue we may run into is
                    x = coreObj.x;
                    figure()
                    histogram(x)
                    hold on
                    xlabel('X Value (binned)')
                    ylabel('N Cycles')
                    title('Histogram for X')
                    
                    % Histogram for y
                    y = coreObj.y;
                    figure()
                    histogram(y)
                    hold on
                    xlabel('Y Value (binned)')
                    ylabel('N Cycles')
                    title('Histogram for Y')
                    
                    % For multidimensional data
                    if all(size(x) > 1) & all(size(y) > 1)
                        % Plot the first pc1x and pc1y against each other
                        % Variability plot for x and y individually 
                        
                        % Obtain relevant PCA values
                        [~, scorex, latentx] = pca(x);
                        [~, scorey, latenty] = pca(y); 
                        
                        % Plot the first components against each other 
                        figure()
                        scatter(scorex(:,1), scorey(:,1), 'x');
                        axis equal; 
                        xlabel('1st X Principle Component')
                        ylabel('1st Y Principle Component')
                        title('PCA Joint Plot')
                        
                        % Plot variability in x components 
                        figure()
                        cumsumx = sum(latentx);
                        perVarx = latentx / cumsumx;
                        bar(perVarx)
                        xlabel('Principle Component')
                        ylabel('Percent Variability Explained')
                        title('Scree Plot for X')
                        
                        % Plot variability in y components 
                        figure()
                        cumsumy = sum(latenty);
                        perVary = latenty / cumsumy;
                        bar(perVary)
                        xlabel('Principle Component')
                        ylabel('Percent Variability Explained')
                        title('Scree Plot for Y')
                    else
                        % Check for discrete data in both variables
                        if all(all(rem(x,1) == 0)) & all(all(rem(y,1) == 0))
                            try
                                % For discrete data, plot a jittered histogram

                                % Add noise for joint histogram
                                x_plot = x + 0.2*rand(size(x));
                                y_plot = y + 0.2*rand(size(y));

                                % Make density map 
                                pts = linspace(0, max(max(x), max(y)) + 1, max(max(x), max(y)) + 2); 
                                N = histcounts2(y(:), x(:), pts, pts);
                                N = log(N);

                                % Plot scattered data:
                                figure()
                                scatter(x_plot, y_plot, 'x');
                                axis equal;
                                set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]));
                                xlabel('Discrete Value: X')
                                ylabel('Discrete Value: Y')
                                title('P(X,Y) Discrete Joint Distribution')

                                % Plot heatmap:
                                figure()
                                imagesc(pts, pts, N);
                                axis equal;
                                set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
                                xlabel('Discrete Value: X')
                                ylabel('Discrete Value: Y')
                                title('P(X,Y) Density Map')
                                cBar = colorbar('Direction', 'normal', 'Ticks', 1:max(N, [], 'all'));
                                cBar.Label.String = "log(# of data points)";
                            catch
                                continue
                            end
                            
                            
                        elseif all(all(rem(x,1) == 0)) | all(all(rem(y,1) == 0))
                            % Try to catch any empty vectors. 
                            try
                                % Add jitter only to the variable that is discrete, which for our data, will always be the second variable.
                                if all(rem(x,1) == 0)
                                    x_plot = x + 0.2*rand(size(x));
                                    x_L = 'Discrete Value: X';
                                else
                                    x_plot = x;
                                    x_L = 'Continuous Value: X';
                                end
                                if all(rem(y,1) == 0)
                                    y_plot = y + 0.2*rand(size(y));
                                    y_L = 'Discrete Value: Y';
                                else
                                    y_plot = y;
                                    y_L = 'Continuous Value: Y';
                                end

                                % Make figure
                                figure()
                                plot(x_plot, y_plot, 'x')
                                hold on
                                xlabel(x_L)
                                ylabel(y_L)
                                title('P(X,Y) Mixed Joint Distribution')
                            catch
                                continue
                            end
                        else
                            % The assumption is that both distributions are continuous if neither of the above if statements are true.
                            try
                                % Make figure
                                figure()
                                plot(x,y, 'x')
                                hold on
                                xlabel('Continuous Variable: X')
                                ylabel('Continuous Variable: Y')
                                title('P(X,Y) Continuous Joint Distribution')
                            
                            catch
                                continue
                            end
                        end
                        
                    end
                end
                
            end
        end
        
        function rasterPlots(mi_analysis)
            for iNeuron = 1:length(mi_analysis.varNames)
                spikeTimes = mi_analysis.objData.get_spikes('name', mi_analysis.varNames{1, iNeuron}, 'format', 'timing', 'cycleTimes', mi_analysis.objBehav.data.cycleTimes.data);
                figure()
                for icycle = 1:length(mi_analysis.objBehav.data.cycleTimes.data)
                    spikes = spikeTimes(icycle, :);
                    plot(spikes, -1*icycle*ones(size(spikes)), 'bo')
                    hold on
                end
            end
        end
        
    end
end
