% Initiate files to use in for loop
currentDirectory = uigetdir;

fileString = fullfile(currentDirectory, '*.mat');

files = dir(fileString);

%% FOR MI_xx files: Initiate for loop

for iFile = 1:length(files)
% Set filename
fileName = files(iFile).name;

figFileName = strcat(fileName(1:end-4),'_dataFracFig.png');

% Load file.
 load(fileName);
 
% Define object  
if exist('mi_isi')
    obj = mi_isi;
    titleStart = 'isi';
elseif exist('mi_cc')
    obj = mi_cc;
    titleStart = 'cc';
elseif exist('mi_tc')
    obj = mi_tc;
    titleStart = 'tc';
elseif exist('mi_tt')
    obj = mi_tt;
    titleStart = 'tt';
elseif exist('mi_cb')
    obj = mi_cb;
    titleStart = 'cb';
elseif exist('mi_tb')
    obj = mi_tb;
    titleStart = 'tb';
elseif exist('mi_tcb')
    obj = mi_tcb;
    titleStart = 'tcb';
elseif exist('mi_ttb')
    obj = mittb;
    titleStart = 'ttb';
end
 
    
 % Next, find the number of subgroups in the calculation: 
 n_subgroups = size(obj.arrMIcore,1);
 
 % Set size of subplot
 
 nRows_plot = ceil(n_subgroups / 5);
 
 newFig = figure;
 
 % NOTE- arrMIcore is supposed to document the k-value used, but it seems
 % that this is currently not happening.We can get the optimal k-value from
 % each mi_ksg_core object in arrMIcore.
 
 % Iterate through each subgroup and make data fraction subplot. 
 for iSubgroup = 1:n_subgroups
     
    % Initiate subplot axes 
    %ax = subplot(nRows_plot,5,iSubgroup);
    ax = subplot(nRows_plot,1,iSubgroup);
    % Identify k value for this subgroup
    coreObj = obj.arrMIcore(iSubgroup,1);
    opt_k = coreObj{1,1}.opt_k;
    
    % Find the data fraction calculations that correspond to the optimized
    % k value
    mi_calcs = cell2mat(coreObj{1,1}.mi_data);
    
    dataFracIdx = find(mi_calcs(:,4) == opt_k);
    
    % Set variables for plot
    xs = mi_calcs(dataFracIdx,3);
    ys = mi_calcs(dataFracIdx,1);
    err = mi_calcs(dataFracIdx,2);
    
    % Make plot
    errorbar(ax,xs,ys,err,'-b','Marker','.','MarkerSize',15);
    
    % Set x limits
    xlim([min(xs)*0.8 max(xs)*1.1]);
    
    % Get rid of x labels
    set(gca,'xtick',[])

    % Set subplot title as the percentage of data in the subgroup
    percentData = num2str(obj.arrMIcore{iSubgroup,2}*100);

    
    figTitle = strcat(percentData,'%');
    
    title(figTitle)
    
 end

 
 set(gca,'xtickMode', 'auto')
 
 saveas(newFig,figFileName);
 
end
%% FOR MI_xx_jitter files: Initiate for loop

for iFile = 1:length(files)
% Set filename
fileName = files(iFile).name;

% Load file.
 load(fileName);
 
 dataStruct = MIs;
 
 dataFields = fields(dataStruct);
 
 for iField = 1:length(dataFields)
     
 
 % Set MI object (NOTE - this will need to be customized for each file)
 obj = dataStruct(1).(dataFields{iField});
 
 
 % Next, find the number of subgroups in the calculation: 
 n_subgroups = size(obj.arrMIcore,1);
 
 % Set size of subplot
 
 nRows_plot = ceil(n_subgroups / 5);
 
 newFig = figure;
 
 % NOTE- arrMIcore is supposed to document the k-value used, but it seems
 % that this is currently not happening.We can get the optimal k-value from
 % each mi_ksg_core object in arrMIcore.
 
 % Iterate through each subgroup and make data fraction subplot. 
 for iSubgroup = 1:n_subgroups
     
    % Initiate subplot axes 
    %ax = subplot(nRows_plot,5,iSubgroup);
    ax = subplot(nRows_plot,1,iSubgroup);
    % Identify k value for this subgroup
    coreObj = obj.arrMIcore(iSubgroup,1);
    opt_k = coreObj{1,1}.opt_k;
    
    % Find the data fraction calculations that correspond to the optimized
    % k value
    mi_calcs = cell2mat(coreObj{1,1}.mi_data);
    
    dataFracIdx = find(mi_calcs(:,4) == opt_k);
    
    % Set variables for plot
    xs = mi_calcs(dataFracIdx,3);
    ys = mi_calcs(dataFracIdx,1);
    err = mi_calcs(dataFracIdx,2);
    
    % Make plot
    errorbar(ax,xs,ys,err,'-b','Marker','.','MarkerSize',15);
    
    % Set x limits
    xlim([min(xs)*0.8 max(xs)*1.1]);
    
    % Get rid of x labels
    set(gca,'xtick',[])
    
    % Set subplot title as the percentage of data in the subgroup
    %percentData = num2str(obj.arrMIcore{iSubgroup,2}*100);

    
    %figTitle = strcat(percentData,'%');
    
    %title(figTitle)

 end

 figFileName = strcat(fileName(1:end-4),dataFields{iField},'_dataFracFig.png');
 
 
 
 set(gca,'xtickMode', 'auto')
 
 saveas(newFig,figFileName);
 
 close all
 
 end
 
end

%% MAKE DATA FRACTION PLOTS FOR EACH K VALUE FOR EACH SUB-GROUP
% Initiate files to use in for loop
currentDirectory = uigetdir;

fileString = fullfile(currentDirectory, '*.mat');

files = dir(fileString);
%%
for iFile = 1:length(files)
% Set filename
 fileName = files(iFile).name;

% Load file.
 load(fileName);
 
 %Set MI object (NOTE - this will need to be customized for each file)
% Define object  
if exist('mi_isi')
    obj = mi_isi;
    titleStart = 'isi';
elseif exist('mi_cc')
    obj = mi_cc;
    titleStart = 'cc';
elseif exist('mi_tc')
    obj = mi_tc;
    titleStart = 'tc';
elseif exist('mi_tt')
    obj = mi_tt;
    titleStart = 'tt';
elseif exist('mi_cb')
    obj = mi_cb;
    titleStart = 'cb';
elseif exist('mi_tb')
    obj = mi_tb;
    titleStart = 'tb';
elseif exist('mi_tcb')
    obj = mi_tcb;
    titleStart = 'tcb';
elseif exist('mi_ccb')
    obj = mi_ccb;
    titleStart = 'ccb';
elseif exist('mittb')
    obj = mittb;
    titleStart = 'ttb';
end


 vizObj = mi_ksg_viz;
% MI_tc
subgroups = obj.arrMIcore;
for iSubgroup = 1:size(subgroups,1)
    %
    % We will make one plot for every value of k for every subgroup
    % Each sub-group gets a single large plot with k + 1 different sub plots 
    
    % First, instantiate the large subplot and plot the k - dependence
    
    newFig = figure();
    
    ax(1) = subplot(3,3,1);  
    
    % Set the core object
    objCore = obj.arrMIcore{iSubgroup,1};
    
    % Plot the k-dependence
    rplot = vizObj.plot_k_dependence(objCore,newFig, ax(1));
    
    % Set subplot title as the percentage of data in the subgroup
    percentData = num2str(obj.arrMIcore{iSubgroup,2}*100);

    figNum = num2str(iSubgroup);  
    figTitle1 = strcat(titleStart,' ',fileName(24:28), '-', figNum);
    figTitle2 = strcat(percentData, '%');

    title({ figTitle1 figTitle2});
    % Iterate through each possible value of k to make different data
    % fraction plots
    
    for k = 1:9
        % Set new axes for each subplot
        ax(k+1) = subplot(4, 3, k + 3);
        
        % Plot data fraction
        rplot = vizObj.plot_data_fraction(objCore,k,newFig,ax(k+1));
    end
    
linkaxes(ax,'y')    
figFileName = strcat(fileName(1:end-4),'_dataFracFig_', figNum ,'.png');

saveas(newFig,figFileName);

close all

end
clearvars -except currentDirectory files fileString
end