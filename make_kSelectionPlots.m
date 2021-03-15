%% GET FILES OF INTEREST

currentDirectory = uigetdir;

fileString = fullfile(currentDirectory, '*.mat');

files = dir(fileString);

%% MAKE K PLOTS FOR EACH FILE OF INTEREST
for iFile = 1:length(files)
% Set filename
 fileName = files(iFile).name;

% Load file.
 load(fileName);
 
 %Set MI object (NOTE - this will need to be customized for each file)
% Define object  
if exist('a_isi')
    obj = a_isi;
    titleStart = 'isi';
elseif exist('a_cc')
    obj = a_cc;
    titleStart = 'cc';
elseif exist('a_tc')
    obj = a_tc;
    titleStart = 'tc';
elseif exist('a_tt')
    obj = a_tt;
    titleStart = 'tt';
elseif exist('a_cb')
    obj = a_cb;
    titleStart = 'cb';
elseif exist('a_tb')
    obj = a_tb;
    titleStart = 'tb';
elseif exist('a_tcb')
    obj = a_tcb;
    titleStart = 'tcb';
elseif exist('a_ccb')
    obj = a_ccb;
    titleStart = 'ccb';
elseif exist('a_ttb')
    obj = a_ttb;
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