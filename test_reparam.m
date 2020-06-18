clear all
close('all')

verbose_level = 5;
load('TestData/20200127_bl21lb21_spikedata.mat');
unit1 = spikedata.unit1;
unit2 = spikedata.unit3;
unit3 = spikedata.unit4;

str_unit1 = 'TestData/20200127_bl21lb21_spikedata.mat/spikedata.unit1';
str_unit2 = 'TestData/20200127_bl21lb21_spikedata.mat/spikedata.unit2';
str_unit3 = 'TestData/20200127_bl21lb21_spikedata.mat/spikedata.unit3';

cycle_times = [spikedata.pressure.Ontime(1:end-1,1) spikedata.pressure.Ontime(2:end,1)]; % Needs to be N x 2 matrix of [on off] x N
str_cycles = 'TestData/20200127_bl21lb21_spikedata.mat/spikedata.pressure.Ontime';

addpath('C:\Users\kthom88\OneDrive - Georgia Institute of Technology\Year 1 PhD\Lab - Sober\ContinuousMIEstimation\kraskovStoegbauerGrassberger')
addpath('C:\Users\Kyle\OneDrive - Georgia Institute of Technology\Year 1 PhD\Lab - Sober\ContinuousMIEstimation\kraskovStoegbauerGrassberger')
%% Run MIs
clear d b cc
d = mi_data_neural('test', 'verbose', verbose_level);

b = mi_data_pressure('test', 'verbose', verbose_level);
add_cycleTimes(b, cycle_times, str_cycles, 30000);

add_spikes(d, unit2, str_unit2, 30000, 'unit1');
add_spikes(d, unit1, str_unit1, 30000, 'unit2');

% for j = 1:4
    cc = test_analysis("count_count",0,d,b,verbose_level);
    base = [cc.arrMIcore{1,1}.x;cc.arrMIcore{1,1}.y];
%%

%     x = reparameterize_data(x+randn(1,length(x)))';
%     y = reparameterize_data(y+randn(1,length(y)))';
% 
%     x = x+randn(1,length(x));
%     y = y+randn(1,length(y));
% 
%     x = x.^0.75;
%     y = y.^0.75;

for j=1:10

    if j == 1
        change_data(cc,base)
    else
        rng shuffle
        randn
        cc_reparam = [reparameterize_data(base(1,:)+randn(1,length(base))),reparameterize_data(base(2,:)+randn(1,length(base)))]';
        change_data(cc,cc_reparam)
    end

    cc.arrMIcore{1,1}.mi_data=[];
    cc.arrMIcore{1,1}.opt_k=1;

    cc.calcMIs();
    cc.getMIs();
    
    cc_output = DistLinearRegression(cc);
    
    mi_est(j,1) = cc.arrMIcore{1,4};
    mi_est(j,2) = cc_output.r_2;
end
%%
figure
subplot(2,1,1)
plot(mi_est(:,1),'.-')
subplot(2,1,2)
plot(mi_est(:,2),'.-')


% end
% tc_noreparam = test_analysis("timing_count",0,d,b,verbose_level);
% tc_reparam = test_analysis("timing_count",2,d,b,verbose_level);
% % 
% tt_noreparam = test_analysis("timing_timing",0,d,b,verbose_level);
% tt_reparam = test_analysis("timing_timing",2,d,b,verbose_level);

% save('test_reparam_results','tc_noreparam','tc_reparam')%,'tt_noreparam','tt_reparam')
%% Testing Regression
close all
cc_output = DistLinearRegression(cc);
% tc_noreparam_output = DistLinearRegression(tc_noreparam);
% tc_reparam_output = DistLinearRegression(tc_reparam);
%% Histogram
figure
hold on
histogram(cc.arrMIcore{1, 1}  .y)
histogram(cc.arrMIcore{1, 1}  .x)
hold off
xlim([0,40])
xlabel('Spike Count')
ylabel('Cycles')
legend('Unit 1','Unit 2')
set(gca,'fontsize',18)

%% Information Ratio
for j=1:14
    i_true(j) = tc_reparam_output(j).i_true
end
figure
hold on
plot(1:14,cell2mat(mi_tc.a.arrMIcore(:,4) )./i_true','.-','linewidth',2)
plot([0,1],[14,1],':')
hold off
xlabel('Subgroup')
ylabel('Ratio')
title('KSG MI to I_t_r_u_e Ratio')
set(gca,'fontsize',18)

%% Testing Regression and Assumptions
clear n r c
for j=1:10000
change_data(cc,base);

% cc_normal = base;
% n(j)=abs(randn);
% cc_pwr = base.^n(j);
% cc_log = log10(base+1);
cc_reparam = [reparameterize_data(base(1,:)+randn(1,length(base))),reparameterize_data(base(2,:)+randn(1,length(base)))]';

% change_data(cc,cc_reparam)
% cc_output = DistLinearRegression(cc); %Make sure plots are turned off
% r(j)=cc_output.r_2;

c=corrcoef(cc_reparam');
r(j) = c(2).^2;
end
figure
histogram(r)
% figure
% plot(n,r,'.')
mean(r)
%%
clear rep_norep
% diff_rep_norep = zeros(14,1)
rep_norep = cell2mat(tc_reparam.arrMIcore(:,4))  - cell2mat(tc_noreparam.arrMIcore(:,4));

figure
hold on
errorbar(1:14,zeros(1,14),cell2mat(tc_reparam.arrMIcore(:,5)),'.','linewidth',2)
errorbar((1:14)+0.1,zeros(1,14),cell2mat(tc_noreparam.arrMIcore(:,5)),'.','linewidth',2)
scatter(cell2mat(tc_reparam.arrMIcore(:,4)),cell2mat(tc_noreparam.arrMIcore(:,4)),'.','linewidth',3)
title('reparameterized - not reparameterized')
xlabel('Subgroup')
ylabel('MI')
% legend('Error - Reparameterized','Error - Not Reparameterized','MI Difference')
% ylim([-
%%
% close all
% display_results(tc_noreparam)

knn_mi = cell2mat(tc_noreparam.arrMIcore(:,4));
r2_mi = cell2mat(tc_noreparam.arrMIcore(:,7));

figure
hold on
scatter(knn_mi,r2_mi,20,'filled')
plot([0,.5],[0,.5],':','linewidth',2)
xlabel('MI from k-nearest neighbors')
ylabel('MI from R^2')
axis('square')
xlim([-0.2,0.8])
ylim([-0.2,0.8])
title('Reparameterized Timing-Count Units 2-3')
%%
function a = test_analysis(analysis,reparam,d,b,verbose_level)

if analysis == "timing_count"
    a = calc_timing_count(d, b, {'unit1' , 'unit2'}, 'verbose', verbose_level,'reparam',reparam);
elseif analysis == "timing_timing"
    a = calc_timing_timing(d, b, {'unit1' , 'unit2'}, 'verbose', verbose_level,'reparam',reparam);
elseif analysis == "count_count"
    a = calc_count_count(d, b, {'unit1' , 'unit2'}, 'verbose', verbose_level,'reparam',reparam);
end
a.buildMIs();
% 
% a.calcMIs();
% a.getMIs();
end

function display_results(a)
    for subgroups = 1:size(a.arrMIcore,1)
        makekplots_kyle(a,subgroups)
        sgtitle(['Subgroup ',num2str(subgroups)])
    end
end

function change_data(cc,values)
    cc.arrMIcore{1,1}.x = values(1,:);
    cc.arrMIcore{1,1}.y = values(2,:);
end