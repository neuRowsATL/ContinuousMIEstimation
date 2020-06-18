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

%% Build cc MI
clear d b cc
d = mi_data_neural('test', 'verbose', verbose_level);
b = mi_data_pressure('test', 'verbose', verbose_level);
add_cycleTimes(b, cycle_times, str_cycles, 30000);
add_spikes(d, unit2, str_unit2, 30000, 'unit1');
add_spikes(d, unit1, str_unit1, 30000, 'unit2');

cc = test_analysis("count_count",0,d,b,verbose_level);
base_values = [cc.arrMIcore{1,1}.x;cc.arrMIcore{1,1}.y]; %Original values
cc.calcMIs();
cc.getMIs();

cc_output = DistLinearRegression(cc);

est_act(1) = cc.arrMIcore{1,4};
est_act(2) = cc.arrMIcore{1,5};
est_act(3) = cc_output.r_2;

%% Reparameterization or Jitter
%Currently set to reparameterize the jittered discrete count data for both
%units
trials = 10;
est_reparam = compare_mi_lr(cc,base_values,"reparameter",trials);
%% Reparam Figure
figure
subplot(2,1,1)

hold on
h=fill([1,trials,trials,1],[est_act(1)-est_act(2),est_act(1)-est_act(2),est_act(1)+est_act(2),est_act(1)+est_act(2)],'b','linestyle','none');
h.FaceAlpha=0.1;
plot([1,trials],[est_act(1),est_act(1)],'linewidth',2)
errorbar(1:trials,est_reparam(:,1),est_reparam(:,2),'.')
hold off

legend('Actual MI Estimate with error')
title(['Mean = ',num2str(round(mean(est_reparam(:,1)),3)),' STD = ',num2str(round(std(est_reparam(:,1)),3))])
ylabel('KSG MI')
set(gca,'fontsize',14)
ylim([0.02,0.12])
xlim([1,trials])

subplot(2,1,2)
hold on
plot([1,trials],[est_act(3),est_act(3)],'linewidth',2)
scatter(1:trials,est_reparam(:,3),14,'filled') 
hold off
legend('Actual R^2')
title(['Mean = ',num2str(round(mean(est_reparam(:,3)),3)),' STD = ',num2str(round(std(est_reparam(:,3)),3))])
xlabel('Trial')
ylabel('R^2 (Linear Correlation)')
set(gca,'fontsize',14)

sgtitle('Jittered and Reparamterized')
%% Power
%Raise both unit counts to the same power
trials = 20;
powers = linspace(0.1,3,trials);
est_power = compare_mi_lr(cc,base_values,"power",trials);
%% Power Figure
figure
subplot(2,1,1)
hold on
plot([powers(1),powers(end)],[est_act(1),est_act(1)],'linewidth',2)
h=fill([powers(1),powers(end),powers(end),powers(1)],[est_act(1)-est_act(2),est_act(1)-est_act(2),est_act(1)+est_act(2),est_act(1)+est_act(2)],'b','linestyle','none');
h.FaceAlpha=0.1;
errorbar(powers,est_power(:,1),est_power(:,2),'.')
hold off
title(['Mean = ',num2str(round(mean(est_power(:,1)),3)),' STD = ',num2str(round(std(est_power(:,1)),3))])
legend('Actual MI Estimate with error')
ylabel('KSG MI')
set(gca,'fontsize',14)
ylim([0.02,0.12])

subplot(2,1,2)
hold on
plot([powers(1),powers(end)],[est_act(3),est_act(3)],'linewidth',2)
scatter(powers,est_power(:,3),14,'filled')
hold off
legend('Actual R^2')
title(['Mean = ',num2str(round(mean(est_power(:,3)),3)),' STD = ',num2str(round(std(est_power(:,3)),3))])
xlabel('Power (data was raised to this power)')
ylabel('R^2 (Linear Correlation)')
set(gca,'fontsize',14)

sgtitle('Raised to a Power')
%% Log10
trials = 10;
est_log = compare_mi_lr(cc,base_values,"log",trials);
%% Log10 Figure
figure
hold on
subplot(2,1,1)
scatter(1:trials,est_log(:,1),14,'filled')
errorbar(1:trials,est_log(:,2))
title(['Mean = ',num2str(round(mean(est_log(:,1)),3)),' STD = ',num2str(round(std(est_log(:,1)),3))])
ylabel('KSG MI')
set(gca,'fontsize',14)

subplot(2,1,2)
scatter(1:trials,est_log(:,3),14,'filled')
title(['Mean = ',num2str(round(mean(est_log(:,3)),3)),' STD = ',num2str(round(std(est_log(:,3)),3))])
xlabel('Trial')
ylabel('R^2 (Linear Correlation)')
set(gca,'fontsize',14)

sgtitle('Log')
%%
function mi_est = compare_mi_lr(cc,base_values,change,trials)
    powers = linspace(0.1,3,trials);
    for j=1:trials

            if change == "reparameter"
                rng shuffle %reinsert variability
                
                random_values = randn(2,length(base_values));
                cc_reparam = [reparameterize_data(base_values(1,:)+random_values(1,:)),reparameterize_data(base_values(2,:)+random_values(2,:))]';

                change_data(cc,cc_reparam)
                
            elseif change == "power"
                change_data(cc,base_values.^powers(j))
                
            elseif change == "log"
                change_data(cc,log10(base_values+1))
            end

        cc.arrMIcore{1,1}.mi_data=[];
        cc.arrMIcore{1,1}.opt_k=1;

        cc.calcMIs();
        cc.getMIs();

        cc_output = DistLinearRegression(cc);

        mi_est(j,1) = cc.arrMIcore{1,4};
        mi_est(j,2) = cc.arrMIcore{1,5};
        mi_est(j,3) = cc_output.r_2;
    end
 
end
    
function change_data(cc,values)
    cc.arrMIcore{1,1}.x = values(1,:);
    cc.arrMIcore{1,1}.y = values(2,:);
end

function a = test_analysis(analysis,reparam,d,b,verbose_level)

if analysis == "timing_count"
    a = calc_timing_count(d, b, {'unit1' , 'unit2'}, 'verbose', verbose_level,'reparam',reparam);
elseif analysis == "timing_timing"
    a = calc_timing_timing(d, b, {'unit1' , 'unit2'}, 'verbose', verbose_level,'reparam',reparam);
elseif analysis == "count_count"
    a = calc_count_count(d, b, {'unit1' , 'unit2'}, 'verbose', verbose_level,'reparam',reparam);
end
a.buildMIs();

end