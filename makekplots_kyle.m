function makekplots_kyle(a,subgroups)

    figure
    for k=1:9
        clear est err frac
        for df=1:10
            est(df) = a.arrMIcore{subgroups,1}.mi_data{df+10*(k-1),1};
            err(df) = a.arrMIcore{subgroups,1}.mi_data{df+10*(k-1),2};
            frac(df) = a.arrMIcore{subgroups,1}.mi_data{df+10*(k-1),3};
        end

        subplot(3,3,k)
        hold on
        scatter(frac,est,20,'filled','b');
        errorbar(frac,est,sqrt(err),'b','linewidth',1)
    %     xlabel('data fracs')
    %     ylabel('MI Estimate')
        ylim([0,0.4])
        title(['k = ',num2str(k)])
    end