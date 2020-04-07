function output = DistLinearRegression(obj)
% figure
    for j=1:size(obj.arrMIcore,1)
        
        unit_1 = obj.arrMIcore{j,1}.x;
        unit_2 = obj.arrMIcore{j,1}.y;
        
        if size(obj.arrMIcore,1)==1
            unit_1 = unit_1';
            unit_2 = unit_2';
        end
        
        x_plot = unit_1 + 0.2*rand(size(unit_1));
        y_plot = unit_2 + 0.2*rand(size(unit_2));
        
        in_r = [ones(size(unit_1,1),1),unit_1];
        [filt,~,r,~,stats] = regress(unit_2,in_r);

        r_2 = stats(1);
        i_true = -0.5 * log2(1-r_2);
        y_pred = in_r*filt;
%         
%         figure
%         for k = 1:size(unit_1,2)
%             subplot(3,5,k)
%             hold on
%             plot(x_plot(:,k),y_plot,'x')
%             plot(x_plot(:,k),y_pred,'x')
%             hold off
%         end
%         figure
%         subplot(3,5,j)
%         hold on
%         plot(y_plot,y_pred,'.','linewidth',2)
%         plot([min([y_plot,y_pred]),max([y_plot,y_pred])],[min([y_plot,y_pred]),max([y_plot,y_pred])],'k--','linewidth',2)
%         hold off
%         title(['Subgroup ',num2str(j), ': I_t_r_u_e = ',num2str(round(i_true,3))])
% %         title(['Count-Count: I_t_r_u_e = ',num2str(round(i_true,2))])
%         legend(['R^2 = ',num2str(round(r_2,3))])
%         axis('square')
% %         xlim([0,60])
% %         ylim([0,60])
%         xlabel('Actual Unit 1 Count')
%         ylabel('Predicted Unit 1 Count')
%         set(gca,'fontsize', 18)
       
        
%         figure
%         subplot(3,1,1)
%         hold on
%         plot(x_plot(:,end),y_plot,'x','linewidth',1)
%         h=plot(x_plot(:,end),y_pred,'linewidth',2)
%         legend(h,['Regression, R^2= ',num2str(round(r_2,2))])
%         hold off
%         title(['R^2 = ',num2str(r_2)])
%         xlabel('Unit 2 Spike Count')
%         ylabel('Unit 1 Spike Count')
%                 set(gca,'fontsize', 18)
%                         xlim([0,60])
%         ylim([0,60])

%                 
%         subplot(3,1,2)
%         plot(r)
%         ylabel('Residual')
%         xlabel('Cycle')
%         
%         subplot(3,1,3)
%         plot(filt)
%         ylabel('Filter Coefficient')
%         xlabel('Spike Number')
        
        output(j).filter = filt;
        output(j).xforfilt = in_r;
        output(j).r_2 = r_2;
        output(j).i_true = i_true;
        output(j).out_y = y_pred;
        output(j).real_x = unit_1;
        output(j).real_y = unit_2;

    end
end