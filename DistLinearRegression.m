function [r_2,out_p] = DistLinearRegression(unit1,unit2)

%Works for count-count
r_2 = zeros(size(unit1,1),1);
% out_p = zeros(size(unit2));

    for j = 1:size(unit1,1)
        in_r = [ones(size(unit1(j,:)));unit1(j,:)];
        [filt,~,~,~,stats] = regress(unit2(j,:)',in_r');

        r_2(j) = stats(1);
        out_p(j,:) = (in_r'*filt)';
    end
   
end