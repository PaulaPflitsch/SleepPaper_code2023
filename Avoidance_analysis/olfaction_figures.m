%% 1) Preference index figure
% calls the function findPIs (preference indices) which calculates the PI for
% control and treatment groups
% this code plots the preference index as time spent on the side with low
% odor concentration vs the side with high odor concentration (binary)
% +1 signals high attraction, -1 signals high aversion
%For statistics, replace PI_post with PIcontrol_post for the control data

% run treatment(filtered) first, then load control (filtered) AllFish2_filtered2 & uncomment
% control and comment treatment

% treatment
%[PI_post,PI_mean_post,PI_SE_post] = findPIs(AllFish2_filtered,617) % 617 is the center point of the lane (pixels)
% control 
[PIcontrol_post,PIcontrol_mean_post,PIcontrol_SE_post] = findPIs(AllFish2_filtered,617) % 617 is the center point of the lane (pixels)

figure
hold on
bar(2,PI_mean_post,'FaceColor',[0.3,0.3,0.3]);
errorbar(2,PI_mean_post,PI_SE_post,'ok');
scatter(2*ones(2,size(PI_post,2)),PI_post,20,[0,0,0])

bar(1,PIcontrol_mean_post,'FaceColor',[0.6,0.6,0.6]);
errorbar(1,PIcontrol_mean_post,PIcontrol_SE_post,'ok');
scatter(ones(1,size(PIcontrol_post,2)),PIcontrol_post,20,[0,0,0])

axis([0 3 -1 1]);
xticks([1, 2]);
xticklabels({'Control', 'Treatment'});
ylabel('Preference Index');
title('Preference Index Comparison');

hold off


% Perform ranksum test
p = ranksum(PIcontrol_post, PI_post);
disp(['ranksum: ', num2str(p)]);

% Save p-value to a CSV file
save_directory = 'E:\Olfaction\Preference_index_cortisol\100ng_cortisol_15min'; % Replace with your desired path
csv_filename = 'PI_ranksum_results.csv';
full_csv_path = fullfile(save_directory, csv_filename);
csvwrite(full_csv_path, p);
disp(['Ranksum p-value saved to ', full_csv_path]);

% Perform two-sided independent t-test
[h, p_ttest, ci, stats] = ttest2(PIcontrol_post, PI_post);

% Optional: display t-test results
disp('Two-sided t-test results:');
disp(['p-value: ', num2str(p_ttest)]);
disp(stats);

%% 2) Bout rate Difference (NOT USED, instead jupyter lab bout rate)!
% calculates the difference in bout rate between low odor concentration and
% high odor concentration for each individual fish, then takes the mean and
% SEM

% run treatment (filtered) first, then load control (filtered) AllFish2_filtered2 & uncomment
% control and comment treatment

% treatment
%[total_boutrate, count_l, count_r, bpm_left, bpm_right,total_bout_left,total_bout_right,timeinterbout_left,timeinterbout_right,timebout_left,timebout_right,sum_timebout_left,sum_timeinterbout_left,sum_timeinterbout_right,sum_timebout_right,time_left,time_right,AllFish2,boutrate_left,boutrate_right,boutrateDiff]=Bout_Rate(AllFish2_filtered,2) % 2: calculates the bout rate within 2 cms of the no-odor-agarose and the odor-agarose
%meanboutrateDiff = mean(boutrateDiff,'omitnan');
% control
[total_boutrate, count_l, count_r, bpm_left, bpm_right,total_bout_left,total_bout_right,timeinterbout_left,timeinterbout_right,timebout_left,timebout_right,sum_timebout_left,sum_timeinterbout_left,sum_timeinterbout_right,sum_timebout_right,time_left,time_right,AllFish2,boutrate_left,boutrate_right,boutrateDiffC]=Bout_Rate(AllFish2_filtered,2)
meanboutRateDiffC = mean(boutrateDiffC,'omitnan');


% define y axis
y = [meanboutRateDiffC meanboutrateDiff];

figure 
hold on
bar(y);

diffSEC = std(boutrateDiffC,'omitnan')/sqrt(length(boutrateDiffC));
errorbar(1, meanboutRateDiffC,diffSEC, 'ok');
diffSE = std(boutrateDiff,'omitnan')/sqrt(length(boutrateDiff));
errorbar(2, meanboutrateDiff,diffSE, 'ok');
scatter(ones(1,size(boutrateDiffC,2)),boutrateDiffC,20,[0,0,0])
scatter(2*ones(1,size(boutrateDiff,2)),boutrateDiff,20,[0,0,0])

axis([0 3 -200 200]);
%axis([0 3 -30 45]); %without the raw data
xticks([1, 2]);
xticklabels({'Control', 'Treatment'});
ylabel('Bout Rate Difference (bouts/min)');
title("Bout rate difference(bouts/min spent within 2cm)");

hold off


% Perform ranksum test
p = ranksum(boutrateDiffC, boutrateDiff);
disp(['ranksum: ', num2str(p)]);

% Save p-value to a CSV file
save_directory = 'E:\Olfaction\Preference_index_cortisol\100ng_cortisol_15min'; % Replace with your desired path
csv_filename = 'BoutRate_ranksum_results.csv';
full_csv_path = fullfile(save_directory, csv_filename);
csvwrite(full_csv_path, p);
disp(['Ranksum p-value saved to ', full_csv_path]);

% Perform two-sided independent t-test
[h, p_ttest, ci, stats] = ttest2(boutrateDiffC, boutrateDiff);

% Optional: display t-test results
disp('Two-sided t-test results:');
disp(['p-value: ', num2str(p_ttest)]);
disp(stats);

%% 3) Turning angle
% load the control (filtered) and treatment (filtered) group seperately and then save them as a
% csv file in the same directory

[AllFish2_filtered2,bin_means_bout_rate_per_fish, mean_IBI_per_fish, bin_means_bout_rate, bin_means_IBI,IBIs_in_bin,bout_rates_in_bin,mean_IBI_per_fish_bin,mean_IBI_per_bin,mean_boutrate_per_bin,mean_boutrate_per_fish_bin,all_IBIs_with_postBout_x,all_bout_rates_with_postBout_x]  = binned_means_BoutRate_perFish_test2(AllFish2_filtered, 46) % 85 bins: 1mm bins in the 8.5 cm long lane; 55 bins for a 550 pixel lane = 10 pixel/lane --> 0.166cm/bin

% save file as csv, CHANGE directory
writematrix(bout_rates_in_bin,'E:\Olfaction\Preference_index_cortisol\100ng_cortisol_15min\bout_rates_control.csv');


% then plot in jupyter lab python

%% 4) plot scatterplot bout rate with linear fit
% plots all bout rates along the x-axis of the lanes (along relative odor
% concentration) and overlays a linear fit

% run seperately for control (filtered) and treatment (filtered) after
% running section 3)

sorted_br =sortrows(all_bout_rates_with_postBout_x,2); % sort based on position in lane
%sorted_IBI_ctrl=sortrows(IBI_matrix_ctrl,2);

% Normalize preBout_x to range from 360 (0%) to 870 (100%)
postBout_x = sorted_br(:, 2);  % Extract postBout_x
x_min = 330;
x_max = 900;
normalized_postBout_x = (postBout_x - x_min) / (x_max - x_min) * 100;


%make figure
figure 
hold on
x1 = sorted_br(:,2);   % x position in lane (2nd col)
y1 = sorted_br(:,1);   % bout rate in min in 1st column
scatter(normalized_postBout_x,y1,2, 'filled'); % Create a scatter plot

% Set x-axis limits from 0% to 100%
%xlim([0 100]);

% Add a title
title('bout rate vs. position melatonin');

% Fit a linear regression model to the data
p = polyfit(normalized_postBout_x, y1, 1); % Linear fit (1st degree polynomial)

% Evaluate the fitted polynomial over the x-range for plotting
%x_range = linspace(0, 100, 100);  % Generate 100 points from 0 to 100
y_fit = polyval(p, normalized_postBout_x);      % Compute the fitted y-values

% Hold the plot to add the regression line
hold on;

% Plot the regression line
plot(normalized_postBout_x, y_fit, 'r-', 'LineWidth', 1); % 'r-' is a red solid line

% add probability density function
%pd = fitdist(x1,'wbl');
%pdfEst = pdf(pd,x1);
%line(x1,pdfEst, 'Color','red');

% add a curve fit using the Weibull model
% custom equation (modelFun) created with Weibull fitter
%modelFun =  @(p,x) p(3) .* (x./p(1)).^(p(2)-1) .* exp(-(x./p(1)).^p(2));
%startingVals = [400 400 100];
%nlModel = fitnlm(x1,y1,modelFun,startingVals);
%line(x1,predict(nlModel,x1),'Color','r');


%ylim([0 20])
% Add a legend
legend('Data points', 'Linear regression', 'Location', 'best');
xlabel('X position in lane ');
ylabel('bout rate (bouts/s)');


%corr coeff
[R,P] = corrcoef(normalized_postBout_x,y1);

%safe coefficient and p value
Corr_coeff=R(2,1);
p_value = P(2,1);

%add text to plot
str = {'Correlation coefficient= ',Corr_coeff,'p= ', p_value};
text(340,18,str)


hold off