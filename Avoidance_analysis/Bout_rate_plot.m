% Bout rate figures

% call all variables
[total_boutrate, count_l, count_r, bpm_left, bpm_right,total_bout_left,total_bout_right,timeinterbout_left,timeinterbout_right,timebout_left,timebout_right,sum_timebout_left,sum_timeinterbout_left,sum_timeinterbout_right,sum_timebout_right,time_left,time_right,AllFish2,boutrate_left,boutrate_right,boutrateDiffC, boutrate_left_for_fish_i]=Bout_Rate(AllFish2_filtered,2)

%make the figures for bout rate per total time
meanboutRateDiffC = mean(boutrateDiffC,'omitnan');
meanboutrateDiff = mean(boutrateDiff,'omitnan');
y = [meanboutRateDiffC meanboutrateDiff];
figure 
hold on
bar(y);
title("Bout rate difference(bouts/min spent within 2cm)");
diffSEC = std(boutrateDiffC,'omitnan')/sqrt(length(boutrateDiffC));
errorbar(1, meanboutRateDiffC,diffSEC, 'ok');
diffSE = std(boutrateDiff,'omitnan')/sqrt(length(boutrateDiff));
errorbar(2, meanboutrateDiff,diffSE, 'ok');
scatter(ones(1,size(boutrateDiffC,2)),boutrateDiffC,20,[0,0,0])
scatter(2*ones(1,size(boutrateDiff,2)),boutrateDiff,20,[0,0,0])
%axis([0 3 -150 150]);

% bout rate over time spent within 2cm far from cadaverine
meanbr_left = mean(boutrate_left,"omitnan");
y=meanbr_left;
figure 
hold on
bar(y);
title("far from to cadaverine");
diffSE_l = std(boutrate_left,'omitnan')/sqrt(length(boutrate_left));
errorbar(1, meanbr_left,diffSE_l, 'ok');
scatter(ones(1,size(boutrate_left,2)),boutrate_left,20,[0,0,0])
%axis([0 2 0 200]);

% bout rate over time spent within 2cm close to cadaverine
meanbr_right = mean(boutrate_right,"omitnan");
y=meanbr_right;
figure 
hold on
bar(y);
title("close to cadaverine");
diffSE_l = std(boutrate_right,'omitnan')/sqrt(length(boutrate_right));
errorbar(1, meanbr_right,diffSE_l, 'ok');
scatter(ones(1,size(boutrate_right,2)),boutrate_right,20,[0,0,0])
%axis([0 2 0 200]);


% statistical test
%pboutRate = ranksum(boutrate_left, boutrate_right);
pboutRate = ranksum(boutrateDiff, boutrateDiffC);

[h,p,ci,stats] = ttest2(boutrateDiff, boutrateDiffC)

