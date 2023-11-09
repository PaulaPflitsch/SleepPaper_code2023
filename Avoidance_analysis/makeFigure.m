%For statistics, replace PI_post with PIcontrol_post for the control data
[PIcontrol_post,PI_mean_post,PI_SE_post] = findPIs(AllFish2,617);%changed center point to be correct 
figure
hold on
bar(1,PI_mean_post,'FaceColor',[0.6,0.6,0.6]);
errorbar(1,PI_mean_post,PI_SE_post,'ok');
scatter(ones(1,size(PI_post,2)),PI_post,20,[0,0,0])
axis([0 2 -1 1])

%[file,path] = uiputfile('*.pdf');
%save(file)

%[filename, pathname] = uiputfile('*.pdf');
%path_file=fullfile(pathname,filename)
%a=fopen(path_file,'wt');
%fprintf(a,'write something.');
%%
%down is away from cadaverine
for n = 1:size(AllFish2,2)
    binnedY = discretize(AllFish2(n).preBout_y,[000:100:1200]);
    for i = unique(binnedY)'
        tmp = (binnedY == i);
        angles = abs(AllFish2(n).delta(tmp));
        MeanAngles(n,i) = mean(angles);
    end
    
end


figure 
hold on
plot(mean(MeanAngles),'LineWidth',2,'Color',[0,0,0])
for i = 1:size(AllFish2,2)
    plot(MeanAngles,'Color',[0.6,0.6,0.6])
end

%9.5 cm/570 pixels, 60 pixels/cm
%For statistics, replace fishSpeedDiff with fishSpeedDiffC for the control data
[lSE, rSE, fishLeft, fishRight, fishSpeedDiffC, meanSpeed, meanLeft, meanRight, r] = Nadine_Speed(AllFish2,2);


meanSpeedDiffC = mean(fishSpeedDiffC,'omitnan');
meanSpeedDiff = mean(fishSpeedDiff,'omitnan');
y = [meanSpeedDiffC, meanSpeedDiff];
figure 
hold on
bar(y);
%boxplot([meanSpeedDiffC,meanSpeedDiff])
diffSEC = std(fishSpeedDiffC,'omitnan')/sqrt(length(fishSpeedDiffC));
errorbar(1, meanSpeedDiffC,diffSEC, 'ok');
diffSE = std(fishSpeedDiff,'omitnan')/sqrt(length(fishSpeedDiff));
errorbar(2, meanSpeedDiff,diffSE, 'ok');
scatter(ones(1,size(fishSpeedDiffC,2)),fishSpeedDiffC,20,[0,0,0],'jitter','on', 'jitterAmount',0.05)
scatter(2*ones(1,size(fishSpeedDiff,2)),fishSpeedDiff,20,[0,0,0],'jitter','on', 'jitterAmount',0.05)
axis([0 3 -2 2.5]);

pspeed = ranksum(fishSpeedDiffC, fishSpeedDiff);
%two sided independent t-test
[h,p,ci,stats] = ttest2(fishSpeedDiffC,fishSpeedDiff);

%y = [meanLeft meanRight];
%figure 
%hold on
%bar(y);
%leftSE = std(fishLeft,'omitnan')/sqrt(length(fishLeft));
%righttSE = std(fishRight,'omitnan')/sqrt(length(fishRight));
%errorbar(1,meanLeft,leftSE,'ok');
%errorbar(2,meanRight,righttSE,'ok');
%axis([0 3 0 4])
%[h,pspeed] = ttest(fishLeft,fishRight);
%controlRight = fishRight;
%rightP = ranksum(controlRight, fishRight);

p = ranksum(PIcontrol_post, PI_post);
%two sided independent t-test
[h,p,ci,stats] = ttest2(PIcontrol_post,PI_post);

% figure 
% hold on
% scatter(PI_post, fishSpeedDiff)
% ylim([-0.1 0.5])
%xlabel('PI');
%ylabel('Speed');

%[speedProbSleep, toCadSleep, awayCadSleep, locProbSleep] = Nadine_SpeedProb(AllFish2)
%[speedProbControl, toCadControl, awayCadControl, locProbControl] = Nadine_SpeedProb(AllFish2)
 vals = [speedProbControl; speedProbSleep];
 %valsProb = [locProbControl; locProbSleep];
% valsTo = [toCadControl; toCadSleep];
% valsAway = [awayCadControl; awayCadSleep];
% x = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5];
 x = [1,2,3,4,5,6,7,8,9];
 figure 
 hold on
 bar(x,vals)
% 
% figure 
% hold on
% bar(x, valsTo)
% ylim([0 7])
% 
% figure 
% hold on
% bar(x, valsAway)
% 
 %figure 
 %hold on
 %bar(x, valsProb)
