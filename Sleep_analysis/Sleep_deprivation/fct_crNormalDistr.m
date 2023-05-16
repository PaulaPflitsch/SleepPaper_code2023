function [meanAct, sumAct, cumsumAct, tg, wt, allWells] = fct_crNormalDistr(a, expName, geno, allWells,  Group1)

% meanAct = mean active seconds per fish
% sumAct = summed active seconds per fish
% cumsumAct = cumulative active seconds per fish
% all = all wells with a fish, exclude wells w/o fish
allWells = allWells
%% make different groups based on genotype data
% allWells
% plot each fish individually in heatmap
figure('Color',[1 1 1]);
colormap('hot'); % default bone gray hot hsv jet winter copper pink
ax=imagesc(a(:,allWells)');
% set(ax,'box','off');
% colorbar('location','southoutside');
xlabel('time in seconds');
ylabel('individual fish activity (active seconds/min)');
title('Activity "heatmap"');
fct_suptitle([' Overview activity allWells: ' expName]);
% datetick('x','HHPM');

% include time in x axis, plot against newtimenum, datetick
% doesnt work to plot against newtimenum....


%% visualize mean activity per well position in heat map
% make vectors for each row in plate
% clear R1 R2 R3 R4 R5 R6 R7 R8 plate
% R1 = linspace(1,89,12);
% R2 = linspace(2,90,12);
% R3 = linspace(3,91,12);
% R4 = linspace(4,92,12);
% R5 = linspace(5,93,12);
% R6 = linspace(6,94,12);
% R7 = linspace(7,95,12);
% R8 = linspace(8,96,12);
% plate = [R1; R2; R3; R4; R5; R6; R7; R8]
% next time something with reshape and linspace 1 to 96

clear meanAct
meanAct = mean(a(:,allWells));

% reshape meanAct according to plate
meanActpl = reshape(meanAct, 8, [])

% visualize mean activity on plate position
figure('Color', [1 1 1]);
colormap('hot');
ax = imagesc(meanActpl)

title('Mean Activity per well position 24well plate');
fct_suptitle([' Overview activity Group1: ' expName]);

% visualize each genotype differently??? 


%% Group1
% plot each fish individually in heatmap

figure('Color',[1 1 1]);
colormap('hot'); % default bone gray hot hsv jet winter copper pink
%ax=imagesc(a(:,Group1)');
%paula: what is group1? Try: all of a
ax=imagesc(a(:,:)');
% set(ax,'box','off');
% colorbar
xlabel('time in seconds');
ylabel('individual fish, activity white');
title('Activity "heatmap"');
fct_suptitle([' Overview activity Group1: ' expName]);
% datetick('x','HHPM');

% include time in x axis, plot against newtimenum, datetick
% doesnt work to plot against newtimenum....

%% Group2
% plot each fish individually in heatmap
% figure('Color',[1 1 1]);
% colormap('hot'); % default bone gray hot hsv jet winter copper pink
% ax=imagesc(a(:,Group2)');
% set(ax,'box','off');
% colorbar
% xlabel('time in seconds');
% ylabel('individual fish, activity white');
% title('Activity "heatmap"');
% fct_suptitle([' Overview activity Group2: ' expName]);
% datetick('x','HHPM');

% include time in x axis, plot against newtimenum, datetick
% doesnt work to plot against newtimenum....

%% 96 subplots ...
% all = [1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96]
% figure('Color',[1 1 1]);
% ax(1)= subplot(96,1,1);
% bar(logical(activity(:,1)));
% ax(2)= subplot(96,1,2);
% bar(logical(activity(:,2)));
% ax(3)= subplot(96,1,3);
% bar(logical(activity(:,3)));
% ax(4)=subplot(96,1,4);
% bar(logical(activity(:,4)));
% ax(5)=subplot(96,1,5);
% bar(logical(activity(:,5)));
% ax(6)=subplot(96,1,6);
% bar(logical(activity(:,6)));
% ax(7)=subplot(96,1,7);
% bar(logical(activity(:,7)));
% ax(8)=subplot(96,1,8);
% bar(logical(activity(:,8)));
% ax(9)=subplot(96,1,9);
% bar(logical(activity(:,9)));
% ax(10)=subplot(96,1,10);
% bar(logical(activity(:,10)));
% ax(11)=subplot(96,1,11);
% bar(logical(activity(:,11)));
% ax(12)=subplot(96,1,12);
% bar(logical(activity(:,12)));
% ax(13)=subplot(96,1,13);
% bar(logical(activity(:,13)));
% ax(14)= subplot(96,1,14);
% bar(logical(activity(:,14)));
% ax(15)= subplot(96,1,15);
% bar(logical(activity(:,15)));
% ax(16)= subplot(96,1,16);
% bar(logical(activity(:,16)));
% ax(17)=subplot(96,1,17);
% bar(logical(activity(:,17)));
% ax(18)=subplot(96,1,18);
% bar(logical(activity(:,18)));
% ax(19)=subplot(96,1,19);
% bar(logical(activity(:,19)));
% ax(20)=subplot(96,1,20);
% bar(logical(activity(:,20)));
% ax(21)= subplot(96,1,21);
% bar(logical(activity(:,21)));
% ax(22)= subplot(96,1,22);
% bar(logical(activity(:,22)));
% ax(23)= subplot(96,1,23);
% bar(logical(activity(:,23)));
% ax(24)= subplot(96,1,24);
% bar(logical(activity(:,24)));
% ax(25)= subplot(96,1,25);
% bar(logical(activity(:,25)));
% ax(26)= subplot(96,1,26);
% bar(logical(activity(:,26)));
% ax(27)= subplot(96,1,27);
% bar(logical(activity(:,27)));
% ax(28)= subplot(96,1,28);
% bar(logical(activity(:,28)));
% ax(29)= subplot(96,1,29);
% bar(logical(activity(:,29)));
% ax(30)= subplot(96,1,30);
% bar(logical(activity(:,30)));
% ax(31)= subplot(96,1,31);
% bar(logical(activity(:,31)));
% ax(32)= subplot(96,1,32);
% bar(logical(activity(:,32)));
% ax(33)= subplot(96,1,33);
% bar(logical(activity(:,33)));
% ax(34)= subplot(96,1,34);
% bar(logical(activity(:,34)));
% ax(35)= subplot(96,1,35);
% bar(logical(activity(:,35)));
% ax(36)= subplot(96,1,36);
% bar(logical(activity(:,36)));
% ax(37)= subplot(96,1,37);
% bar(logical(activity(:,37)));
% ax(38)= subplot(96,1,38);
% bar(logical(activity(:,38)));
% ax(39)= subplot(96,1,39);
% bar(logical(activity(:,39)));
% ax(40)= subplot(96,1,40);
% bar(logical(activity(:,40)));
% ax(41)= subplot(96,1,41);
% bar(logical(activity(:,41)));
% ax(42)= subplot(96,1,42);
% bar(logical(activity(:,42)));
% ax(43)= subplot(96,1,43);
% bar(logical(activity(:,43)));
% ax(44)= subplot(96,1,44);
% bar(logical(activity(:,44)));
% ax(45)= subplot(96,1,45);
% bar(logical(activity(:,45)));
% ax(46)= subplot(96,1,46);
% bar(logical(activity(:,46)));
% ax(47)= subplot(96,1,47);
% bar(logical(activity(:,47)));
% ax(48)= subplot(96,1,48);
% bar(logical(activity(:,48)));
% ax(49)= subplot(96,1,49);
% bar(logical(activity(:,49)));
% ax(50)= subplot(96,1,50);
% bar(logical(activity(:,50)));
% ax(51)= subplot(96,1,51);
% bar(logical(activity(:,51)));
% ax(52)= subplot(96,1,54);
% bar(logical(activity(:,52)));
% ax(53)= subplot(96,1,53);
% bar(logical(activity(:,53)));
% ax(54)= subplot(96,1,54);
% bar(logical(activity(:,54)));
% ax(55)= subplot(96,1,55);
% bar(logical(activity(:,55)));
% ax(56)= subplot(96,1,56);
% bar(logical(activity(:,56)));
% ax(57)= subplot(96,1,57);
% bar(newtimenum,logical(activity(:,57)));
% 
% datetick('x','HHPM');
% % Turn off the X-tick labels in the top axes
% 
% set(ax(1:56),'XTickLabel','')
% 
% % Set the color of the X-axis in the top axes
% % to the axes background color
% 
% set(ax(1:56),'XColor',get(gca,'Color'))
% 
% % Turn off the box so that only the left 
% % vertical axis and bottom axis are drawn
% 
% set(ax,'box','off')

%% mean activity etc.
% allWells
format long;
clear meanAct xlabel ylabel title sumAct; 
meanAct = mean(a(:,allWells));  % mean of active seconds for each fish for whole experiment, single row vector
sumAct=sum(a(:,allWells)); % summed active seconds each fish for whole experiment, single row vect.

% Mean Fish activities overview
figure('Color',[1 1 1])  % make figure with white backround
subplot(2,3,1) % make figure with 6 panels
plot(mean(a(:,allWells),2)); % mean(a,2) gives mean of rows = mean for each time point

% Create xlabel
xlabel('Minutes','Interpreter','none');
% Create ylabel
ylabel('Active seconds/minute','Interpreter','none');
% Create title
title('1. Mean activity','Interpreter','none','FontSize',10);


% Mean Fish activities normal distributed?
subplot(2,3,2)
histfit(meanAct); % normal distribution of Mean activity change

% Create xlabel
xlabel('Mean seconds per experiment','Interpreter','none');
% Create ylabel
ylabel('Number of fish?','Interpreter','none');
% Create title
title('2. Distribution of mean activity','Interpreter','none','FontSize',10);

% Mean Fish activities normal probability?

subplot(2,3,3)
probplot('normal',meanAct); % probability distribution?


xlabel('Summed active seconds','Interpreter','none');
ylabel('Probability','Interpreter','none');
title('3. Probability of summed activity','Interpreter','none','FontSize',10);


% What are Frequencies of lazy and active fish? willie vs biene maya
clear nelements xcenters xvalues nbins
xvalues=0:0.2:4.5;
% nelements is number of fish per bin, (xcenter is center value of each
% bin? but has values of xvalues?)
[nelements,xcenters]=hist(meanAct,xvalues); 

subplot(2,3,4)
bar(xvalues, nelements);  % make bar plot for activity frequencies 

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[1 0 1],'EdgeColor','w');
% Create xlabel
xlabel('Active seconds per minute','Interpreter','none');
% Create ylabel
ylabel('Number of fish','Interpreter','none');
% Create title
title('4. Frequency of activity','Interpreter','none','FontSize',10);

% Probability for a fish to be lazy/willie or active/maja?

% subplot(2,3,4)
% probplot('normal',meanAct)
% % Create xlabel
% xlabel('Mean activity per fish','Interpreter','none');
% % Create ylabel
% ylabel('Probability','Interpreter','none','Fontsize',10);
% % Create title
% title('4. Probability of mean activity','Interpreter','none','FontSize',10);

% Ranked activity
clear sortmeanAct
[sortmeanAct]=sort(meanAct, 'descend');
subplot(2,3,5)
bar(sortmeanAct)

xlabel('Fish','Interpreter','none');
ylabel('Mean active seconds','Interpreter','none');
title('5. Sorted mean activity','Interpreter','none','FontSize',10);

% Cummulative activity

clear cumsumAct
cumsumAct=cumsum(a(:,allWells));
subplot(2,3,6)
plot(cumsumAct)

xlabel('Minutes','Interpreter','none');
ylabel('Cumulative active seconds','Interpreter','none');
title('6. Cumulative activity','Interpreter','none','FontSize',10);

fct_suptitle([' Overview activity allWells: ' expName]);


%% Group1

format long;
clear meanAct xlabel ylabel title sumAct; 
%meanAct = mean(a(:,Group1));  % mean of active seconds for each fish for whole experiment, single row vector
meanAct = mean(a(:,:)); %paula: what is group1? try all of a
%sumAct=sum(a(:,Group1)); % summed active seconds each fish for whole experiment, single row vect.
sumAct=sum(a(:,:)); %Paula: what is Group1? try all of a

% Mean Fish activities overview
figure('Color',[1 1 1])  % make figure with white backround
subplot(2,3,1) % make figure with 6 panels
%plot(mean(a(:,Group1),2)); % mean(a,2) gives mean of rows = mean for each time point
plot(mean(a(:,:),2)); %Paula: all of a, instead of only Group1
plot(mean(a(:,:),2)); %paula: all of a instead of until Group1

% Create xlabel
xlabel('Minutes','Interpreter','none');
% Create ylabel
ylabel('Active seconds/minute','Interpreter','none');
% Create title
title('1. Mean activity','Interpreter','none','FontSize',10);


% Mean Fish activities normal distributed?
subplot(2,3,2)
histfit(meanAct); % normal distribution of Mean activity change

% Create xlabel
xlabel('Mean seconds per experiment','Interpreter','none');
% Create ylabel
ylabel('Number of fish?','Interpreter','none');
% Create title
title('2. Distribution of mean activity','Interpreter','none','FontSize',10);

% Mean Fish activities normal probability?

subplot(2,3,3)
probplot('normal',meanAct); % probability distribution?


xlabel('Summed active seconds','Interpreter','none');
ylabel('Probability','Interpreter','none');
title('3. Probability of summed activity','Interpreter','none','FontSize',10);


% What are Frequencies of lazy and active fish? willie vs biene maya
clear nelements xcenters xvalues nbins
xvalues=0:0.2:4.5;
% nelements is number of fish per bin, (xcenter is center value of each
% bin? but has values of xvalues?)
[nelements,xcenters]=hist(meanAct,xvalues); 

subplot(2,3,4)
bar(xvalues, nelements);  % make bar plot for activity frequencies 

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[1 0 1],'EdgeColor','w');
% Create xlabel
xlabel('Active seconds per minute','Interpreter','none');
% Create ylabel
ylabel('Number of fish','Interpreter','none');
% Create title
title('4. Frequency of activity','Interpreter','none','FontSize',10);

% Probability for a fish to be lazy/willie or active/maja?

% subplot(2,3,4)
% probplot('normal',meanAct)
% % Create xlabel
% xlabel('Mean activity per fish','Interpreter','none');
% % Create ylabel
% ylabel('Probability','Interpreter','none','Fontsize',10);
% % Create title
% title('4. Probability of mean activity','Interpreter','none','FontSize',10);

% Ranked activity
clear sortmeanAct
[sortmeanAct]=sort(meanAct, 'descend');
subplot(2,3,5)
bar(sortmeanAct)

xlabel('Fish','Interpreter','none');
ylabel('Mean active seconds','Interpreter','none');
title('5. Sorted mean activity','Interpreter','none','FontSize',10);

% Cummulative activity

clear cumsumAct
%cumsumAct=cumsum(a(:,Group1));
cumsumAct=cumsum(a(:,:)); %Paula: all of a
subplot(2,3,6)
plot(cumsumAct)

xlabel('Minutes','Interpreter','none');
ylabel('Cumulative active seconds','Interpreter','none');
title('6. Cumulative activity','Interpreter','none','FontSize',10);

fct_suptitle([' Overview activity Group1: ' expName]);


%% Group2
%{
format long;
clear meanAct xlabel ylabel title sumAct; 
meanAct = mean(a(:,Group2));  % mean of active seconds for each fish for whole experiment, single row vector
sumAct=sum(a(:,Group2)); % summed active seconds each fish for whole experiment, single row vect.

% Mean Fish activities overview
figure('Color',[1 1 1])  % make figure with white backround
subplot(2,3,1) % make figure with 6 panels
plot(mean(a(:,Group2),2)); % mean(a,2) gives mean of rows = mean for each time point

% Create xlabel
xlabel('Minutes','Interpreter','none');
% Create ylabel
ylabel('Active seconds/minute','Interpreter','none');
% Create title
title('1. Mean activity','Interpreter','none','FontSize',10);


% Mean Fish activities normal distributed?
subplot(2,3,2)
histfit(meanAct); % normal distribution of Mean activity change

% Create xlabel
xlabel('Mean seconds per experiment','Interpreter','none');
% Create ylabel
ylabel('Number of fish?','Interpreter','none');
% Create title
title('2. Distribution of mean activity','Interpreter','none','FontSize',10);

% Mean Fish activities normal probability?

subplot(2,3,3)
probplot('normal',meanAct); % probability distribution?


xlabel('Summed active seconds','Interpreter','none');
ylabel('Probability','Interpreter','none');
title('3. Probability of summed activity','Interpreter','none','FontSize',10);


% What are Frequencies of lazy and active fish? willie vs biene maya
clear nelements xcenters xvalues nbins
xvalues=0:0.2:4.5;
% nelements is number of fish per bin, (xcenter is center value of each
% bin? but has values of xvalues?)
[nelements,xcenters]=hist(meanAct,xvalues); 

subplot(2,3,4)
bar(xvalues, nelements);  % make bar plot for activity frequencies 

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[1 0 1],'EdgeColor','w');
% Create xlabel
xlabel('Active seconds per minute','Interpreter','none');
% Create ylabel
ylabel('Number of fish','Interpreter','none');
% Create title
title('4. Frequency of activity','Interpreter','none','FontSize',10);

% Probability for a fish to be lazy/willie or active/maja?

% subplot(2,3,4)
% probplot('normal',meanAct)
% % Create xlabel
% xlabel('Mean activity per fish','Interpreter','none');
% % Create ylabel
% ylabel('Probability','Interpreter','none','Fontsize',10);
% % Create title
% title('4. Probability of mean activity','Interpreter','none','FontSize',10);

% Ranked activity
clear sortmeanAct
[sortmeanAct]=sort(meanAct, 'descend');
subplot(2,3,5)
bar(sortmeanAct)

xlabel('Fish','Interpreter','none');
ylabel('Mean active seconds','Interpreter','none');
title('5. Sorted mean activity','Interpreter','none','FontSize',10);

% Cummulative activity

clear cumsumAct
cumsumAct=cumsum(a(:,Group2));
subplot(2,3,6)
plot(cumsumAct)

xlabel('Minutes','Interpreter','none');
ylabel('Cumulative active seconds','Interpreter','none');
title('6. Cumulative activity','Interpreter','none','FontSize',10);

fct_suptitle([' Overview activity Group2: ' expName]);


end
%}
