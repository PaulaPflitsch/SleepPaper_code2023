%% 1. Experiment parameters.

clear b theshold NumFish FrameRate Interval expName

expName = 'LarvalSD_pilot_210121';



b = 1; % binlength (in minutes) to convert pixel data to seconds-with-motion (see below)
threshold = 2; % threshold how many pixels represents a real movement, rather than noise in the data?
NumFish = 96; % the number of fish in the box
%FrameRate = 30; % NEW NOTE: Now derive this from actual timestamp data below
Interval = 60; % the number of seconds over which we sum the bouts

%% 2. Load motion file (DELTA PIXELS)
%get filename and path, go to folder, open file
tic
[filename_motion, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_motion);             
clear data96
data96 = fread(fopen(filename_motion), inf, 'uint16=>uint16', 'b'); % returns one column with fish1 t1, fish2 t1, ...., fish96, fish1 t2, fish2 t2...
data96 = reshape(data96, 96, []);  %96 % number of rows m NumFish x number of columns n [] of data
data96 = rot90(data96, 3); % rotate matrix counterclockw by k*90 degree
data96 = fliplr(data96); % flip matrix left-right (vertical axis)
d96 = double(data96);
toc

%% 3. Load timestamp file
%get filename and path, go to folder, open file
tic
[filename_time, PathName, FilterIndex]=uigetfile('*');
cd(PathName);
fopen(filename_time);
clear timeArray           
timeArrayCell = importdata(filename_time); % import timestamp txt file to matlab
startTime = timeArrayCell(1,:);
endTime = timeArrayCell(end,:);
% determine length of time vector
[time_m, time_n] = size(timeArrayCell); 
toc

%% 4. Automatically find frame # for day/night transitions

% MANUALLY INPUT dates first according to timestamp file:
day0 = '5/3/2023'; % start date, no "0" in date!
day1 = '5/4/2023'; % Lights on at 9am (for ctrl); Run ends 11pm this day
%day2 = '21/4/2021'; % Run stopped this day
%day3 = '2/15/2019'; 

% output frame numbers for desired transitions 
% Below, specified a few phases based on your test file, but can be
% changed/renamed as desired. 

% beginning of first night 
frame_night1 = strfind(timeArrayCell, (strcat(day0,'11:30:00 PM')));
frame_night1 = find(not(cellfun('isempty', frame_night1)));
size_frame_night1 = size(frame_night1,1);
frame_night1 = frame_night1(1,1);

% beginning of light phase the next morning
frame_day1 = strfind(timeArrayCell, (strcat(day1,'9:00:00 AM')));
frame_day1 = find(not(cellfun('isempty', frame_day1)));
size_frame_day1 = size(frame_day1,1);
frame_day1 = frame_day1(1,1);

% beginning of lights off rebound test phase 
frame_testphase = strfind(timeArrayCell, (strcat(day1,'12:00:00 PM')));
frame_testphase = find(not(cellfun('isempty', frame_testphase)));
frame_testphase = frame_testphase(1,1);

% frame for 3pm, to split 9-12pm vs 12-3pm
%frame_3pm = strfind(timeArrayCell, (strcat(day1,'3:00:00 PM')));
%frame_3pm = find(not(cellfun('isempty', frame_3pm)));
%frame_3pm = frame_3pm(1,1);

% Time near the end of the experiment: for now use night 2
frame_night2 = strfind(timeArrayCell, (strcat(day1,'11:00:00 PM')));
frame_night2 = find(not(cellfun('isempty', frame_night2)));
size_frame_night2 = size(frame_night2,1);
frame_night2 = frame_night2(1,1);


% Mean framerate calculated based on desired time interval
% 210820 test dataset: looks like framerate is closer to 13-15??
frame_mean = (frame_day1-frame_night1)/(24*60*60);
% Set "official" framerate to nearest 0.5
framerate = round(frame_mean);
csvwrite('framerate.csv', framerate); %print this in case you need it later


%% 5. Clip data for consistent plotting using timestamp/frame:

% Using time interval from above: beginning of night until end of experiment
d96_clip = d96(frame_night1:frame_night2,:); % last one originally ":"
%d96_clip = d96(frame_testphase:frame_night2,:);
%% 6. Bin activity/rest
% According to Prober 2006 and Rihel 2010.
% Activity = seconds with motion, 
% Convert pixel data to seconds-with-motion data
% and bin into 1 minute intervals (because 'rest bout' as defined in these
% papers is a 1-minute period of continuous rest). Will process individual
% movement data in later section
clear b t a 
b = 1; % binlength in minutes to convert pixel data to seconds with motion
t = 2; % threshold how many pixels represent a real movement, rather than noise
%a = fct_BinActivity(d96,b,threshold,framerate);
a_clip = fct_BinActivity(d96_clip,b,threshold,framerate);
% activity is active seconds per minute
% For a standard 3day run, expect ~3480 bins.

%% 7. ENTER GENOTYPE INFORMATION
%!!!!!!!!!!!!!!!NOTE ORIENTATION OF NUMBER SYSTEM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Numbering here runs by column

% USE GUI OR ENTER MANUALLY:
%run geno2.m;

% SD experiments: Since you want to analyze only 1 group per plate, can
% use the 'ctrl' variable below as stand-in

% 210121_lightSD_ctrl_box9
%ctrl = [11  12  13  14  15  19  20  22  23  27  28  29  31  35  36  37  38  39  43  44  46  47  51  52  54  55  59  61  63  67  68  70  71  75  76  77  78  79  83  84  85  86  87];

% 210121_lightSD_3am_box2
%ctrl = [11  12  13  14  19  20  21  22  27  28  29  30  35  36  37  38  43  44  45  46  51  52  54  59  62  67  68  70  75  76  77  78  83  84  85  86];

% 210121_lightSD_6am_box1
%ctrl = [11  12  13  14  19  20  21  22  27  28  29  35  37  38  43  44  45  46  51  52  53  59  60  61  62  67  69  75  76  77  83  84  85];

% 210121_lightSD_continuous_box7
%ctrl = [11  12  13  14  19  22  30  35  36  37  38  45  46  51  52  53  54  59  61  62  68  70  75  76  78  83  84  85  86];

% 210121_lightSD_pulse_box8
%ctrl = [11  12  13  14  15  19  20  21  22  27  28  29  30  31  35  36  43  46  51  52  53  54  59  60  67  70  75  76  84  85  86];


% Rebound
% 210422 ctrl box 9
%ctrl = [10  11  12  13  14  18  19  20  21  22  26  27  29  30  34  35  36  37  38  43  44  45  46  51  52  53  54  59  60  61  62  63  67  68  69  70  71  75  76  77  78  79  83  84  85  86  87  94]; % Define fish #s here manually or with GUI above, etc

% 210422 LightSD 3am
%ctrl = [10  11  12  13  14  18  19  20  21  22  26  27  28  29  30  34  35  36  37  38  43  44  45  46  51  52  53  54  59  60  61  62  63  67  68  69  70  71  75  76  77  78  79  83  84  85  86  87];

% 210504 LightSD ctrl box 9
%ctrl = [10  11  12  13  14  15  18  19  20  21  22  23  26  27  28  29  30  31  34  35  36  37  38  39  42  43  44  45  46  47  50  51  53  54  55  58  59  60  61  62  63  66  67  68  69  70  71  74  75  76  77  78  79  82  83  84  85  87];

% 210504 LightSD pulse box 8
%ctrl = [10  11  12  13  14  15  18  19  20  21  22  23  26  27  28  29  30  31  34  35  36  37  38  39  42  43  44  45  46  47  50  51  53  54  55  58  59  60  61  62  63  66  67  68  69  70  71  74  75  76  77  78  79  82  83  84  85  86  87];

% 210504 SoundSD ctrl box10
%ctrl = [10  11  12  13  14  15  18  19  20  21  22  23  26  27  28  29  30  31  34  35  36  37  38  39  42  43  44  45  46  50  51  53  54  55  58  59  60  61  62  63  66  67  68  69 71  74  75  76  77  78  79  82  83  84  85  86  87];

% 20230503 LightSD control box7
%ctrl = [4  5  7  9  10  11  12  13  14  15  17  18  20  21  22  23  24  25  26 30  31  32  33  34  35  38  39  40  41  42  43  46  47  49  50  51  52  53  54  55   56  57  58  59  60  62  63  64  65  67  69  70  71  72  73  74  77  78  79  80  81  82  83  85   87  88  89  90  91  94  95  96]

% 20230503 LightSD continuous box8
ctrl = [2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18  19  21  22  23  24  25  26  27  29  30  31  32  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  72  73  74  75  77 78  79  80  81  82  83  84  89  90  92  93  94  95  96]


%230410 LightSD control box8
%ctrl = [20  21  22  27  28  30  35  36  37  38  43  44  45  46  51  52  53  54  61  62  67  68  69  70  75  76  77  78 ] % based on activity plot I removed: 19, 59

%230410 LightSD continuous box7
%ctrl = [19  20  21  22  27  28  29  30  35  37  38  43  44  45  51  52  53  54  59  60  61  62  67  68  69  70  78 ] % based on activity plot I removed: 46, 75, 76

%230419 LightSD control box7
%ctrl = [10  11  12  14  15  18  19  20  21  22  23  26  27  28  29  30  31  35  36  37  38  39  43  44  45  46  47 50  51  52  53  54  55  59  60  62  66  67  68  69  70  71  74  75  76  77  78  79  82  84  85  86  87]  

%230419 LightSD continuous box7
%ctrl = [10  11  14  18  19  20  21  22  23  26  27  28  29  30  34  35  36  38  39  42  43  44  46  47  50  51  52  53  54  55  58  59  60  61  62  67  68  69  70  71  74  76  78  79  82  83  84  85  86  87]

%230424 LightSD control box9 ONLY
%ctrl = [14  15  16  17  18  19  20  22  23  26  28  29  30  32  34  35  38  39  40  41  43  44  45  50  51  52  53  55  57  58  59  63  65  66  67  68  69  70  71  74  75  77  78  79  80  81  83]

%230501 LightSD control box7
%ctrl = [1  4  5  7  10  11  17 18  19  20  21  22  25  26  27  28  29  30  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  50  51  52  53  54  55  56  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79  80  81  82  85  86  87  88  90  93  94  95  96]

%230501 LightSD contiuous box8
%ctrl = [1  2  3  4  9  10  11  12  13  17  18  19  21  23  24  25  26  27  28  29  30  31  33  34  36  37  38  39  40  41  42  44  45  46  47  48  49  50  52  53  54  55  56  57  58  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  78  79  80  82  85  90  91  92  93  94  95]


%24 well plate LightSD rebound control
%ctrl = [1:24]

%24 well plate LightSD rebound continuous
%ctrl = [1:24]

% 211201 LightSD continuous box 8
%ctrl = [1:96]; %96
%ctrl(ctrl==4) = [];
%ctrl(ctrl==13) = [];
allWells = [1:96]; %96

% Save the fish #s in case you want to reanalyze later
fish_ids = fopen('fish_ids_group_1.txt','w');
fprintf(fish_ids,'%.0f,', ctrl);
fclose(fish_ids);

%% 8. OVERVIEW STATISTICS
clear subplot

geno = {'ctrl'; ctrl};

%Overall entire run info:
fct_crNormalDistr(a_clip, expName, geno, allWells, ctrl) %change to a for whole recording time
%fct_crNormalDistr(a, expName, geno, allWells, ctrl)

%% 9. ACTIVITY/REST PLOTS
% this calculates:
% waking activity (seconds of motion per 'waking minute')
%   where a 'waking minute' is a 
% minute with any activity
% rest (number of resting minutes in each 10 minute interval)
% activity (number of active seconds per 10 minute interval)

geno = {'ctrl'; ctrl};

%To plot the designated timeframe based on d96_clip:
% Save as svg for best results in Illustrator, etc
% MANUALLY SAVE from window to preserve display dimensions
a = a_clip
f = fct_crSleepPlots(a,geno,expName);
%f = fct_crSleepPlots(a,geno,expName);  %uncomment for unclipped excel
%files and plots

% Export structure f by genotype in case you want to replot outside Matlab
T_f = struct2table(f, 'AsArray', 1);
csvwrite('output_binned_rest_ctrl.csv', T_f.rest{1}); % ctrl
csvwrite('output_binned_wakeAct_ctrl.csv', T_f.wakeAct{1}); % ctrl
csvwrite('output_binned_act_ctrl.csv', T_f.activity{1}); % ctrl

%% 10. Split motion data into phases

d_night1 = d96(frame_night1:(frame_day1-1), :);
d_day1 = d96(frame_day1:(frame_night2-1), :);
d_testphase = d96(frame_testphase:(frame_night2-1), :);


%% 10a. Paula: Repeat all steps of d96 with d_night, d_day1 and d_testphase
% and save data for further analysis in python

% Bin activity & rest
clear b t a 
b = 1; % binlength in minutes to convert pixel data to seconds with motion
t = 2; % threshold how many pixels represent a real movement, rather than noise
a_night1 = fct_BinActivity(d_night1,b,threshold,framerate); % might have to change cut offs
a_day1 = fct_BinActivity(d_day1,b,threshold,framerate);
a_testphase = fct_BinActivity(d_testphase,b,threshold,framerate);

% use function fct_crSleepPlots to get desired parameters
% save data as csv
a = a_night1
f = fct_crSleepPlots(a,geno,expName);

T_f = struct2table(f, 'AsArray', 1);
csvwrite('output_binned_rest_night1.csv', T_f.rest{1}); 
csvwrite('output_binned_wakeAct_night1.csv', T_f.wakeAct{1}); 
csvwrite('output_binned_act_night1.csv', T_f.activity{1});


a = a_day1
f = fct_crSleepPlots(a,geno,expName);

T_f = struct2table(f, 'AsArray', 1);
csvwrite('output_binned_rest_day1.csv', T_f.rest{1}); 
csvwrite('output_binned_wakeAct_day1.csv', T_f.wakeAct{1}); 
csvwrite('output_binned_act_day1.csv', T_f.activity{1});


a = a_testphase
f = fct_crSleepPlots(a,geno,expName);

T_f = struct2table(f, 'AsArray', 1);
csvwrite('output_binned_rest_testphase.csv', T_f.rest{1}); 
csvwrite('output_binned_wakeAct_testphase.csv', T_f.wakeAct{1}); 
csvwrite('output_binned_act_testphase.csv', T_f.activity{1});



%% 11. BIN DAY/NIGHT ACTIVITY
%Active seconds per minute
a_night1 = fct_BinActivity(d_night1,b,threshold,framerate);
a_day1 = fct_BinActivity(d_day1,b,threshold,framerate);
a_testphase = fct_BinActivity(d_testphase,b,threshold,framerate);


%% 12. REST LATENCY + Plot
% find indices of first zeros in each column
% Find first zero: This works because data in 'a' is already binned per
% minute (depending on your settings);

% night1
[logi, la]= max(a_night1 == 0, [], 1);   % stackoverflow, c gives index maxima, if several identical returns first
lat_night1 = [la(geno{2,1})]; % extract values for genotypes

% day1
[logi la]= max(a_day1 == 0, [], 1);   % stackoverflow, c gives index maxima, if several identical returns first
lat_day1 = [la(geno{2,1})];

% testphase
[logi la]= max(a_testphase == 0, [], 1);   % stackoverflow, c gives index maxima, if several identical returns first
lat_testphase = [la(geno{2,1})];


% combine into user-friendly table
lat = catpad(2,lat_night1.', lat_day1.', lat_testphase.');
T_lat = table(lat(:,1), lat(:,2), lat(:,3),...
        'VariableNames',{'night1','day1', 'testphase'});

% Export latency data:
writetable(T_lat, 'output_latency.csv');

% Plot average rest latency
% figure
% subplot(1,3,1)
% CategoricalScatterplot(T_lat{:,{'ctrl_n1', 'hom_n1', 'het_n1'}}, 'Labels',{'ctrl_n1', 'hom_n1', 'het_n1'});
% title('Night1')
% ylabel('Latency (min)')
% ylim([0,inf])
% set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
% subplot(1,3,2)
% CategoricalScatterplot(T_lat{:,{'ctrl_n2', 'hom_n2', 'het_n2'}}, 'Labels',{'ctrl_n2', 'hom_n2', 'het_n2'});
% title('Night2')
% ylim([0,inf])
% set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
% subplot(1,3,3)
% CategoricalScatterplot(T_lat{:,{'ctrl_n3', 'hom_n3', 'het_n3'}}, 'Labels',{'ctrl_n3', 'hom_n3', 'het_n3'});
% title('Night3')
% ylim([0,inf])
% set(gca,'FontSize', 12, 'XTickLabelRotation',45, 'TickLabelInterpreter','none')
% saveas(gcf,'mean_RestLatency_night', 'svg'); % save as svg
% saveas(gcf,'mean_RestLatency_night', 'png'); % save as png


figure
CategoricalScatterplot(T_lat{:,{'night1', 'day1', 'testphase'}}, 'Labels',{'night1', 'day1', 'testphase'});
title('Day1')
ylabel('Latency (min)')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
%saveas(gcf,'mean_RestLatency_day', 'svg'); % save as svg
saveas(gcf,'mean_RestLatency_day', 'png'); % save as png

%% 13. AVERAGE REST BOUT NUMBER AND LENGTH
%First, split a by time of day and genotype.
%Then use fct_SearchRuns to calculate RestBout number and lengths 
%based on runs of 0s in binned activity data 'a'.

%Night1
a_night1_ctrl = a_night1(:,ctrl);
[night1_ctrl_RestBoutNumber, night1_ctrl_RestBoutLength] = fct_SearchRuns(a_night1_ctrl);

%Testphase
a_testphase_ctrl = a_testphase(:,ctrl);
[testphase_ctrl_RestBoutNumber, testphase_ctrl_RestBoutLength] = fct_SearchRuns(a_testphase_ctrl);

%Day1
a_day1_ctrl = a_day1(:,ctrl);
[day1_ctrl_RestBoutNumber, day1_ctrl_RestBoutLength] = fct_SearchRuns(a_day1_ctrl);

% Also separate overall 'a_clip' by genotype and extract overall value
% across entire time:
a_clip_ctrl = a_clip(:,ctrl);
[total_ctrl_RestBoutNumber, total_ctrl_RestBoutLength] = fct_SearchRuns(a_clip_ctrl);


% combine into user-friendly tables:
RestBoutLength_ctrl = catpad(2,night1_ctrl_RestBoutLength.', day1_ctrl_RestBoutLength.', testphase_ctrl_RestBoutLength.');
RestBoutNumber_ctrl = catpad(2,night1_ctrl_RestBoutNumber.', day1_ctrl_RestBoutNumber.', testphase_ctrl_RestBoutNumber.');

T_RestBoutLength = table(RestBoutLength_ctrl(:,1), RestBoutLength_ctrl(:,2), RestBoutLength_ctrl(:,3),...
       'VariableNames',{'Night1','Day1', 'Testphase'});
   
T_RestBoutNumber = table(RestBoutNumber_ctrl(:,1), RestBoutNumber_ctrl(:,2), RestBoutNumber_ctrl(:,3),...
       'VariableNames',{'Night1','Day1', 'Testphase'});
          
% Export rest bout data:
writetable(T_RestBoutLength, 'output_RestBoutLength.csv');
writetable(T_RestBoutNumber, 'output_RestBoutNumber.csv');



%% 14. Plot average rest bout number/length
% rest bout number
figure
CategoricalScatterplot(T_RestBoutNumber{:,{'Night1', 'Day1', 'Testphase'}}, 'Labels',{'Night1', 'Day1', 'Testphase'});
ylabel('# Rest Bouts')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_RestBoutNumber', 'png'); % save as png

% rest bout length
figure
CategoricalScatterplot(T_RestBoutLength{:,{'Night1', 'Day1', 'Testphase'}}, 'Labels',{'Night1', 'Day1', 'Testphase'});
ylabel('Rest Bout Length (min)')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_RestBoutLength_night', 'png'); % save as png



%% 15. Plot cumulative rest over time for defined time windows
%Same as plotting all the 0s in a_rebound, for example. 
%Turn a_rebound_hom into logical values. 
%a_rebound_hom_logical = logical(a_rebound_hom);
%a_rebound_hom_logical_rest = ~a_rebound_hom_logical;
%a_rebound_hom_logical_rest = double(a_rebound_hom_logical_rest);

%cumulative_rest_rebound = cumsum(a_rebound_hom_logical_rest);
%cumulative_rest_rebound_flip = double(rot90(cumulative_rest_rebound));
%plot_areaerrorbar(cumulative_rest_rebound_flip);

%Overview as specified by a_clip time interval above
a_clip_ctrl_logical = logical(a_clip_ctrl);
a_clip_ctrl_logical_rest = ~a_clip_ctrl_logical;
a_clip_ctrl_logical_rest = double(a_clip_ctrl_logical_rest);
cumulative_rest_ctrl = cumsum(a_clip_ctrl_logical_rest);
plot_areaerrorbar(double(rot90(cumulative_rest_ctrl)));
%%%%%%%%%%%%%%%%%%%%
title('Overview')
xlabel('Minutes')
ylabel('Cumulative Rest (min)')
ylim([0,inf])
set(gca,'FontSize', 12)
saveas(gcf,'cumulative_rest_overview', 'png'); % save as png
hold off
close gcf
%%%%%%%%%%%%%%%%%%%%

%night1
a_night1_ctrl_logical = logical(a_night1_ctrl);
a_night1_ctrl_logical_rest = ~a_night1_ctrl_logical;
a_night1_ctrl_logical_rest = double(a_night1_ctrl_logical_rest);
cumulative_rest_night1_ctrl = cumsum(a_night1_ctrl_logical_rest);
plot_areaerrorbar(double(rot90(cumulative_rest_night1_ctrl)));
title('Night1')
xlabel('Minutes')
ylabel('Cumulative Rest (min)')
ylim([0,inf])
set(gca,'FontSize', 12)
saveas(gcf,'cumulative_rest_night1', 'png'); % save as png

% day1
a_testphase_ctrl_logical = logical(a_day1_ctrl);
a_day1_testphase_logical_rest = ~a_testphase_ctrl_logical;
a_day1_testphase_logical_rest = double(a_day1_testphase_logical_rest);
cumulative_rest_day1_ctrl = cumsum(a_day1_testphase_logical_rest);
plot_areaerrorbar(double(rot90(cumulative_rest_day1_ctrl)));
title('Day1')
xlabel('Minutes')
ylabel('Cumulative Rest (min)')
ylim([0,inf])
set(gca,'FontSize', 12)
saveas(gcf,'cumulative_rest_day1', 'png'); % save as png


% testphase
a_testphase_ctrl_logical = logical(a_testphase_ctrl);
a_day1_testphase_logical_rest = ~a_testphase_ctrl_logical;
a_day1_testphase_logical_rest = double(a_day1_testphase_logical_rest);
cumulative_rest_testphase_ctrl = cumsum(a_day1_testphase_logical_rest);
plot_areaerrorbar(double(rot90(cumulative_rest_testphase_ctrl)));
title('Test phase')
xlabel('Minutes')
ylabel('Cumulative Rest (min)')
ylim([0,inf])
set(gca,'FontSize', 12)
saveas(gcf,'cumulative_rest_testphase', 'png'); % save as png



% Export all cumulative rest bout data:
csvwrite('output_cumulative_rest_night1.csv', cumulative_rest_night_ctrl);
csvwrite('output_cumulative_rest_day1.csv', cumulative_rest_day_ctrl);
csvwrite('output_cumulative_rest_testphase.csv', cumulative_rest_testphase_ctrl);

%% Average activity bout number/length

%First, re-bin activity by SECOND rather than MINUTE, if defining one
%movement bout as activity separated by a second? This is only useful to
%get ACTIVITY BOUT length/number, rather than individual movements. 

b_sec = (1/60); % binlength in minutes to convert pixel data to seconds with motion
threshold = 2; % threshold how many pixels represent a real movement, rather than noise

a_sec_clip = fct_BinActivity(d96_clip,b_sec,threshold,framerate);

a_sec_day1 = fct_BinActivity(d_day1,b_sec,threshold,framerate);
a_sec_night1 = fct_BinActivity(d_night1,b_sec,threshold,framerate);
a_sec_testphase = fct_BinActivity(d_testphase,b_sec,threshold,framerate);

%Then use fct_SearchRuns_activity to calculate bout number and lengths 
%based on runs of ones in binned activity data 'a'.

%Night1
a_sec_night1_ctrl = a_sec_night1(:,ctrl);
[night1_ctrl_ActivityBoutNumber, night1_ctrl_ActivityBoutLength] = fct_SearchRuns_activity(a_sec_night1_ctrl);

%Day1
a_sec_day1_ctrl = a_sec_day1(:,ctrl);
[day1_ctrl_ActivityBoutNumber, day1_ctrl_ActivityBoutLength] = fct_SearchRuns_activity(a_sec_day1_ctrl);

%Testphase
a_sec_testphase_ctrl = a_sec_testphase(:,ctrl);
[testphase_ctrl_ActivityBoutNumber, testphase_ctrl_ActivityBoutLength] = fct_SearchRuns_activity(a_sec_testphase_ctrl);

%Total
a_sec_clip_ctrl = a_sec_clip(:,ctrl);
[total_ctrl_ActivityBoutNumber, total_ctrl_ActivityBoutLength] = fct_SearchRuns_activity(a_sec_clip_ctrl);

% Save user-friendly tables:
night_ActivityBoutNumber_ctrl = night1_ctrl_ActivityBoutNumber.';
day_ActivityBoutNumber_ctrl = day1_ctrl_ActivityBoutNumber.';
testphase_ActivityBoutNumber_ctrl = testphase_ctrl_ActivityBoutNumber.';

night_ActivityBoutLength_ctrl = night1_ctrl_ActivityBoutLength.';
day_ActivityBoutLength_ctrl = day1_ctrl_ActivityBoutLength.';
testphase_ActivityBoutLength_ctrl = testphase_ctrl_ActivityBoutLength.';


% activity bout number
ActivityBoutNumber = catpad(2,night_ActivityBoutNumber_ctrl,...
                                day_ActivityBoutNumber_ctrl,...
                                testphase_ActivityBoutNumber_ctrl);
% convert to table
T_ActivityBoutNumber = table(ActivityBoutNumber(:,1),...
    ActivityBoutNumber(:,2), ActivityBoutNumber(:,3),...
    'VariableNames',{'night','day', 'testphase'});


% activity bout length
ActivityBoutLength = catpad(2,night_ActivityBoutLength_ctrl,...
                                day_ActivityBoutLength_ctrl,...
                                testphase_ActivityBoutLength_ctrl);
% convert to table
T_ActivityBoutLength = table(ActivityBoutLength(:,1),...
    ActivityBoutLength(:,2), ActivityBoutLength(:,3),...
    'VariableNames',{'night','day', 'testphase'});
     
            
% total overview               
total_ActivityBoutNumber = total_ctrl_ActivityBoutNumber.';
T_ActivityBoutNumber_total = table(total_ActivityBoutNumber(:,1),...
                               'VariableNames', {'total'});
total_ActivityBoutLength = total_ctrl_ActivityBoutLength.';
T_ActivityBoutLength_total = table(total_ActivityBoutLength(:,1),...
                               'VariableNames', {'total'});

disp('Saving output data...')

% Export all activity bout data:
writetable(T_ActivityBoutNumber, 'output_ActivityBoutNumber.csv');
writetable(T_ActivityBoutNumber_total, 'output_ActivityBoutNumber_total.csv');
writetable(T_ActivityBoutLength, 'output_ActivityBoutLength.csv');
writetable(T_ActivityBoutLength_total, 'output_ActivityBoutLength_total.csv');

disp('Average activity bout section complete')

%% Activity bout mean plots

% Activity bout number
figure
subplot(1,3,1)
CategoricalScatterplot(T_ActivityBoutNumber_night{:,{'ctrl_n1', 'hom_n1', 'het_n1'}}, 'Labels',{'ctrl_n1', 'hom_n1', 'het_n1'});
title('Night1')
ylabel('# Activity Bouts')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,2)
CategoricalScatterplot(T_ActivityBoutNumber_night{:,{'ctrl_n2', 'hom_n2', 'het_n2'}}, 'Labels',{'ctrl_n2', 'hom_n2', 'het_n2'});
title('Night2')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,3)
CategoricalScatterplot(T_ActivityBoutNumber_night{:,{'ctrl_n3', 'hom_n3', 'het_n3'}}, 'Labels',{'ctrl_n3', 'hom_n3', 'het_n3'});
title('Night3')
ylim([0,inf])
set(gca,'FontSize', 12, 'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_ActivityBoutNumber_night', 'svg'); % save as svg
saveas(gcf,'mean_ActivityBoutNumber_night', 'png'); % save as png

figure
subplot(1,2,1)
CategoricalScatterplot(T_ActivityBoutNumber_day{:,{'ctrl_d1', 'hom_d1', 'het_d1'}}, 'Labels',{'ctrl_d1', 'hom_d1', 'het_d1'});
title('Day1')
ylabel('# Activity Bouts')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,2,2)
CategoricalScatterplot(T_ActivityBoutNumber_day{:,{'ctrl_d2', 'hom_d2', 'het_d2'}}, 'Labels',{'ctrl_d2', 'hom_d2', 'het_d2'});
title('Day2')
ylim([0,inf])
set(gca,'FontSize', 12, 'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_ActivityBoutNumber_day', 'svg'); % save as svg
saveas(gcf,'mean_ActivityBoutNumber_day', 'png'); % save as png

figure
CategoricalScatterplot(T_ActivityBoutNumber_total{:,{'ctrl', 'hom', 'het'}}, 'Labels',{'ctrl', 'hom', 'het'});
title('Tota1')
ylabel('# Activity Bouts')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_ActivityBoutNumber_total', 'svg'); % save as svg
saveas(gcf,'mean_ActivityBoutNumber_total', 'png'); % save as png

% Activity bout Length
figure
subplot(1,3,1)
CategoricalScatterplot(T_ActivityBoutLength_night{:,{'ctrl_n1', 'hom_n1', 'het_n1'}}, 'Labels',{'ctrl_n1', 'hom_n1', 'het_n1'});
title('Night1')
ylabel('Activity Bout Length (sec)')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,2)
CategoricalScatterplot(T_ActivityBoutLength_night{:,{'ctrl_n2', 'hom_n2', 'het_n2'}}, 'Labels',{'ctrl_n2', 'hom_n2', 'het_n2'});
title('Night2')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,3)
CategoricalScatterplot(T_ActivityBoutLength_night{:,{'ctrl_n3', 'hom_n3', 'het_n3'}}, 'Labels',{'ctrl_n3', 'hom_n3', 'het_n3'});
title('Night3')
ylim([0,inf])
set(gca,'FontSize', 12, 'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_ActivityBoutLength_night', 'svg'); % save as svg
saveas(gcf,'mean_ActivityBoutLength_night', 'png'); % save as png

figure
subplot(1,2,1)
CategoricalScatterplot(T_ActivityBoutLength_day{:,{'ctrl_d1', 'hom_d1', 'het_d1'}}, 'Labels',{'ctrl_d1', 'hom_d1', 'het_d1'});
title('Day1')
ylabel('Activity Bout Length (sec)')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,2,2)
CategoricalScatterplot(T_ActivityBoutLength_day{:,{'ctrl_d2', 'hom_d2', 'het_d2'}}, 'Labels',{'ctrl_d2', 'hom_d2', 'het_d2'});
title('Day2')
ylim([0,inf])
set(gca,'FontSize', 12, 'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_ActivityBoutLength_day', 'svg'); % save as svg
saveas(gcf,'mean_ActivityBoutLength_day', 'png'); % save as png

figure
CategoricalScatterplot(T_ActivityBoutLength_total{:,{'ctrl', 'hom', 'het'}}, 'Labels',{'ctrl', 'hom', 'het'});
title('Total')
ylabel('Activity Bout Length (sec)')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_ActivityBoutLength_total', 'svg'); % save as svg
saveas(gcf,'mean_ActivityBoutLength_total', 'png'); % save as png

%% Cumulative plots of active sec/min

% Overview
a_sec_clip_ctrl_logical = logical(a_sec_clip_ctrl);
a_clip_ctrl_logical = double(a_sec_clip_ctrl_logical);
cumulative_activity_ctrl = cumsum(a_clip_ctrl_logical);
plot_areaerrorbar(double(rot90(cumulative_activity_ctrl)));
%%%%%%%%%%%%%%%%%%%%
title('Overview')
xlabel('Seconds')
ylabel('Cumulative Activity (sec)')
ylim([0,inf])
set(gca,'FontSize', 12)
saveas(gcf,'cumulative_activity_overview', 'png'); % save as png
hold off
close gcf

% Night
a_sec_night1_ctrl_logical = logical(a_sec_night1_ctrl);
a_night1_ctrl_logical = double(a_sec_night1_ctrl_logical);
cumulative_activity_night1_ctrl = cumsum(a_night1_ctrl_logical);
plot_areaerrorbar(double(rot90(cumulative_activity_night1_ctrl)));
%%%%%%%%%%%%%%%%%%%%
title('Night')
xlabel('Seconds')
ylabel('Cumulative Activity (sec)')
ylim([0,inf])
set(gca,'FontSize', 12)
saveas(gcf,'cumulative_activity_night2', 'png'); % save as png
hold off
close gcf

% day1
a_sec_day1_ctrl_logical = logical(a_sec_day1_ctrl);
a_testphase_ctrl_logical = double(a_sec_day1_ctrl_logical);
cumulative_activity_day1_ctrl = cumsum(a_testphase_ctrl_logical);
plot_areaerrorbar(double(rot90(cumulative_activity_day1_ctrl)));
title('Day')
xlabel('Seconds')
ylabel('Cumulative Activity (sec)')
ylim([0,inf])
set(gca,'FontSize', 12)
saveas(gcf,'cumulative_activity_day1', 'png'); % save as png
hold off
close gcf

% testphase
a_sec_testphase_ctrl_logical = logical(a_sec_testphase_ctrl);
a_testphase_ctrl_logical = double(a_sec_testphase_ctrl_logical);
cumulative_activity_testphase_ctrl = cumsum(a_testphase_ctrl_logical);
plot_areaerrorbar(double(rot90(cumulative_activity_testphase_ctrl)));
%%%%%%%%%%%%%%%%%%%%
title('Testphase')
xlabel('Seconds')
ylabel('Cumulative Activity (sec)')
ylim([0,inf])
set(gca,'FontSize', 12)
saveas(gcf,'cumulative_activity_day2', 'png'); % save as png
hold off
close gcf

%% Export cumulative activity bout (sec/min) data
csvwrite('output_cumulative_ActivityBout_total.csv', cumulative_activity_ctrl);
csvwrite('output_cumulative_ActivityBout_night.csv', cumulative_activity_night1_ctrl);
csvwrite('output_cumulative_ActivityBout_day.csv', cumulative_activity_day1_ctrl);
csvwrite('output_cumulative_ActivityBout_testphase.csv', cumulative_activity_testphase_ctrl);

disp('Done exporting cumulative activity bout data')

%% Analyze individual movements
%remove motion below a certain threshold of delta-pixels
tic
disp('Filtering movements...')
m = d96_clip;
m(m<=threshold)=0; % use threshold to remove background
m_night1 = d_night1;
m_day1 = d_day1;
m_testphase = d_testphase;

m_night1(m_night1<=threshold)=0;
m_day1(m_day1<=threshold)=0;
m_testphase(m_testphase<=threshold)=0;

% Separate movement data by genotype
m_night1 = m_night1(:,ctrl);
m_night1_cell = num2cell(m_night1, 1);
m_day1 = m_day1(:,ctrl);
m_day1_cell = num2cell(m_day1, 1);
m_testphase = m_testphase(:,ctrl);
m_testphase_cell = num2cell(m_testphase, 1);

% Transpose to prepare for movement extraction function
m_night1_cell = cellfun(@transpose, m_night1_cell, 'UniformOutput', 0);
m_day1_cell = cellfun(@transpose, m_day1_cell, 'UniformOutput', 0);
m_testphase_cell = cellfun(@transpose, m_testphase_cell, 'UniformOutput', 0);


% Extract profiles of individual movements:
disp('Extracting individual movements...')
clear r
mvmnts_night1 = cellfun(@(r) regionprops(logical(r),r, 'PixelValues'), m_night1_cell, 'UniformOutput', 0);
mvmnts_day1 = cellfun(@(r) regionprops(logical(r),r, 'PixelValues'), m_day1_cell, 'UniformOutput', 0);
mvmnts_testphase = cellfun(@(r) regionprops(logical(r),r, 'PixelValues'), m_testphase_cell, 'UniformOutput', 0);

% Convert each array of structures into array of cells:
mvmnts_night1 = cellfun(@struct2cell, mvmnts_night1,'UniformOutput', 0);
mvmnts_night1 = cellfun(@transpose, mvmnts_night1, 'UniformOutput', 0);
mvmnts_day1 = cellfun(@struct2cell, mvmnts_day1,'UniformOutput', 0);
mvmnts_day1 = cellfun(@transpose, mvmnts_day1, 'UniformOutput', 0);
mvmnts_testphase = cellfun(@struct2cell, mvmnts_testphase,'UniformOutput', 0);
mvmnts_testphase = cellfun(@transpose, mvmnts_testphase, 'UniformOutput', 0);

% Sum up each run; Output is amplitude of each movement for each fish
mvmnt_amp_night1 = cell(size(mvmnts_night1)); %pre-allocate
for k = 1:length(mvmnts_night1)
	mvmnt_amp_night1{k} = cellfun(@sum, mvmnts_night1{k}, 'UniformOutput', 0);
end

mvmnt_amp_day1 = cell(size(mvmnts_day1)); %pre-allocate
for k = 1:length(mvmnts_day1)
	mvmnt_amp_day1{k} = cellfun(@sum, mvmnts_day1{k}, 'UniformOutput', 0);
end

mvmnt_amp_testphase = cell(size(mvmnts_testphase)); %pre-allocate
for k = 1:length(mvmnts_testphase)
	mvmnt_amp_testphase{k} = cellfun(@sum, mvmnts_testphase{k}, 'UniformOutput', 0);
end

disp('Gathering movement amplitudes...')
% concat all the mvmnts/fish into one array each:
mvmnt_amp_night1 = cellfun(@cell2mat, mvmnt_amp_night1, 'UniformOutput', 0);
mvmnt_amp_day1 = cellfun(@cell2mat, mvmnt_amp_day1, 'UniformOutput', 0);
mvmnt_amp_testphase = cellfun(@cell2mat, mvmnt_amp_testphase, 'UniformOutput', 0);

toc
disp('Extract movements step complete.')

%% Get # and amplitude of all movements <10pix and >10pix separately:
disp('Classifying movements...')
tic
%night1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% small
indices_mvmnt_small_night1 = cellfun(@(x) x<10, mvmnt_amp_night1, 'UniformOutput', 0); %indices
mvmnt_amp_night1_small = cell(size(mvmnt_amp_night1)); 
for k = 1:length(mvmnt_amp_night1)
    mvmnt_amp_night1_small{k} = mvmnt_amp_night1{k}(indices_mvmnt_small_night1{k});
end
% big
indices_mvmnt_big_night1 = cellfun(@(x) x>10, mvmnt_amp_night1, 'UniformOutput', 0); %indices
mvmnt_amp_night1_big = cell(size(mvmnt_amp_night1));
for k = 1:length(mvmnt_amp_night1)
    mvmnt_amp_night1_big{k} = mvmnt_amp_night1{k}(indices_mvmnt_big_night1{k});
end

%day1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% small
indices_mvmnt_small_day1 = cellfun(@(x) x<10, mvmnt_amp_day1, 'UniformOutput', 0); %indices
mvmnt_amp_day1_small = cell(size(mvmnt_amp_day1)); 
for k = 1:length(mvmnt_amp_day1)
    mvmnt_amp_day1_small{k} = mvmnt_amp_day1{k}(indices_mvmnt_small_day1{k});
end
% big
indices_mvmnt_big_day1 = cellfun(@(x) x>10, mvmnt_amp_day1, 'UniformOutput', 0); %indices
mvmnt_amp_day1_big = cell(size(mvmnt_amp_day1));
for k = 1:length(mvmnt_amp_day1)
    mvmnt_amp_day1_big{k} = mvmnt_amp_day1{k}(indices_mvmnt_big_day1{k});
end

%testphase%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% small
indices_mvmnt_small_testphase = cellfun(@(x) x<10, mvmnt_amp_testphase, 'UniformOutput', 0); %indices
mvmnt_amp_testphase_small = cell(size(mvmnt_amp_testphase)); 
for k = 1:length(mvmnt_amp_testphase)
    mvmnt_amp_testphase_small{k} = mvmnt_amp_testphase{k}(indices_mvmnt_small_testphase{k});
end
% big
indices_mvmnt_big_testphase = cellfun(@(x) x>10, mvmnt_amp_testphase, 'UniformOutput', 0); %indices
mvmnt_amp_testphase_big = cell(size(mvmnt_amp_testphase));
for k = 1:length(mvmnt_amp_testphase)
    mvmnt_amp_testphase_big{k} = mvmnt_amp_testphase{k}(indices_mvmnt_big_testphase{k});
end

toc
disp('Done classifying movements.')


%% How long is each type of movement? Use small/big indices to count frames
disp('Quantifying movement lengths...')
tic
%night1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% small
mvmnt_night1_small = cell(size(mvmnts_night1));
for k = 1:length(mvmnts_night1)
    mvmnt_night1_small{k} = mvmnts_night1{k}(indices_mvmnt_small_night1{k});
    mvmnt_len_night1_small{k} = cellfun(@length, mvmnt_night1_small{k});
end
% big
mvmnt_night1_big = cell(size(mvmnts_night1));
for k = 1:length(mvmnts_night1)
    mvmnt_night1_big{k} = mvmnts_night1{k}(indices_mvmnt_big_night1{k});
    mvmnt_len_night1_big{k} = cellfun(@length, mvmnt_night1_big{k});
end

%day1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% small
mvmnt_day1_small = cell(size(mvmnts_day1));
for k = 1:length(mvmnts_day1)
    mvmnt_day1_small{k} = mvmnts_day1{k}(indices_mvmnt_small_day1{k});
    mvmnt_len_day1_small{k} = cellfun(@length, mvmnt_day1_small{k});
end
% big
mvmnt_day1_big = cell(size(mvmnts_day1));
for k = 1:length(mvmnts_day1)
    mvmnt_day1_big{k} = mvmnts_day1{k}(indices_mvmnt_big_day1{k});
    mvmnt_len_day1_big{k} = cellfun(@length, mvmnt_day1_big{k});
end

%testphase%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% small
mvmnt_testphase_small = cell(size(mvmnts_testphase));
for k = 1:length(mvmnts_testphase)
    mvmnt_testphase_small{k} = mvmnts_testphase{k}(indices_mvmnt_small_testphase{k});
    mvmnt_len_testphase_small{k} = cellfun(@length, mvmnt_testphase_small{k});
end
% big
mvmnt_testphase_big = cell(size(mvmnts_testphase));
for k = 1:length(mvmnts_testphase)
    mvmnt_testphase_big{k} = mvmnts_testphase{k}(indices_mvmnt_big_testphase{k});
    mvmnt_len_testphase_big{k} = cellfun(@length, mvmnt_testphase_big{k});
end

toc
disp('Done quantifying movement lengths')


%% Intermovement_Intervals_Section1
% Extract all intermovement zeros and reshape data
disp('Extracting intermovement intervals...')
tic 
clear r
imi_night1 = cellfun(@(r) regionprops(~logical(r),r, 'PixelValues'), m_night1_cell, 'UniformOutput', 0);
imi_day1 = cellfun(@(r) regionprops(~logical(r),r, 'PixelValues'), m_day1_cell, 'UniformOutput', 0);
imi_testphase = cellfun(@(r) regionprops(~logical(r),r, 'PixelValues'), m_testphase_cell, 'UniformOutput', 0);

% Convert each array of structures into array of cells:
imi_night1 = cellfun(@struct2cell, imi_night1,'UniformOutput', 0);
imi_day1 = cellfun(@struct2cell, imi_day1,'UniformOutput', 0);
imi_testphase = cellfun(@struct2cell, imi_testphase,'UniformOutput', 0);

% Transpose
imi_night1 = cellfun(@transpose, imi_night1, 'UniformOutput', 0);
imi_day1 = cellfun(@transpose, imi_day1, 'UniformOutput', 0);
imi_testphase = cellfun(@transpose, imi_testphase, 'UniformOutput', 0);

toc
disp('Intermovement_Intervals_Section1 complete')


%% Intermovement_Intervals_Section2
tic
% Get the lengths of all the intermovement intervals:
imi_len_night1 = cell(size(imi_night1)); % preallocate
for k = 1:length(imi_night1)
    imi_night1{k} = imi_night1{k};
    imi_len_night1{k} = cellfun(@length, imi_night1{k});
end

imi_len_day1 = cell(size(imi_day1)); % preallocate
for k = 1:length(imi_day1)
    imi_day1{k} = imi_day1{k};
    imi_len_day1{k} = cellfun(@length, imi_day1{k});
end

imi_len_testphase = cell(size(imi_testphase)); % preallocate
for k = 1:length(imi_testphase)
    imi_testphase{k} = imi_testphase{k};
    imi_len_testphase{k} = cellfun(@length, imi_testphase{k});
end

% OUTSTANDING ISSUES:
% Do these estimates of rest diverge from binned data above?
%test_rest_indices = imi_len_day1_ctrl{2}>(framerate*60);
%test_rest = imi_len_day1_ctrl{2}(test_rest_indices);

toc
disp('Intermovement_Intervals_Section2 complete')


%% Compile movement data and means for plotting/export

% Mean amplitude of movements/fish
mvmnt_amp_night1_small_mean = cellfun(@mean, mvmnt_amp_night1_small);
mvmnt_amp_night1_big_mean = cellfun(@mean, mvmnt_amp_night1_big);
mvmnt_amp_day1_small_mean = cellfun(@mean, mvmnt_amp_day1_small);
mvmnt_amp_day1_big_mean = cellfun(@mean, mvmnt_amp_day1_big);
mvmnt_amp_testphase_small_mean = cellfun(@mean, mvmnt_amp_testphase_small);
mvmnt_amp_testphase_big_mean = cellfun(@mean, mvmnt_amp_testphase_big);

% Compile amplitude data into tables
mean_mvmnt_amp_small = catpad(2, mvmnt_amp_night1_small_mean.',...
                                    mvmnt_amp_day1_small_mean.',...
                                    mvmnt_amp_testphase_small_mean.');

T_mean_mvmnt_amp_small = table(mean_mvmnt_amp_small(:,1),...
                                mean_mvmnt_amp_small(:,2),...
                                mean_mvmnt_amp_small(:,3),...
    'VariableNames',{'night','day', 'testphase'});


mean_mvmnt_amp_big = catpad(2, mvmnt_amp_night1_big_mean.',...
                                    mvmnt_amp_day1_big_mean.',...
                                    mvmnt_amp_testphase_big_mean.');

T_mean_mvmnt_amp_big = table(mean_mvmnt_amp_big(:,1),...
                                mean_mvmnt_amp_big(:,2),...
                                mean_mvmnt_amp_big(:,3),...
    'VariableNames',{'night','day', 'testphase'});
                 

 
% Number of movements/fish
mvmnt_num_night1_small = cellfun(@numel, mvmnt_amp_night1_small);
mvmnt_num_day1_small = cellfun(@numel, mvmnt_amp_day1_small);
mvmnt_num_testphase_small = cellfun(@numel, mvmnt_amp_testphase_small);

mvmnt_num_night1_big = cellfun(@numel, mvmnt_amp_night1_big);
mvmnt_num_day1_big = cellfun(@numel, mvmnt_amp_day1_big);
mvmnt_num_testphase_big = cellfun(@numel, mvmnt_amp_testphase_big);

% Compile movement number data into table
mvmnt_num_small = catpad(2, mvmnt_num_night1_small.',...
                                  mvmnt_num_day1_small.',...
                                  mvmnt_num_testphase_small.');

T_mvmnt_num_small = table(mvmnt_num_small(:,1), mvmnt_num_small(:,2), mvmnt_num_small(:,3),...
    'VariableNames',{'night','day', 'testphase'});

mvmnt_num_big = catpad(2, mvmnt_num_night1_big.',...
                                  mvmnt_num_day1_big.',...
                                  mvmnt_num_testphase_big.');

T_mvmnt_num_big = table(mvmnt_num_big(:,1), mvmnt_num_big(:,2), mvmnt_num_big(:,3),...
    'VariableNames',{'night','day', 'testphase'});

                  
% Mean movement length and intermovement interval length/fish
mvmnt_len_night1_small_mean = cellfun(@mean, mvmnt_len_night1_small);
mvmnt_len_night1_big_mean = cellfun(@mean, mvmnt_len_night1_big);
mvmnt_len_day1_small_mean = cellfun(@mean, mvmnt_len_day1_small);
mvmnt_len_day1_big_mean = cellfun(@mean, mvmnt_len_day1_big);
mvmnt_len_testphase_small_mean = cellfun(@mean, mvmnt_len_testphase_small);
mvmnt_len_testphase_big_mean = cellfun(@mean, mvmnt_len_testphase_big);


mean_mvmnt_len_small = catpad(2, mvmnt_len_night1_small_mean.',...
                                       mvmnt_len_day1_small_mean.',...
                                       mvmnt_len_testphase_small_mean.');

mean_mvmnt_len_big = catpad(2, mvmnt_len_night1_big_mean.',...
                                       mvmnt_len_day1_big_mean.',...
                                       mvmnt_len_testphase_big_mean.');


% Convert all the lengths to seconds, using average framerate                              
mean_mvmnt_len_small = mean_mvmnt_len_small/framerate;                     
mean_mvmnt_len_big = mean_mvmnt_len_big/framerate;

T_mean_mvmnt_len_small = table(mean_mvmnt_len_small(:,1),...
                                mean_mvmnt_len_small(:,2),...
                                mean_mvmnt_len_small(:,3),...
    'VariableNames',{'night','day', 'testphase'});

T_mean_mvmnt_len_big = table(mean_mvmnt_len_big(:,1),...
                                mean_mvmnt_len_big(:,2),...
                                mean_mvmnt_len_big(:,3),...
    'VariableNames',{'night','day', 'testphase'});

                 
% Intermovement Intervals           
% Remove all intermovement intervals of only 1 frame:
clear x
imi_len_night1 = cellfun(@(x) x(x > 1), imi_len_night1, 'UniformOutput', 0);
imi_len_day1 = cellfun(@(x) x(x > 1), imi_len_day1, 'UniformOutput', 0);
imi_len_testphase = cellfun(@(x) x(x > 1), imi_len_testphase, 'UniformOutput', 0);

% Also remove all sleep bouts. Theoretically this is covered by rest bout
% analysis above. 
imi_len_night1 = cellfun(@(x) x(x < framerate*60), imi_len_night1, 'UniformOutput', 0);
imi_len_day1 = cellfun(@(x) x(x < framerate*60), imi_len_day1, 'UniformOutput', 0);
imi_len_testphase = cellfun(@(x) x(x < framerate*60), imi_len_testphase, 'UniformOutput', 0);

% Mean intermovement interval/fish
imi_len_night1_mean = cellfun(@mean, imi_len_night1);
imi_len_day1_mean = cellfun(@mean, imi_len_day1);
imi_len_testphase_mean = cellfun(@mean, imi_len_testphase);

mean_imi_len = catpad(2, imi_len_night1_mean.',...
                               imi_len_day1_mean.',...
                               imi_len_testphase_mean.');
             
% Convert all the length tables to seconds, using average framerate
mean_imi_len = mean_imi_len/framerate;

% Gather tables
T_mean_imi_len = table(mean_imi_len(:,1), mean_imi_len(:,2), mean_imi_len(:,3),...
    'VariableNames',{'night','day', 'testphase'});
    
disp('Gathering movement bout data...')
                 
% Export all movement bout data:
writetable(T_mean_imi_len, 'output_mean_imi_len.csv');
writetable(T_mean_mvmnt_amp_big, 'output_mean_mvmnt_amp_big.csv');
writetable(T_mean_mvmnt_amp_small, 'output_mean_mvmnt_amp_small.csv');
writetable(T_mean_mvmnt_len_big, 'output_mean_mvmnt_len_big.csv');
writetable(T_mean_mvmnt_len_small, 'output_mean_mvmnt_len_small.csv');
writetable(T_mvmnt_num_big, 'output_mvmnt_num_big.csv');
writetable(T_mvmnt_num_small, 'output_mvmnt_num_small.csv');
disp('Done')

%% Movement plots: Mean Data

% Mean Intermovement Interval
figure
subplot(1,3,1)
CategoricalScatterplot(T_mean_imi_len_night{:,{'ctrl_n1', 'hom_n1', 'het_n1'}}, 'Labels',{'ctrl_n1', 'hom_n1', 'het_n1'});
title('Night1')
ylabel('Intermovement Interval (sec)')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,2)
CategoricalScatterplot(T_mean_imi_len_night{:,{'ctrl_n2', 'hom_n2', 'het_n2'}}, 'Labels',{'ctrl_n2', 'hom_n2', 'het_n2'});
title('Night2')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,3)
CategoricalScatterplot(T_mean_imi_len_night{:,{'ctrl_n3', 'hom_n3', 'het_n3'}}, 'Labels',{'ctrl_n3', 'hom_n3', 'het_n3'});
title('Night3')
ylim([0,inf])
set(gca,'FontSize', 12, 'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_imi_len_night', 'svg'); % save as svg
saveas(gcf,'mean_imi_len_night', 'png'); % save as png

figure
subplot(1,2,1)
CategoricalScatterplot(T_mean_imi_len_day{:,{'ctrl_d1', 'hom_d1', 'het_d1'}}, 'Labels',{'ctrl_d1', 'hom_d1', 'het_d1'});
title('Day1')
ylabel('Intermovement Interval (sec)')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,2,2)
CategoricalScatterplot(T_mean_imi_len_day{:,{'ctrl_d2', 'hom_d2', 'het_d2'}}, 'Labels',{'ctrl_d2', 'hom_d2', 'het_d2'});
title('Day2')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_imi_len_day', 'svg'); % save as svg
saveas(gcf,'mean_imi_len_day', 'png'); % save as png
savefig('mean_imi_len_day.fig') % save matlab figure in case you need to reformat


% Mean Movement amplitude
figure
subplot(1,3,1)
CategoricalScatterplot(T_mean_mvmnt_amp_night_big{:,{'ctrl_n1', 'hom_n1', 'het_n1'}}, 'Labels',{'ctrl_n1', 'hom_n1', 'het_n1'});
title('night1')
ylabel('Movement Amplitude (pixels)')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,2)
CategoricalScatterplot(T_mean_mvmnt_amp_night_big{:,{'ctrl_n2', 'hom_n2', 'het_n2'}}, 'Labels',{'ctrl_n2', 'hom_n2', 'het_n2'});
title('night2')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,3)
CategoricalScatterplot(T_mean_mvmnt_amp_night_big{:,{'ctrl_n3', 'hom_n3', 'het_n3'}}, 'Labels',{'ctrl_n3', 'hom_n3', 'het_n3'});
title('night3')
ylim([0,inf])
set(gca,'FontSize', 12, 'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_mvmnt_amp_night_big', 'svg'); % save as svg
saveas(gcf,'mean_mvmnt_amp_night_big', 'png'); % save as png
savefig('mean_mvmnt_amp_night_big.fig') % save matlab figure in case you need to reformat

figure
subplot(1,3,1)
CategoricalScatterplot(T_mean_mvmnt_amp_night_small{:,{'ctrl_n1', 'hom_n1', 'het_n1'}}, 'Labels',{'ctrl_n1', 'hom_n1', 'het_n1'});
title('night1')
ylabel('Movement Amplitude (pixels)')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,2)
CategoricalScatterplot(T_mean_mvmnt_amp_night_small{:,{'ctrl_n2', 'hom_n2', 'het_n2'}}, 'Labels',{'ctrl_n2', 'hom_n2', 'het_n2'});
title('night2')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,3)
CategoricalScatterplot(T_mean_mvmnt_amp_night_small{:,{'ctrl_n3', 'hom_n3', 'het_n3'}}, 'Labels',{'ctrl_n3', 'hom_n3', 'het_n3'});
title('night3')
ylim([0,inf])
set(gca,'FontSize', 12, 'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_mvmnt_amp_night_small', 'svg'); % save as svg
saveas(gcf,'mean_mvmnt_amp_night_small', 'png'); % save as png
savefig('mean_mvmnt_amp_night_small.fig') % save matlab figure in case you need to reformat

figure
subplot(1,2,1)
CategoricalScatterplot(T_mean_mvmnt_amp_day_big{:,{'ctrl_d1', 'hom_d1', 'het_d1'}}, 'Labels',{'ctrl_d1', 'hom_d1', 'het_d1'});
title('Day1')
ylabel('Movement Amplitude (pixels)')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,2,2)
CategoricalScatterplot(T_mean_mvmnt_amp_day_big{:,{'ctrl_d2', 'hom_d2', 'het_d2'}}, 'Labels',{'ctrl_d2', 'hom_d2', 'het_d2'});
title('Day2')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_mvmnt_amp_day_big', 'svg'); % save as svg
saveas(gcf,'mean_mvmnt_amp_day_big', 'png'); % save as png
savefig('mean_mvmnt_amp_day_big.fig') % save matlab figure in case you need to reformat

figure
subplot(1,2,1)
CategoricalScatterplot(T_mean_mvmnt_amp_day_small{:,{'ctrl_d1', 'hom_d1', 'het_d1'}}, 'Labels',{'ctrl_d1', 'hom_d1', 'het_d1'});
title('Day1')
ylabel('Movement Amplitude (pixels)')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,2,2)
CategoricalScatterplot(T_mean_mvmnt_amp_day_small{:,{'ctrl_d2', 'hom_d2', 'het_d2'}}, 'Labels',{'ctrl_d2', 'hom_d2', 'het_d2'});
title('Day2')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_mvmnt_amp_day_small', 'svg'); % save as svg
saveas(gcf,'mean_mvmnt_amp_day_small', 'png'); % save as png
savefig('mean_mvmnt_amp_day_small.fig') % save matlab figure in case you need to reformat


% Mean movement length
figure
subplot(1,3,1)
CategoricalScatterplot(T_mean_mvmnt_len_night_big{:,{'ctrl_n1', 'hom_n1', 'het_n1'}}, 'Labels',{'ctrl_n1', 'hom_n1', 'het_n1'});
title('night1')
ylabel('Movement Length (sec)')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,2)
CategoricalScatterplot(T_mean_mvmnt_len_night_big{:,{'ctrl_n2', 'hom_n2', 'het_n2'}}, 'Labels',{'ctrl_n2', 'hom_n2', 'het_n2'});
title('night2')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,3)
CategoricalScatterplot(T_mean_mvmnt_len_night_big{:,{'ctrl_n3', 'hom_n3', 'het_n3'}}, 'Labels',{'ctrl_n3', 'hom_n3', 'het_n3'});
title('night3')
ylim([0,inf])
set(gca,'FontSize', 12, 'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_mvmnt_len_night_big', 'svg'); % save as svg
saveas(gcf,'mean_mvmnt_len_night_big', 'png'); % save as png
savefig('mean_mvmnt_len_night_big.fig') % save matlab figure in case you need to reformat

figure
subplot(1,3,1)
CategoricalScatterplot(T_mean_mvmnt_len_night_small{:,{'ctrl_n1', 'hom_n1', 'het_n1'}}, 'Labels',{'ctrl_n1', 'hom_n1', 'het_n1'});
title('night1')
ylabel('Movement Length (sec)')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,2)
CategoricalScatterplot(T_mean_mvmnt_len_night_small{:,{'ctrl_n2', 'hom_n2', 'het_n2'}}, 'Labels',{'ctrl_n2', 'hom_n2', 'het_n2'});
title('night2')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,3)
CategoricalScatterplot(T_mean_mvmnt_len_night_small{:,{'ctrl_n3', 'hom_n3', 'het_n3'}}, 'Labels',{'ctrl_n3', 'hom_n3', 'het_n3'});
title('night3')
ylim([0,inf])
set(gca,'FontSize', 12, 'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_mvmnt_len_night_small', 'svg'); % save as svg
saveas(gcf,'mean_mvmnt_len_night_small', 'png'); % save as png
savefig('mean_mvmnt_len_night_small.fig') % save matlab figure in case you need to reformat

figure
subplot(1,2,1)
CategoricalScatterplot(T_mean_mvmnt_len_day_big{:,{'ctrl_d1', 'hom_d1', 'het_d1'}}, 'Labels',{'ctrl_d1', 'hom_d1', 'het_d1'});
title('Day1')
ylabel('Movement Length (sec)')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,2,2)
CategoricalScatterplot(T_mean_mvmnt_len_day_big{:,{'ctrl_d2', 'hom_d2', 'het_d2'}}, 'Labels',{'ctrl_d2', 'hom_d2', 'het_d2'});
title('Day2')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_mvmnt_len_day_big', 'svg'); % save as svg
saveas(gcf,'mean_mvmnt_len_day_big', 'png'); % save as png
savefig('mean_mvmnt_len_day_big.fig') % save matlab figure in case you need to reformat

figure
subplot(1,2,1)
CategoricalScatterplot(T_mean_mvmnt_len_day_small{:,{'ctrl_d1', 'hom_d1', 'het_d1'}}, 'Labels',{'ctrl_d1', 'hom_d1', 'het_d1'});
title('Day1')
ylabel('Movement Length (sec)')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,2,2)
CategoricalScatterplot(T_mean_mvmnt_len_day_small{:,{'ctrl_d2', 'hom_d2', 'het_d2'}}, 'Labels',{'ctrl_d2', 'hom_d2', 'het_d2'});
title('Day2')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_mvmnt_len_day_small', 'svg'); % save as svg
saveas(gcf,'mean_mvmnt_len_day_small', 'png'); % save as png
savefig('mean_mvmnt_len_day_small.fig') % save matlab figure in case you need to reformat


% Movement number
figure
subplot(1,3,1)
CategoricalScatterplot(T_mvmnt_num_night_big{:,{'ctrl_n1', 'hom_n1', 'het_n1'}}, 'Labels',{'ctrl_n1', 'hom_n1', 'het_n1'});
title('night1')
ylabel('# Movements ')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,2)
CategoricalScatterplot(T_mvmnt_num_night_big{:,{'ctrl_n2', 'hom_n2', 'het_n2'}}, 'Labels',{'ctrl_n2', 'hom_n2', 'het_n2'});
title('night2')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,3)
CategoricalScatterplot(T_mvmnt_num_night_big{:,{'ctrl_n3', 'hom_n3', 'het_n3'}}, 'Labels',{'ctrl_n3', 'hom_n3', 'het_n3'});
title('night3')
ylim([0,inf])
set(gca,'FontSize', 12, 'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_mvmnt_num_night_big', 'svg'); % save as svg
saveas(gcf,'mean_mvmnt_num_night_big', 'png'); % save as png
savefig('mean_mvmnt_num_night_big.fig') % save matlab figure in case you need to reformat

figure
subplot(1,3,1)
CategoricalScatterplot(T_mvmnt_num_night_small{:,{'ctrl_n1', 'hom_n1', 'het_n1'}}, 'Labels',{'ctrl_n1', 'hom_n1', 'het_n1'});
title('night1')
ylabel('# Movements')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,2)
CategoricalScatterplot(T_mvmnt_num_night_small{:,{'ctrl_n2', 'hom_n2', 'het_n2'}}, 'Labels',{'ctrl_n2', 'hom_n2', 'het_n2'});
title('night2')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,3,3)
CategoricalScatterplot(T_mvmnt_num_night_small{:,{'ctrl_n3', 'hom_n3', 'het_n3'}}, 'Labels',{'ctrl_n3', 'hom_n3', 'het_n3'});
title('night3')
ylim([0,inf])
set(gca,'FontSize', 12, 'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_mvmnt_num_night_small', 'svg'); % save as svg
saveas(gcf,'mean_mvmnt_num_night_small', 'png'); % save as png
savefig('mean_mvmnt_num_night_small.fig') % save matlab figure in case you need to reformat

figure
subplot(1,2,1)
CategoricalScatterplot(T_mvmnt_num_day_big{:,{'ctrl_d1', 'hom_d1', 'het_d1'}}, 'Labels',{'ctrl_d1', 'hom_d1', 'het_d1'});
title('Day1')
ylabel('# Movements')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,2,2)
CategoricalScatterplot(T_mvmnt_num_day_big{:,{'ctrl_d2', 'hom_d2', 'het_d2'}}, 'Labels',{'ctrl_d2', 'hom_d2', 'het_d2'});
title('Day2')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_mvmnt_num_day_big', 'svg'); % save as svg
saveas(gcf,'mean_mvmnt_num_day_big', 'png'); % save as png
savefig('mean_mvmnt_num_day_big.fig') % save matlab figure in case you need to reformat

figure
subplot(1,2,1)
CategoricalScatterplot(T_mvmnt_num_day_small{:,{'ctrl_d1', 'hom_d1', 'het_d1'}}, 'Labels',{'ctrl_d1', 'hom_d1', 'het_d1'});
title('Day1')
ylabel('# Movements')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
subplot(1,2,2)
CategoricalScatterplot(T_mvmnt_num_day_small{:,{'ctrl_d2', 'hom_d2', 'het_d2'}}, 'Labels',{'ctrl_d2', 'hom_d2', 'het_d2'});
title('Day2')
ylim([0,inf])
set(gca,'FontSize', 12,'XTickLabelRotation',45, 'TickLabelInterpreter','none')
saveas(gcf,'mean_mvmnt_num_day_small', 'svg'); % save as svg
saveas(gcf,'mean_mvmnt_num_day_small', 'png'); % save as png
savefig('mean_mvmnt_num_day_small.fig') % save matlab figure in case you need to reformat


% Plotting with default matlab boxplot function
% figure
% subplot(1,2,1)
% boxplot(T_mean_imi_len_day{:,{'ctrl_d1', 'hom_d1', 'het_d1'}}, 'Labels',{'ctrl_d1', 'hom_d1', 'het_d1'});
% title('Day1')
% ylabel('Intermovement Interval (sec)')
% set(gca,'XTickLabelRotation',45)
% subplot(1,2,2)
% boxplot(T_mean_imi_len_day{:,{'ctrl_d2', 'hom_d2', 'het_d2'}}, 'Labels',{'ctrl_d2', 'hom_d2', 'het_d2'});
% title('Day2')
% set(gca,'XTickLabelRotation',45)
% 
% figure
% subplot(1,3,1)
% boxplot(T_mean_imi_len_night{:,{'ctrl_n1', 'hom_n1', 'het_n1'}}, 'Labels',{'ctrl_n1', 'hom_n1', 'het_n1'});
% title('Night1')
% set(gca,'XTickLabelRotation',45)
% ylabel('Intermovement Interval (sec)')
% subplot(1,3,2)
% boxplot(T_mean_imi_len_night{:,{'ctrl_n2', 'hom_n2', 'het_n2'}}, 'Labels',{'ctrl_n2', 'hom_n2', 'het_n2'});
% title('Night2')
% set(gca,'XTickLabelRotation',45)
% subplot(1,3,3)
% boxplot(T_mean_imi_len_night{:,{'ctrl_n3', 'hom_n3', 'het_n3'}}, 'Labels',{'ctrl_n3', 'hom_n3', 'het_n3'});
% title('Night3')
% set(gca,'XTickLabelRotation',45)


%% Movement histograms: Movement Amplitudes
% histogram(test_amplitude, 'Normalization', 'cumcount', 'DisplayStyle', 'stairs');
% histogram(test_amplitude, 'Normalization', 'cdf','DisplayStyle', 'stairs');

% Histogram of movement amplitudes, combined across all fish/genotype
mvmnt_amp_all_night1_big = vertcat(mvmnt_amp_night1_big{:});
mvmnt_amp_all_day1_big = vertcat(mvmnt_amp_day1_big{:});
mvmnt_amp_all_testphase_big = vertcat(mvmnt_amp_testphase_big{:});

figure
h1=histogram(mvmnt_amp_all_night1_big, 'Normalization', 'probability', 'FaceColor', '[0.4,0.4,0.4]', 'BinWidth', 20);
title('mvmnt_amp_night_big', 'Interpreter', 'none')
saveas(gcf,'mvmnt_amp_night_big', 'png'); % save as png

figure
h1=histogram(mvmnt_amp_all_day1_big, 'Normalization', 'probability', 'FaceColor', '[0.4,0.4,0.4]', 'BinWidth', 20);
title('mvmnt_amp_day_big', 'Interpreter', 'none')
saveas(gcf,'mvmnt_amp_day_big', 'png'); % save as png

figure
h1=histogram(mvmnt_amp_all_testphase_big, 'Normalization', 'probability', 'FaceColor', '[0.4,0.4,0.4]', 'BinWidth', 20);
title('mvmnt_amp_testphase_big', 'Interpreter', 'none')
saveas(gcf,'mvmnt_amp_testphase_big', 'png'); % save as png

figure
h1=histogram(mvmnt_amp_all_day1_big, 'Normalization', 'probability', 'FaceColor', '[0.4,0.4,0.4]', 'BinWidth', 20);
hold on
h2=histogram(mvmnt_amp_all_testphase_big, 'Normalization', 'probability', 'FaceColor', '[1,0.5,0.5]', 'BinWidth', 20);



%% Movement histograms: Movement Lengths

% Histogram of overall movement lengths, converted to seconds
mvmnt_len_all_night1_ctrl_big = vertcat(mvmnt_len_night1_ctrl_big{:})/framerate;
mvmnt_len_all_night1_hom_big = vertcat(mvmnt_len_night1_hom_big{:})/framerate;
mvmnt_len_all_night1_het_big = vertcat(mvmnt_len_night1_het_big{:})/framerate;
mvmnt_len_all_night2_ctrl_big = vertcat(mvmnt_len_night2_ctrl_big{:})/framerate;
mvmnt_len_all_night2_hom_big = vertcat(mvmnt_len_night2_hom_big{:})/framerate;
mvmnt_len_all_night2_het_big = vertcat(mvmnt_len_night2_het_big{:})/framerate;
mvmnt_len_all_night3_ctrl_big = vertcat(mvmnt_len_night3_ctrl_big{:})/framerate;
mvmnt_len_all_night3_hom_big = vertcat(mvmnt_len_night3_hom_big{:})/framerate;
mvmnt_len_all_night3_het_big = vertcat(mvmnt_len_night3_het_big{:})/framerate;
mvmnt_len_all_day1_ctrl_big = vertcat(mvmnt_len_day1_ctrl_big{:})/framerate;
mvmnt_len_all_day1_hom_big = vertcat(mvmnt_len_day1_hom_big{:})/framerate;
mvmnt_len_all_day1_het_big = vertcat(mvmnt_len_day1_het_big{:})/framerate;
mvmnt_len_all_day2_ctrl_big = vertcat(mvmnt_len_day2_ctrl_big{:})/framerate;
mvmnt_len_all_day2_hom_big = vertcat(mvmnt_len_day2_hom_big{:})/framerate;
mvmnt_len_all_day2_het_big = vertcat(mvmnt_len_day2_het_big{:})/framerate;

figure
h1=histogram(mvmnt_len_all_night1_ctrl_big, 'Normalization','probability', 'FaceColor', '[0.4,0.4,0.4]', 'BinWidth', 1/framerate);
hold on
h2=histogram(mvmnt_len_all_night1_hom_big, 'Normalization','probability', 'FaceColor', '[1,0.5,0.5]', 'BinWidth', 1/framerate);
hold on
h3=histogram(mvmnt_len_all_night1_het_big, 'Normalization','probability', 'FaceColor', '[0.85,0.85,0.85]', 'BinWidth', 1/framerate);
title('mvmnt_len_all_night1_big', 'Interpreter', 'none')
saveas(gcf,'mvmnt_len_all_night1_big', 'svg'); % save as svg
saveas(gcf,'mvmnt_len_all_night1_big', 'png');
savefig('mvmnt_len_all_night1_big.fig'); % save matlab figure in case you need to reformat

figure
h1=histogram(mvmnt_len_all_night2_ctrl_big, 'Normalization','probability', 'FaceColor', '[0.4,0.4,0.4]', 'BinWidth', 1/framerate);
hold on
h2=histogram(mvmnt_len_all_night2_hom_big, 'Normalization','probability', 'FaceColor', '[1,0.5,0.5]', 'BinWidth', 1/framerate);
hold on
h3=histogram(mvmnt_len_all_night2_het_big, 'Normalization','probability', 'FaceColor', '[0.85,0.85,0.85]', 'BinWidth', 1/framerate);
title('mvmnt_len_all_night2_big', 'Interpreter', 'none')
saveas(gcf,'mvmnt_len_all_night2_big', 'svg'); % save as svg
saveas(gcf,'mvmnt_len_all_night2_big', 'png');
savefig('mvmnt_len_all_night2_big.fig'); % save matlab figure in case you need to reformat

figure
h1=histogram(mvmnt_len_all_night3_ctrl_big, 'Normalization','probability', 'FaceColor', '[0.4,0.4,0.4]', 'BinWidth', 1/framerate);
hold on
h2=histogram(mvmnt_len_all_night3_hom_big, 'Normalization','probability', 'FaceColor', '[1,0.5,0.5]', 'BinWidth', 1/framerate);
hold on
h3=histogram(mvmnt_len_all_night3_het_big, 'Normalization','probability', 'FaceColor', '[0.85,0.85,0.85]', 'BinWidth', 1/framerate);
title('mvmnt_len_all_night3_big', 'Interpreter', 'none')
saveas(gcf,'mvmnt_len_all_night3_big', 'svg'); % save as svg
saveas(gcf,'mvmnt_len_all_night3_big', 'png');
savefig('mvmnt_len_all_night3_big.fig'); % save matlab figure in case you need to reformat

figure
h1=histogram(mvmnt_len_all_day1_ctrl_big, 'Normalization','probability', 'FaceColor', '[0.4,0.4,0.4]', 'BinWidth', 1/framerate);
hold on
h2=histogram(mvmnt_len_all_day1_hom_big, 'Normalization','probability', 'FaceColor', '[1,0.5,0.5]', 'BinWidth', 1/framerate);
hold on
h3=histogram(mvmnt_len_all_day1_het_big, 'Normalization','probability', 'FaceColor', '[0.85,0.85,0.85]', 'BinWidth', 1/framerate);
title('mvmnt_len_all_day1_big', 'Interpreter', 'none')
saveas(gcf,'mvmnt_len_all_day1_big', 'svg'); % save as svg
saveas(gcf,'mvmnt_len_all_day1_big', 'png');
savefig('mvmnt_len_all_day1_big.fig'); % save matlab figure in case you need to reformat

figure
h1=histogram(mvmnt_len_all_day2_ctrl_big, 'Normalization','probability', 'FaceColor', '[0.4,0.4,0.4]', 'BinWidth', 1/framerate);
hold on
h2=histogram(mvmnt_len_all_day2_hom_big, 'Normalization','probability', 'FaceColor', '[1,0.5,0.5]', 'BinWidth', 1/framerate);
hold on
h3=histogram(mvmnt_len_all_day2_het_big, 'Normalization','probability', 'FaceColor', '[0.85,0.85,0.85]', 'BinWidth', 1/framerate);
title('mvmnt_len_all_day2_big', 'Interpreter', 'none')
saveas(gcf,'mvmnt_len_all_day2_big', 'svg'); % save as svg
saveas(gcf,'mvmnt_len_all_day2_big', 'png');
savefig('mvmnt_len_all_day2_big.fig'); % save matlab figure in case you need to reformat


%% Movement histograms: intermovement intervals
% Get all imis across all fish/genotype, in seconds:
imi_all_night1_ctrl = vertcat(imi_len_night1_ctrl{:})/framerate;
imi_all_night1_hom = vertcat(imi_len_night1_hom{:})/framerate;
imi_all_night1_het = vertcat(imi_len_night1_het{:})/framerate;
imi_all_night2_ctrl = vertcat(imi_len_night2_ctrl{:})/framerate;
imi_all_night2_hom = vertcat(imi_len_night2_hom{:})/framerate;
imi_all_night2_het = vertcat(imi_len_night2_het{:})/framerate;
imi_all_night3_ctrl = vertcat(imi_len_night3_ctrl{:})/framerate;
imi_all_night3_hom = vertcat(imi_len_night3_hom{:})/framerate;
imi_all_night3_het = vertcat(imi_len_night3_het{:})/framerate;
imi_all_day1_ctrl = vertcat(imi_len_day1_ctrl{:})/framerate;
imi_all_day1_hom = vertcat(imi_len_day1_hom{:})/framerate;
imi_all_day1_het = vertcat(imi_len_day1_het{:})/framerate;
imi_all_day2_ctrl = vertcat(imi_len_day2_ctrl{:})/framerate;
imi_all_day2_hom = vertcat(imi_len_day2_hom{:})/framerate;
imi_all_day2_het = vertcat(imi_len_day2_het{:})/framerate;

figure
h1=histogram(log(imi_all_night1_ctrl), 'Normalization','probability', 'FaceColor', '[0.4,0.4,0.4]', 'BinWidth', 0.25);
hold on
h2=histogram(log(imi_all_night1_hom), 'Normalization','probability', 'FaceColor', '[1,0.5,0.5]', 'BinWidth', 0.25);
hold on
h3=histogram(log(imi_all_night1_het), 'Normalization','probability', 'FaceColor', '[0.85,0.85,0.85]', 'BinWidth', 0.25);
title('imi_all_night1', 'Interpreter', 'none')
saveas(gcf,'imi_all_night1', 'svg'); % save as svg
saveas(gcf,'imi_all_night1', 'png');
savefig('imi_all_night1.fig'); % save matlab figure in case you need to reformat

figure
h1=histogram(log(imi_all_night2_ctrl), 'Normalization','probability', 'FaceColor', '[0.4,0.4,0.4]', 'BinWidth', 0.25);
hold on
h2=histogram(log(imi_all_night2_hom), 'Normalization','probability', 'FaceColor', '[1,0.5,0.5]', 'BinWidth', 0.25);
hold on
h3=histogram(log(imi_all_night2_het), 'Normalization','probability', 'FaceColor', '[0.85,0.85,0.85]', 'BinWidth', 0.25);
title('imi_all_night2', 'Interpreter', 'none')
saveas(gcf,'imi_all_night2', 'svg'); % save as svg
saveas(gcf,'imi_all_night2', 'png');
savefig('imi_all_night2.fig'); % save matlab figure in case you need to reformat

figure
h1=histogram(log(imi_all_night3_ctrl), 'Normalization','probability', 'FaceColor', '[0.4,0.4,0.4]', 'BinWidth', 0.25);
hold on
h2=histogram(log(imi_all_night3_hom), 'Normalization','probability', 'FaceColor', '[1,0.5,0.5]', 'BinWidth', 0.25);
hold on
h3=histogram(log(imi_all_night3_het), 'Normalization','probability', 'FaceColor', '[0.85,0.85,0.85]', 'BinWidth', 0.25);
title('imi_all_night3', 'Interpreter', 'none')
saveas(gcf,'imi_all_night3', 'svg'); % save as svg
saveas(gcf,'imi_all_night3', 'png');
savefig('imi_all_night3.fig'); % save matlab figure in case you need to reformat

figure
h1=histogram(log(imi_all_day1_ctrl), 'Normalization','probability', 'FaceColor', '[0.4,0.4,0.4]', 'BinWidth', 0.25);
hold on
h2=histogram(log(imi_all_day1_hom), 'Normalization','probability', 'FaceColor', '[1,0.5,0.5]', 'BinWidth', 0.25);
hold on
h3=histogram(log(imi_all_day1_het), 'Normalization','probability', 'FaceColor', '[0.85,0.85,0.85]', 'BinWidth', 0.25);
title('imi_all_day1', 'Interpreter', 'none')
saveas(gcf,'imi_all_day1', 'svg'); % save as svg
saveas(gcf,'imi_all_day1', 'png');
savefig('imi_all_day1.fig'); % save matlab figure in case you need to reformat

figure
h1=histogram(log(imi_all_day2_ctrl), 'Normalization','probability', 'FaceColor', '[0.4,0.4,0.4]', 'BinWidth', 0.25);
hold on
h2=histogram(log(imi_all_day2_hom), 'Normalization','probability', 'FaceColor', '[1,0.5,0.5]', 'BinWidth', 0.25);
hold on
h3=histogram(log(imi_all_day2_het), 'Normalization','probability', 'FaceColor', '[0.85,0.85,0.85]', 'BinWidth', 0.25);
title('imi_all_day2', 'Interpreter', 'none')
saveas(gcf,'imi_all_day2', 'svg'); % save as svg
saveas(gcf,'imi_all_day2', 'png');
savefig('imi_all_day2.fig'); % save matlab figure in case you need to reformat


% If plotting by frame, not seconds:
% figure
% h=histogram(log(imi_all_day1_ctrl));
% h.NumBins = 70;
% hold on;
% line([log(framerate), log(framerate)], ylim, 'LineWidth', 0.5, 'LineStyle','--', 'Color', 'k'); %line at 1sec
% hold on;
% line([log(framerate*60), log(framerate*60)], ylim, 'LineWidth', 0.5, 'LineStyle','--', 'Color', 'k'); %line at 1min
% title('imi_all_day1_ctrl', 'Interpreter', 'none')

