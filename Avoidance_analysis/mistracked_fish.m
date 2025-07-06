%% 1) this code removes all fish that were not moving 70% of the time and fish that were not tracked properly
% run this code section by section, using "Ctrl" + "Enter"
% a) Filter out fish with low movement
noMove = [];

AllFish2_filtered = AllFish2;  % Initialize the filtered variable

% Filter fish with low movement (more than 70% time spent with zero movement)
for i = 1:size(AllFish2_filtered, 2)
    dist = sqrt(sum(diff(AllFish2_filtered(i).data(:, [1, 2])).^2, 2));
    dist = smooth(dist, 20);
    noMove(i) = sum(dist == 0) / length(dist);
end
AllFish2_filtered(noMove > 0.7) = [];  % Apply filtering to remove low movement fish

% b) check again if there are fish with less than 500 bouting events. Too
%little

% Initialize an empty array to store fish that pass the filter
%FilteredFish = [];

%for i = 1:size(AllFish2_filtered, 2)
    % Check if the fish has at least 500 rows in the data field
    %if size(AllFish2_filtered(i).preBout_x, 1) >= 500
        %FilteredFish = [FilteredFish, AllFish2_filtered(i)];
    %end
%end

% Replace AllFish2_filtered with the filtered version (optional)
%AllFish2_filtered = FilteredFish;

%save the file: paste the path
pathname = 'E:\Olfaction\Preference_index_cortisol\100ng_cortisol_15min';
save(strcat(pathname,'\cortisol_AllFish_noMove.mat'),'AllFish2_filtered','-v7.3');


%% 2) check the X and Y positions of each fish in a plot to determine whether fish have
% been mistracked, i.e. if fluff or a part of the wall have been tracked
% instead of the fish

% Plot Heading angle over time
% prebout_x and prebout_y
% left side: low relative odor concentration. right side: high relative
% odor concentration

u=36;
% change fish number in brackets
x=AllFish2_filtered(u).postBout_x; % check tracking of every single fish.
y=AllFish2_filtered(u).postBout_y; % check tracking of every single fish
figure;
plot(x,y);
xlim([330 900]);
ylim([260 777]);

xlabel('X position');
ylabel('Y position');
title('X and y position in the swimming lane');

% go through every single fish and not down fish with tracking errors (long
% straight jumps along x or y axis) below for removal

%% 3) this function removes fish from the data analysis where the tracking was
% visibly bad
% where fish were not tracked properly for more than 70% of the time
% (indicated by long straight line,s where the tracking jumps visibly on
% the X or Y axis)

% remove fish with bad tracking by hand
% uncommennt one by one
AllFish2_filtered2=AllFish2_filtered

% 0mM cadaverine NEW
%AllFish2_filtered2([20]) = []; %control
%AllFish2_filtered2([3,4,9,13,20,23]) = [];  %sleep deprived

% 1mM cadaverine NEW
%AllFish2_filtered2([6,7,14,21,22,23,28]) = []; %control
%AllFish2_filtered2([5,15,17,18]) = [];  %sleep deprived

% 2.5mM cadaverine NEW
%AllFish2_filtered2([4,8,13,16,17,20,31,32,37,38,39]) = []; %control
%AllFish2_filtered2([2,10,11,12,18,19,29,30,33,37,39]) = [];  %sleep deprived

% 5mM cadaverine NEW
%AllFish2_filtered2([3,5,21,26,27,33,]) = []; %control
%AllFish2_filtered2([2,6,8,9,20,33,40,]) = [];  %sleep deprived


% 10mM putrescine
%AllFish2_filtered2([2,24,28,29,31,35,36,44,]) = []; %control
%AllFish2_filtered2([13,16,19,29,36]) = [];  %sleep deprived

% 100mM salt NEW
%AllFish2_filtered2([1,6,9,25,28,32,33]) = []; %control
%AllFish2_filtered2([3,4,15,18,25,26,28,29,33,]) = [];  %sleep deprived

% 5 ng/ml cortisol (in 2.5mM cadaverine) NEW
%AllFish2_filtered2([4,11,18,19,30]) = []; %control
%AllFish2_filtered2([21,31,34]) = [];  %cortisol treated

% 100 ng/ml cortisol (in 2.5mM cadaverine) NEW
%AllFish2_filtered2([30,33]) = []; %control
%AllFish2_filtered2([7,31]) = [];  %cortisol treated

% 5ug/ml cortisol (in 2.5mM cadaverine) NEW
%AllFish2_filtered2([23,29,30,33,34,37,38,39,41,45,47,48]) = []; %control
%AllFish2_filtered2([5,7,29,30,33,37,38,39]) = [];  %cortisol treated


% complete starvation (in 2.5mM cadaverine)
%AllFish2_filtered2([4,]) = []; %control
%AllFish2_filtered2([]) = [];  %starved

% partial starvation (in 2.5mM cadaverine)
%AllFish2_filtered2([9,10,14,24,25,30,31,32]) = []; %control
%AllFish2_filtered2([8,11,16,17,20,21,28,29,]) = [];  %starved

% 100nM melatonin (in 2.5mM cadaverine) NEW
%AllFish2_filtered2([12,15,]) = []; %control
%AllFish2_filtered2([5,10,15]) = [];  %melatonin

% alanine SD  NEW
%AllFish2_filtered2([8,11,14,15,27,]) = []; %control
AllFish2_filtered2([3,4,6,9,13,33,34,]) = [];  %SD


%save the file: paste the path
pathname = 'E:\Olfaction\Preference_index_alanine\compiles';
save(strcat(pathname,'\SD_AllFish_filtered.mat'),'AllFish2_filtered2','-v7.3');

