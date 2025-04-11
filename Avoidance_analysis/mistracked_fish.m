%% 1) this code removes all fish that were not moving 70% of the time and fish that were not tracked properly
% run this code section by section, using "Ctrl" + "Enter"
% a) Filter out fish with low movement
noMove = [];

AllFish2_filtered = AllFish2_filtered;  % Initialize the filtered variable

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

%% 2) check the X and Y positions of each fish in a plot to determine whether fish have
% been mistracked, i.e. if fluff or a part of the wall have been tracked
% instead of the fish

% Plot Heading angle over time
% prebout_x and prebout_y
% left side: low relative odor concentration. right side: high relative
% odor concentration

u=6;
% change fish number in brackets
x=AllFish2_filtered(u).preBout_x; % check tracking of every single fish.
y=AllFish2_filtered(u).preBout_y; % check tracking of every single fish
figure;
plot(x,y);

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
%AllFish2_filtered2([1, 15,16,17,18,19,20,22,24,26,27,28,30,32]) = []; %control
%AllFish2_filtered2([1,3,5,15,16,18,19,20,21,22,23,25,26,32,35,36,37]) = [];  %sleep deprived

% 1mM cadaverine NEW
%AllFish2_filtered2([22,23,24,25,30,31,32,33,36,37,38,39,40,41,42,43,47,48,49]) = []; %control
%AllFish2_filtered2([26,28,29,30,33,39,40,41,42,43,44,48,49]) = [];  %sleep deprived

% 2.5mM cadaverine NEW
%AllFish2_filtered2([12,18,21,24,25,28,32,33,39,40,43,44,45,46,47]) = []; %control
%AllFish2_filtered2([10,12,15,18,19,20,24,26,27,30,37,38,39,41,45,46,47]) = [];  %sleep deprived

% 5mM cadaverine NEW
%AllFish2_filtered2([22,23,25,34,35,36,37,39,41,42,44,45,46,47,49,50,51,53,54]) = []; %control
%AllFish2_filtered2([21,22,24,26,28,29,31,32,35,36,37,38,39,40,41,42,43,44,46,47,48,49,50,51,53,57,60]) = [];  %sleep deprived


% 10mM putrescine
%AllFish2_filtered2([15,19,24,28,29,30,31,34,35,36,37,40]) = []; %control
%AllFish2_filtered2([13,16,20,21,26,27,29,33,36,39,40]) = [];  %sleep deprived

% 100mM salt NEW
%AllFish2_filtered2([9,12,13,14,17,20,22,29,31,33,36,39,40,41,42]) = []; %control
%AllFish2_filtered2([11,12,15,17,20,23,26,27,28,29,30,31,34,35,38,41,43]) = [];  %sleep deprived

% 5 ng/ml cortisol (in 2.5mM cadaverine)
%AllFish2_filtered2([2,4,8,11,18,28,29]) = []; %control
%AllFish2_filtered2([1,2,3,4,8,12,14,15,20,21,26,27,30]) = [];  %cortisol treated

% 100 ng/ml cortisol (in 2.5mM cadaverine) NEW
%AllFish2_filtered2([30,33]) = []; %control
%AllFish2_filtered2([10,31]) = [];  %cortisol treated

% 5ug/ml cortisol (in 2.5mM cadaverine) NEW
%AllFish2_filtered2([23,29,33,34,37,38,40,41,44,45,47,48,49]) = []; %control
%AllFish2_filtered2([5,7,27,29,30,33,34,37,38]) = [];  %cortisol treated


% complete starvation (in 2.5mM cadaverine)
%AllFish2_filtered2([4]) = []; %control
%AllFish2_filtered2([]) = [];  %starved

% partial starvation (in 2.5mM cadaverine)
%AllFish2_filtered2([2,9,10,11,17,18,24,25,29,30,31,32]) = []; %control
%AllFish2_filtered2([8,9,11,16,17,18,20,21,22,24,28,29]) = [];  %starved

% 100nM melatonin (in 2.5mM cadaverine) NEW
%AllFish2_filtered2([13,16,17,18,20,21,23,27,28,31]) = []; %control
%AllFish2_filtered2([16,17,24,26,27,29,30]) = [];  %melatonin


%save the file: paste the path
pathname = 'E:\Harvard_Olfaction\5ugml_cortisol_2,5mMCadaverine\morning_vs_afternoon\afternoon';
save(strcat(pathname,'\control_filtered.mat'),'AllFish2_filtered2','-v7.3');

