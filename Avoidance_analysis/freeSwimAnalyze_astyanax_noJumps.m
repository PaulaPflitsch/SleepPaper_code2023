clear variables %added to clear variables from previous runs, make sure there is no impact from left over variables
[filename,pathname] = uigetfile('*.*','C:\Users\Engert\Dropbox\FreeSwimming','MultiSelect','on');

disp(filename);
disp(size(filename));
disp(class(filename));

if isequal (filename, 0)
    error ('No file selected');
end

if iscell(filename)
    noFiles = numel(filename);
else
    noFiles = 1;
    filename={filename};
end
AllFish = {[]};
Fno = 4;
top = 777;
bottom = 260;
for i = 1:4
    Lane(i,:) = [(bottom + (i-1)*(top-bottom)/4),(260 + i*(top-bottom)/4)];
end
parfor file = 1:noFiles(1) %change to for if you want to trouble shoot if it is flipping values for left and right, does not seem to flip if the whole group is not from the same side
    Fish = [];
    if ~iscellstr(filename(1,1))
        fullname = fullfile(pathname,filename{file});
        disp(fullname)
    else
    fullnamet = fullfile(pathname,filename(file));
    fullname = fullnamet{1};
    disp(fullname);
    end

        
    Data = dlmread(fullname);
    toRemove = sum(Data,2) == 0;
    Data(toRemove,:) = [];
    TPs = Data(:,7);
    delt = [0,diff(TPs)'];
    change = delt < -10000;
    Switch = find(change);
    offset = zeros(size(TPs));
    offset(Switch(1):end) = TPs(Switch(1)-1)+1;
    Data(:,7) = Data(:,7) + offset
    if strcmp(fullname((end-7):(end-4)),'left')
        ang = Data(:,3);
        Data(ang <= 180,3) = 180 - Data(ang <= 180,3);
        Data(ang > 180,3) = 180 + (360 - Data(ang>180,3));
        Data(:,1) = 900 - Data(:,1) + 330;%max(Data(:,1)) - Data(:,1) + min(Data((Data(:,1)>0),1));
        hi = file
    end
    Frames = unique(Data(:,7));
    noFrames = length(Frames);
    Fish(1).data = zeros(noFrames,3);
    Fish(2).data = zeros(noFrames,3);
    Fish(3).data = zeros(noFrames,3);
    Fish(4).data = zeros(noFrames,3);
    time = zeros(noFrames,1);
    v = 1:4;
    h = waitbar(0,sprintf('File %d progress', file));
    for f = 1:Fno
        Fish(f).data(1,:) = Data(f,1:3);
    end
    for t = 2:noFrames
        waitbar(t/noFrames,h)
        currFrame = Frames(t,:);
        frameData = Data(Data(:,7) == currFrame,:);
        time(t) = frameData(1,5);
        ys = frameData(:,2);
        for f = 1:Fno
            PossibleFish = ((ys < Lane(f,2)) & (ys > Lane(f,1)));
            if sum(PossibleFish) == 1
                Fish(f).data(t,:) = frameData(PossibleFish,1:3);
            else
                if sum(PossibleFish) == 0
                    Fish(f).data(t,:) = Fish(f).data(t-1,:);
                else % if there are more fish in one lane do this:
                    last_loc = vertcat(Fish(f).data(t-1,2));
                    curr_loc = frameData(PossibleFish,2);
                    [~,loc] = min(abs(last_loc-curr_loc));
                    PFNo = find(PossibleFish);
                    Fish(f).data(t,:) = frameData(PFNo(loc),1:3);
                end
            end
        end
                    
    
    end
    waitbar(1,h,sprintf('File %d done!!', file)); %moved out of loop so only says done when actually done
    RealTime = cumsum(time)/1000;
    for f = 1:Fno
        disp(f) %prints fish number

        noise = filtSTD(Fish(f).data(:,3),5);
        Bouting = [0,0,0,0,0,noise > 2] > 0;
        changeBout = [0,diff(Bouting)];
        Fish(f).BS = find(changeBout == 1);
        Fish(f).BE = find(changeBout == -1);
        
        if length(Fish(f).BS) > length(Fish(f).BE)
            Fish(f).BS(end) = [];
            
        end
        if ~isempty(Fish(f).BS)
            Fish(f).BS(end) = [];
            Fish(f).BE(end) = [];
            Fish(f).preBout_time = Fish(f).BS - 2;
            Fish(f).postBout_time = Fish(f).BE + 2;
            
            
            Fish(f).preBout_x = Fish(f).data(Fish(f).preBout_time,1);
            Fish(f).preBout_y = Fish(f).data(Fish(f).preBout_time,2);
            Fish(f).preBout_o = Fish(f).data(Fish(f).preBout_time,3);
            Fish(f).preBout_t = RealTime(Fish(f).preBout_time);
            Fish(f).postBout_x = Fish(f).data(Fish(f).postBout_time,1);
            Fish(f).postBout_y = Fish(f).data(Fish(f).postBout_time,2);
            Fish(f).postBout_o = Fish(f).data(Fish(f).postBout_time,3);
            Fish(f).postBout_t = RealTime(Fish(f).postBout_time);
            

            % new test 2025: Removal of mistrackings of the lane boarders.
            % cut 20 pixels of the borders (y axis) of all lanes to remove mistrackings of the
            % borders.
            % Initialize variables to store the corrected positions
            corrected_x = Fish(f).preBout_x; % Start with preBout_x values
            corrected_y = Fish(f).preBout_y; % Start with preBout_y values

            % Repeat the same process for post-bout coordinates
            corrected_post_x = Fish(f).postBout_x;
            corrected_post_y = Fish(f).postBout_y;

            % Define invalid y-ranges for each fish
            invalid_ranges = {
                [260, 280; 370, 390],        % Fish 1
                [390, 410; 499.25, 519.25],        % Fish 2
                [519.25, 539.25; 627.5, 647.5],  % Fish 3
                [647.5, 667.5; 757, 777]     % Fish 4
            };
            
            % Filter out preBout values that fall into invalid y-ranges
            keep_mask = true(size(corrected_y));
            for r = 1:size(invalid_ranges{f}, 1)
                range = invalid_ranges{f}(r, :);
                keep_mask = keep_mask & ~(corrected_y >= range(1) & corrected_y <= range(2));
            end

            %Remove based on X corrodinate threhsold: cut of 10 pixels of
            %the x axis. leaving the threshold at 340 to 880 (oiginal
            %length: 330-900)
            keep_mask=keep_mask & (corrected_x >= 340 & corrected_x <=890);
            keep_mask=keep_mask & (corrected_post_x >=340 & corrected_post_x <=890);

            % Apply the same mask to all pre-bout and post-bout fields
            corrected_x = corrected_x(keep_mask);
            corrected_y = corrected_y(keep_mask);
            Fish(f).preBout_time = Fish(f).preBout_time(keep_mask);
            Fish(f).preBout_x = corrected_x;
            Fish(f).preBout_y = corrected_y;
            
            corrected_post_x = corrected_post_x(keep_mask);
            corrected_post_y = corrected_post_y(keep_mask);
            Fish(f).postBout_time = Fish(f).postBout_time(keep_mask);
            Fish(f).postBout_x = corrected_post_x;
            Fish(f).postBout_y = corrected_post_y;


            % Removal of long distances   
            % Redefine x and y values after the removal of wall trackings.
            corrected_x = Fish(f).preBout_x; % Start with preBout_x values
            corrected_y = Fish(f).preBout_y; % Start with preBout_y values

            % Save the first valid coordinates as the initial reference point
            last_valid_x = corrected_x(1);
            last_valid_y = corrected_y(1);
            
            % Loop through all bout points starting from the second one
            for i = 2:length(corrected_x)
                % Calculate the Euclidean distance from the last valid position
                distance_to_last_valid = sqrt((corrected_x(i) - last_valid_x)^2 + (corrected_y(i) - last_valid_y)^2);
                
                % If the distance is greater than 150 pixels, reset the position to the last valid position
                if distance_to_last_valid > 150
                    corrected_x(i) = last_valid_x;
                    corrected_y(i) = last_valid_y;
                else
                    % Update the last valid position if within the 150-pixel threshold
                    last_valid_x = corrected_x(i);
                    last_valid_y = corrected_y(i);
                end
            end
            
            % Replace the original pre- and post-bout coordinates with the corrected values
            Fish(f).preBout_x = corrected_x;
            Fish(f).preBout_y = corrected_y;
            
            % Repeat the same process for post-bout coordinates
            corrected_post_x = Fish(f).postBout_x;
            corrected_post_y = Fish(f).postBout_y;
            
            % Reset the initial last valid positions for post-bout coordinates
            last_valid_x = corrected_post_x(1);
            last_valid_y = corrected_post_y(1);
            
            for i = 2:length(corrected_post_x)
                distance_to_last_valid = sqrt((corrected_post_x(i) - last_valid_x)^2 + (corrected_post_y(i) - last_valid_y)^2);
                
                if distance_to_last_valid > 150
                    corrected_post_x(i) = last_valid_x;
                    corrected_post_y(i) = last_valid_y;
                else
                    last_valid_x = corrected_post_x(i);
                    last_valid_y = corrected_post_y(i);
                end
            end
            
            % Save the corrected post-bout positions
            Fish(f).postBout_x = corrected_post_x;
            Fish(f).postBout_y = corrected_post_y;

            
            
            % Apply the filter to keep only valid bouts based on both within-bout and between-bout criteria
            %Fish(f).preBout_x = Fish(f).data(Fish(f).preBout_time,1);
            %Fish(f).preBout_y = Fish(f).data(Fish(f).preBout_time,2);
            Fish(f).preBout_o = Fish(f).data(Fish(f).preBout_time,3);
            Fish(f).preBout_t = RealTime(Fish(f).preBout_time);
            %Fish(f).postBout_x = Fish(f).data(Fish(f).postBout_time,1);
            %Fish(f).postBout_y = Fish(f).data(Fish(f).postBout_time,2);
            Fish(f).postBout_o = Fish(f).data(Fish(f).postBout_time,3);
            Fish(f).postBout_t = RealTime(Fish(f).postBout_time);


            % Recalculate delta after filtering to ensure the filtered bouts match expectations
            delta = Fish(f).postBout_o - Fish(f).preBout_o;
            wrongleft = delta > 180;
            wrongright = delta < -180;
            delta(wrongleft) = Fish(f).postBout_o(wrongleft) - (360 - Fish(f).preBout_o(wrongleft));
            delta(wrongright) = (360 + Fish(f).postBout_o(wrongright)) - Fish(f).preBout_o(wrongright);
            Fish(f).delta = delta;
            
            % Optional: Display remaining valid deltas for debugging
            %disp('Filtered deltas:');
            %disp(Fish(f).delta);
        end
    end
    AllFish{file} = Fish;
end



AllFish2 = [];
for i = 1:size(AllFish,2)
    AllFish2 = [AllFish2,AllFish{i}];
end
save(strcat(pathname,'\Results2.mat'),'AllFish2','-v7.3');
