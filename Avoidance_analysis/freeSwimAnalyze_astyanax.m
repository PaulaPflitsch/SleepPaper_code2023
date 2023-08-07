clear variables %added to clear variables from previous runs, make sure there is no impact from left over variables
[filename,pathname] = uigetfile('*.*','C:\Users\Engert\Dropbox\FreeSwimming','MultiSelect','on');
if iscell(filename)
    noFiles = size(filename,2);
else
    noFiles = 1;
end
AllFish = {[]};
Fno = 4;
top = 777;
bottom = 260;
for i = 1:4
    Lane(i,:) = [(bottom + (i-1)*(top-bottom)/4),(260 + i*(top-bottom)/4)];
end
parfor file = 1:noFiles(1)%change to for if you want to trouble shoot if it is flipping values for left and right, does not seem to flip if the whole group is not from the same side
    Fish = [];
    if ~iscellstr(filename(1,1))
        fullname = fullfile(pathname,filename);
    else
    fullnamet = fullfile(pathname,filename(file));
    fullname = fullnamet{1};
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
            delta = Fish(f).postBout_o - Fish(f).preBout_o;
            wrongleft = delta > 180;
            wrongright = delta < -180;
            delta(wrongleft) = Fish(f).postBout_o(wrongleft) - (360 - Fish(f).preBout_o(wrongleft));
            delta(wrongright) = (360 + Fish(f).postBout_o(wrongright)) - Fish(f).preBout_o(wrongright);
            Fish(f).delta = delta;
        end
    end
    AllFish{file} = Fish;
end



AllFish2 = [];
for i = 1:size(AllFish,2)
    AllFish2 = [AllFish2,AllFish{i}];
end
save(strcat(pathname,'\Results2.mat'),'AllFish2','-v7.3');
