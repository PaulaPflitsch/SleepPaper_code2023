function [fish] = fct_crSleepPlots_comparison(ctrl,treatment,geno,expName)

%% d typically comes from binActivity
% from binned data in 1 min intervals

% geno is a cell array of genotypes and individuals of each genotype
% for example geno = {'wt','mut';1:48,49:96};

% make 'classic' charts for
%   Activity (sec / 10 min)
%   Sleep (min / 10 min)
%   Waking Activity (secs / waking minute)

% remove NaN if present: isnan finds NaNs in a and marks them as logical 1
% (true), other as logical 0 (false), in this case NaNs replaced with zero
%a(isnan(a))=0; 
binSize = 10; % in minutes 
numBins_ctrl = floor(size(ctrl,1) / binSize); % floor rounds to smaller value
numBins_treat = floor(size(treatment,1) / binSize);

% make 10 min bins vector
tRange_ctrl = 1:binSize:numBins_ctrl*binSize; 
tRange_treat = 1:binSize:numBins_treat*binSize; 

%convert to binary minutes = fraction of day (MATLAB minutes?)
tRange_ctrl = tRange_ctrl .* 1/(24*60);
tRange_treat = tRange_treat .* 1/(24*60);
 
% first, extract Activity, Rest, and Waking Activity from ctrl and
% treatment and put in bins of appropriate size

i=1;
% make empty vectors for activity, rest, wakeAct
activity_ctrl = zeros(numBins_ctrl,size(ctrl,2));
rest_ctrl = zeros(numBins_ctrl,size(ctrl,2));
wakeAct_ctrl = zeros(numBins_ctrl,size(ctrl,2));

activity_treat = zeros(numBins_treat,size(treatment,2));
rest_treat = zeros(numBins_treat,size(treatment,2));
wakeAct_treat = zeros(numBins_treat,size(treatment,2));

for b = 1:binSize:numBins_ctrl*binSize
    theseRows = b:b+binSize-1;
    sliver_ctrl = ctrl(theseRows,:);
    
    % get Activity data (add the number of seconds in 10 minutes)
    activity_ctrl(i,:) = sum(sliver_ctrl);
    
    % get Rest data (add the number of minutes with zero in 10 min)
    rest_ctrl(i,:) = sum(logical(~sliver_ctrl));
    
    % get Waking Activity 
    % (number of seconds in 10 min / number of minutes with activity) 
    wakeAct_ctrl(i,:) = sum(sliver_ctrl) ./ sum(logical(sliver_ctrl));  
    
    i=i+1;
end

for c = 1:1:binSize:numBins_treat*binSize
    theseRows = c:c+binSize-1;
    sliver_treat = treatment(theseRows,:);
    
    % get Activity data (add the number of seconds in 10 minutes)
    activity_treat(i,:) = sum(sliver_treat);
    
    % get Rest data (add the number of minutes with zero in 10 min)
    rest_treat(i,:) = sum(logical(~sliver_treat));
    
    % get Waking Activity 
    % (number of seconds in 10 min / number of minutes with activity)  
    wakeAct_treat(i,:) = sum(sliver_treat) ./ sum(logical(sliver_treat)); 
    
    i=i+1;
end

% wakeAct has some NaN's because there are some minutes with no activity
% set these to zero
wakeAct_ctrl(isnan(wakeAct_ctrl))=0;
wakeAct_treat(isnan(wakeAct_treat))=0;
 
% now for the plots
figh = figure('Position',[5 300 1440 350],'Color',[1 1 1]);

wakeh = axes('Position',[0.07 0.1 0.25 0.8]); hold on;
resth = axes('Position',[0.39 0.1 0.25 0.8]); hold on;
acth = axes('Position',[0.71 0.1 0.25 0.8]); hold on;
 
colorvals = [0 0 0; ... %  gray; 0 0 0 black
    0.8 0 0;...       % purple; 1 0 1magenta; 0.8 0 0 red
    0.5 0.5 0.5];     % magenta? google colorspec or color values
% green 0 0.5 0.5
% yellow  1 1 0
% gray 0.5 0.5 0.5

for f = 1:size(geno,2) % 
    
    fish(f).genotype = geno{1,f};
    fish(f).color = colorvals(f,:);
    fish(f).fish = geno{2,f};
    
    fish(f).activity_ctrl = activity_ctrl(:,fish(f).fish);
    fish(f).rest_ctrl = rest_ctrl(:,fish(f).fish);
    fish(f).wakeAct_ctrl = wakeAct_ctrl(:,fish(f).fish);
   
    cdark = colorvals(f,:);
    
    % Activity plot:
    plotme(acth,fish(f).activity_ctrl,tRange,cdark,cdark);
    
    % Rest Plot
    plotme(resth,fish(f).rest_ctrl,tRange,cdark,cdark);
    
    % Waking Activity Plot
    plotme(wakeh,fish(f).wakeAct_ctrl,tRange,cdark,cdark);
   
end

for f = 1:size(geno,2) % 
    
    fish(f).genotype = geno{1,f};
    fish(f).color = colorvals(f,:);
    fish(f).fish = geno{2,f};
    
    fish(f).activity_treat = activity_treat(:,fish(f).fish);
    fish(f).rest_treat = rest_treat(:,fish(f).fish);
    fish(f).wakeAct_treat = wakeAct_treat(:,fish(f).fish);
   
    cdark = colorvals(f,:);
    
    % Activity plot:
    plotme(acth,fish(f).activity_treat,tRange,cdark,cdark);
    
    % Rest Plot
    plotme(resth,fish(f).rest_treat,tRange,cdark,cdark);
    
    % Waking Activity Plot
    plotme(wakeh,fish(f).wakeAct_treat,tRange,cdark,cdark);
end
ylabel(acth,'Activity (sec / 10 min)');
ylabel(resth,'Rest (bouts / 10 min)');
ylabel(wakeh,'Waking Activity (secs / waking minute)');
fct_suptitle([' Sleep-Wake Plots: ' expName]);

end

function plotme(h,plotdata,tRange,clight,cdark)

    axes(h);
    ytoplot = mean(plotdata,2);
    err = std(plotdata,0,2); % std
    err = err / sqrt(length(plotdata)); % SEM   
    [x y] = fct_makeRibbon(tRange,ytoplot-err,ytoplot+err);
    fill(x,y,clight,'FaceColor',clight,'EdgeColor',clight); 
    plot(tRange,ytoplot,'color',cdark,'lineWidth',1);
    datetick(h,'x','HH');
    set(h,'xlim',[tRange(1) tRange(end)]);
    set(h,'FontSize',16);
    
end


