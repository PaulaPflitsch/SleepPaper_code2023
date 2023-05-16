function [fish] = fct_crSleepPlots(a,geno,expName)

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
a(isnan(a))=0; 
binSize = 10; % in minutes 
numBins = floor(size(a,1) / binSize); % floor rounds to smaller value

% make 10 min bins vector
tRange = 1:binSize:numBins*binSize; 

%convert to binary minutes = fraction of day (MATLAB minutes?)
tRange = tRange .* 1/(24*60); 
 
% first, extract Activity, Rest, and Waking Activity from a
% and put in bins of appropriate size

i=1;
% make empty vectors for activity, rest, wakeAct
activity = zeros(numBins,size(a,2));
rest = zeros(numBins,size(a,2));
wakeAct = zeros(numBins,size(a,2));

for b = 1:binSize:numBins*binSize
    theseRows = b:b+binSize-1;
    
    sliver = a(theseRows,:);
    
    % get Activity data (add the number of seconds in 10 minutes)
    activity(i,:) = sum(sliver);
    
    % get Rest data (add the number of minutes with zero in 10 min)
    rest(i,:) = sum(logical(~sliver));
    
    % get Waking Activity 
    % (number of seconds in 10 min / number of minutes with activity) 
    wakeAct(i,:) = sum(sliver) ./ sum(logical(sliver));  
    
    i=i+1;
end

% wakeAct has some NaN's because there are some minutes with no activity
% set these to zero
wakeAct(isnan(wakeAct))=0;
 
% now for the plots
figh = figure('Position',[5 300 1440 350],'Color',[1 1 1]);

wakeh = axes('Position',[0.07 0.1 0.25 0.8]); hold on;
resth = axes('Position',[0.39 0.1 0.25 0.8]); hold on;
acth = axes('Position',[0.71 0.1 0.25 0.8]); hold on;

% Assume geno order as ctrl, hom, het
colorvals = [0 0 0; ... % black
             0.8 0 0;... % red
             0.5 0.5 0.5; % gray
             0.4, 0, 0.8]; %purple
colorvals_face = [0.4 0.4 0.4; ... % dark gray
             0.8 0.45 0.45;... % pink
             0.7 0.7 0.7; % light gray
             0.75, 0.65, 0.8]; %light purple
% 0.9290, 0.6940, 0.1250; %yellow
% green 0 0.5 0.5
% yellow  1 1 0
% gray 0.5 0.5 0.5
% orange 0.8500, 0.3250, 0.0980

for f = 1:size(geno,2) % 
    
    fish(f).genotype = geno{1,f};
    fish(f).color = colorvals(f,:);
    fish(f).fish = geno{2,f};
    
    fish(f).activity = activity(:,fish(f).fish);
    fish(f).rest = rest(:,fish(f).fish);
    fish(f).wakeAct = wakeAct(:,fish(f).fish);
   
    cdark = colorvals(f,:);
    clight = colorvals_face(f,:);
    
    % Activity plot:
    plotme(acth,fish(f).activity,tRange,clight,cdark);
    
    % Rest Plot
    plotme(resth,fish(f).rest,tRange,clight,cdark);
    
    % Waking Activity Plot
    plotme(wakeh,fish(f).wakeAct,tRange,clight,cdark);
   
end

ylabel(acth,'Activity (sec / 10 min)');
ylabel(resth,'Rest (min / 10 min)');
ylabel(wakeh,'Waking Activity (secs / minute)');
fct_suptitle([' Sleep-Wake Plots: ' expName]);

end

function plotme(h,plotdata,tRange,clight,cdark)

    axes(h);
    ytoplot = mean(plotdata,2);
    err = std(plotdata,0,2); % std
    err = err / sqrt(length(plotdata)); % SEM   
    [x, y] = fct_makeRibbon(tRange,ytoplot-err,ytoplot+err);
    fill(x,y,clight,'FaceColor',clight,'EdgeColor',clight); 
    plot(tRange,ytoplot,'color',cdark,'lineWidth',0.5);
    %datetick(h,'x','MM');
    datetick(h,'x','HH');
    set(h,'xlim',[tRange(1) tRange(end)]);
    set(h,'FontSize',14);
    
end


