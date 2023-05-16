function [activity] = fct_crBinActivity(d,b,threshold, framerate)

% f = frame by frame data
% b = bin length in minutes
% threshold = minimum delta pixels to be considered a movement

% output = a matrix containing number of active seconds per minute
% this can be used to do the 'classic' measures 
% of waking activity and rest bouts


binLength = b; % minutes (typically 1 for wakeActivity, etc.)
frameRate = framerate; % Now deriving this from actual timestamp data!

binInterval = binLength * frameRate * 60; % number of frames in one Bin e.g. 1800frames/1min

% vector of #frames per bin until end experiment: 1, 1801, 3601, 5401 (why
% not start with zero?)
startBins = 1:binInterval:length(d); % frame number of first frame in each bin
startBins = startBins(1:end-1); % leave off last fragment
numSlices = length(startBins);  % number of bins (= minutes)

% matrix with zeros of dimensions #bins x #fish, size(d, 2)-> iw: 2 if
% nocibox (cr: gives 96 = #columns)CR-use 2 for sleepbox!; -> iw: 1 if
% sleepbox (#frames=#rows)
activity=zeros(numSlices,size(d,2)); 

for i = 1:length(startBins)
    % show progress report in command window
    fprintf('... doing bin %d of %d\n',i,numSlices);  
    
    % define slices/steps for analysis?
    sliver = d(startBins(i):startBins(i)+binInterval-1,1:size(d,2));
    
    % remove motion below a certain threshold of delta-pixels
    sliver(sliver<=threshold)=0;

    % determine number of frames that have any activity in this bin and sum
    % up for each bin (logical fct: converts numeric input into array of
    % logical values, nonzero of input is converted to logical 1 (true) and
    % zeros are converted to logical 0 (false), complex values and NaNs not
    % converted, will cause error, essentially make all movements = 1 and
    % sum ones up sum(A,1): sums each column (each fish); essentially # of
    % frames with movement per bin(minute)
    
    vals = sum(logical(sliver),1);

    % convert frames to seconds (30_frames in 1_sec = 4_frames in x_sec -->
    % 4/30 = sec)
    vals = vals / frameRate;
    
    % add data for this bin to the activity matrix (fill up matrix of zeros
    % bin by bin)
    activity(i,:) = vals;

end

end

% activity is active seconds per minute


