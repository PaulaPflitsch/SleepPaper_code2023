% this function dispalys the sleep and waking times over experiment time in
% a binary (black and white) plot
function [a_binary_matrix] = fct_Binary_activity(a_clip,allWells,expName)

% converts a_clip into a binary activity matrix
% Sleep is defined as 1 minute or more of no movement. Bins in a_clip have 
% a length of 1 minute, thus, if a bin contains a
% number other than 0,the fish were not asleep. In turn, bins with '0'
% signal that the fish was asleep.
% By converting the matrix into a binary matrix we can plot sleep vs
% non-sleep.
% a black bar indicated at elast 1 minute of sleep

a_binary = a_clip>0;

%plot the activity in a heatmap
%make different groups based on genotype data
% allWells
% plot each fish individually in heatmap
figure('Color',[1 1 1]);
colormap('bone'); % default bone gray
ax=imagesc(a_binary(:,allWells)');
% set(ax,'box','off');
% colorbar('location','southoutside');
xlabel('time in seconds');
ylabel('Individual Fish Sleep (non-active seconds/min)');
title('Sleep "heatmap"');
fct_suptitle([' Overview Sleep allWells: ' expName]);

a_binary_matrix = a_binary;

end