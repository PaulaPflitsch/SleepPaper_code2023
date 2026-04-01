%% pick individual fish traces
 i=30; % change number between 1 and the total number of fish

% change fish number in brackets
x=AllFish2_filtered2(i).postBout_x; % check tracking of every single fish.
y=AllFish2_filtered2(i).postBout_y; % check tracking of every single fish
figure;
plot(x,y);
xlim([330 900]);
ylim([260 777]);

xlabel('X position');
ylabel('Y position');
title('X and y position in the swimming lane');


%% plot 4 individual fish traces
% plot multiple fish of the same group together
% plot with a color coded time series

figure;

%2.5mM Cadaverine, control
idxs = [25 2 18 30];

for i = 1:4
    subplot(4,1,i)

    x = AllFish2_filtered2(idxs(i)).postBout_x;
    y = AllFish2_filtered2(idxs(i)).postBout_y;
    t = 1:length(x);

    surface([x(:) x(:)], [y(:) y(:)], [t(:) t(:)], ...
        'EdgeColor','interp', ...
        'FaceColor','none', ...
        'LineWidth',1);

    title(['Fish ' num2str(idxs(i))])
    xlim([340,890]) %length of the lane (in pixels)
    colorbar;
end

colormap(turbo);
%colorbar;

%% Heatmap for individual fish
% Make a heatmap of the pixeled swimming lane showing where in the lane fish spend
% their time.

% group X pixel together for the heatmap.
% x values: lane length in pixels: 550 pixels (340-890) (after removing
% walls due to mistracings (in freeSwimAnalyze_astyanax_noJumps.m)
% y values: lane height in pixels: 109.25 pixels (after cutting)

% created after: https://de.mathworks.com/matlabcentral/answers/298669-how-to-bin-2d-data
for j = 1:size(AllFish2_filtered2, 2)
    %define x and y values to plot in heatmap
    x_hist = AllFish2_filtered2(j).postBout_x; 
    y_hist= AllFish2_filtered2(j).postBout_y;
    
    % create a new matrix for these values
    xy_hist=[x_hist, y_hist]
    %disp(xy_hist);
    sz=size(xy_hist);
    disp(sz);


    % pick bins:
    % 550px/~ 11 bins = 50px per x-bin
    % 109px/ ~5 bins = 22px per bin
    nBins_x=11;
    nBins_y=5;
    
    %make the binned figure for individual fish: one plot per fish
    
    [counts, centers]=hist3(xy_hist, [nBins_x, nBins_y]); %plots 3D bar chart
    
    figure('Units','centimeters','Position',[2 2 20 8]); %wide and short like lane shape
    
    imagesc(centers{1}, centers{2}, counts'); %transpose counts: hist3 output is (x,y), imagesc expects (row=y), col=x)
    axis xy; %flip y so origin is bottom left(like swim lane)
    colormap('hot');
    colorbar;
    
    xlabel('X position (px)');
    ylabel('Y position (px)');
    title('Fish position heatmap');

end;



%% Combined Heatmap across all fish, all lanes normalized

% Lane borders
lane_borders = [260, 389.25, 518.5, 647.75, 777];

% Pool all fish data
all_x = [];
all_y_norm = [];

for j = 1:size(AllFish2_filtered2, 2)
    x_fish = AllFish2_filtered2(j).postBout_x;
    y_fish = AllFish2_filtered2(j).postBout_y;

    % Remove NaNs
    valid = ~isnan(x_fish) & ~isnan(y_fish);
    x_fish = x_fish(valid);
    y_fish = y_fish(valid);

    % Assign lane by median y
    med_y = median(y_fish);
    lane_idx = find(med_y >= lane_borders(1:end-1) & med_y < lane_borders(2:end), 1);

    if isempty(lane_idx)
        warning('Fish %d: median y=%.1f outside all lanes, skipping.', j, med_y);
        continue;
    end

    % Normalize y to 0–129.25 using lane bottom border
    y_norm = y_fish - lane_borders(lane_idx);

    % Accumulate
    all_x = [all_x; x_fish];
    all_y_norm = [all_y_norm; y_norm];
end

fprintf('Total pooled data points: %d across %d fish\n', numel(all_x), size(AllFish2_filtered2, 2));

% Bin settings (same as individual plots)
nBins_x = 11;  % ~50px per bin along x (550px range)
nBins_y = 5;   % ~26px per bin along y (129.25px range)

xy_all = [all_x, all_y_norm];
[counts, centers] = hist3(xy_all, [nBins_x, nBins_y]);

% Plot
figure('Units', 'centimeters', 'Position', [2 2 20 8]);
imagesc(centers{1}, centers{2}, counts');
axis xy;
colormap('hot');
colorbar;

xlabel('X position (px)');
ylabel('Y position (normalized, 0–129.25 px)');
title(sprintf('Combined fish position heatmap — %d fish', size(AllFish2_filtered2, 2)));




