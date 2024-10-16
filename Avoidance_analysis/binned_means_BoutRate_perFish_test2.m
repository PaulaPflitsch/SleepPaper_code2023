function [AllFish2_filtered,bin_means_bout_rate_per_fish, mean_IBI_per_fish, bin_means_bout_rate, bin_means_IBI,IBIs_in_bin,bout_rates_in_bin,mean_IBI_per_fish_bin,mean_IBI_per_bin,mean_boutrate_per_bin,mean_boutrate_per_fish_bin,all_IBIs_with_postBout_x,all_bout_rates_with_postBout_x]  = binned_means_BoutRate_perFish_test2(AllFish2_filtered, num_bins)
    % Filter out fish with no movement
    %noMove = [];
    %AllFish2_filtered = AllFish2;  % Initialize the filtered variable
    %for i = 1:size(AllFish2_filtered, 2)
        %dist = sqrt(sum(diff(AllFish2_filtered(i).data(:, [1, 2])).^2, 2));
        %dist = smooth(dist, 20);
        %noMove(i) = sum(dist == 0) / length(dist);
    %end
    %AllFish2_filtered(noMove > 0.7) = [];  % Apply filtering to remove low movement fish
    
    IBIs_in_bin = [];
    bout_rates_in_bin = [];% Initialize storage for results
    bin_means_bout_rate_per_fish = cell(length(AllFish2_filtered), num_bins);
    mean_IBI_per_fish = cell(length(AllFish2_filtered), num_bins);
    mean_IBI_per_fish_bin=[];
    mean_boutrate_per_fish_bin=[];
    bin_means_bout_rate = zeros(1, num_bins);
    bin_means_IBI = zeros(1, num_bins);
    x_positions = 330:1:900;  % Define possible x positions (in pixels)
    
    % New matrices to track IBIs and bout rates with their postBout_x
    all_IBIs_with_postBout_x = [];  % To track IBIs and their corresponding postBout_x values
    all_bout_rates_with_postBout_x = [];  % To track bout rates and their corresponding postBout_x values

    % Define bin edges for 'num_bins' equal-width bins based on x positions
    bin_edges = linspace(min(x_positions), max(x_positions), num_bins + 1);

    % Loop through each bin
    for b = 1:num_bins
    %b=2
        %bin_IBIs_per_fish = cell(length(AllFish2_filtered), 1);  % IBIs per fish for the current bin
        %bin_bout_rates_per_fish = cell(length(AllFish2_filtered), 1);  % Bout rates per fish for the current bin
        
        IBIs_in_bin = [];  % Reset IBIs_in_bin for each bin to prevent overwriting issues
        bout_rates_in_bin = [];  % Reset bout_rates_in_bin for each bin

        % Now iterate through the remaining fish to calculate IBI
        for i = 1:size(AllFish2_filtered, 2)
        %for i=1 %small test run 
            % Extract preBout and postBout times for fish i
            %preBout = AllFish2(i).preBout;  % Assuming preBout is a field in AllFish2
            %postBout = AllFish2(i).postBout;  % Assuming postBout is a field in AllFish2
            
            % calculates the number of bouts that we will be using in the
            n_bouts=length(AllFish2_filtered(i).postBout_t) - 1;
            
            % Find bouts that fall within the current bin based on preBout_x
            bin_indices = AllFish2_filtered(i).preBout_x(1:n_bouts) >= bin_edges(b) & AllFish2_filtered(i).preBout_x(1:n_bouts) < bin_edges(b+1);
            
            if ~any(bin_indices)  % If no preBout_x values fall within the current bin, continue to next fish
                fprintf('No preBout_x values in bin %d for fish %d, skipping to next fish.\n', b, i);
                continue;  % Skip to the next fish
            end

            % Gate for bout removal to remove bouts that happened within
            % 50ms of each othe. Kristian code
            toRemove =  AllFish2_filtered(i).preBout_t(2:end)-AllFish2_filtered(i).postBout_t(1:end-1) < 0.05 %removes bouts that happen within 50ms of each other
            fprintf('fish i= %d, to remove: %d\n',i, toRemove);

            % Loop through postBout times, subtracting preBout(n+1) - postBout(n)
            for n = 1:size(AllFish2_filtered(i).postBout_t, 1)-1
            %for n=1:5  %small test run
                fprintf('Bin %d, Fish %d, n = %d: AllFish2_filtered(i).preBout_t(n+1) = %.3f, AllFish2_filtered(i).postBout_t(n) = %.3f\n', ...
                    b,i, n, AllFish2_filtered(i).preBout_t(n+1), AllFish2_filtered(i).postBout_t(n));
                % Skip if this bout was marked for removal
                if toRemove(n)
                    fprintf('Skipping bout at n = %d: preBout %.3f, postBout %.3f (within 50ms)\n', ...
                        n, AllFish2_filtered(i).preBout_t(n+1), AllFish2_filtered(i).postBout_t(n));
                    continue;  % Skip this iteration
                end
                
                if bin_indices(n)

                    % calculate IBI
                    IBI = AllFish2_filtered(i).preBout_t(n+1) - AllFish2_filtered(i).postBout_t(n);  % Calculate IBI (next preBout - current postBout)
                    fprintf('IBI= %d\n',IBI)
                    % Append the calculated IBI to the result matrix and keep
                    % track of the bin fish number (2nd col) and (3rd col) 
                    IBIs_in_bin = [IBIs_in_bin; IBI, i,b];  
                    
                    % Calculate bout rate as 1/IBI
                    bout_rate = 1 ./ IBI;
                    fprintf('bout rate = %d\n', bout_rate);
                    bout_rates_in_bin=[bout_rates_in_bin;bout_rate, i,b];
                else
                    continue
                end
                % Track all IBIs with corresponding postBout_x position
                all_IBIs_with_postBout_x = [all_IBIs_with_postBout_x; IBI, AllFish2_filtered(i).postBout_x(n), i];

                % Track all bout rates with corresponding postBout_x position
                all_bout_rates_with_postBout_x = [all_bout_rates_with_postBout_x; bout_rate, AllFish2_filtered(i).postBout_x(n), i];
            end  
        end

        % Calculate mean IBI and bout rate for this bin
        if ~isempty(IBIs_in_bin)  % If the matrix is not empty
            % Compute the mean IBI per fish for the current bin
            fish_combinations = IBIs_in_bin(:, 2);  % Extract fish indices (column 2)
            mean_IBI_per_fish_bin_b = accumarray(fish_combinations, IBIs_in_bin(:, 1), [], @mean);

            % bout rate mean per fish
            fish_combinations_br=bout_rates_in_bin(:,2); %bout rates
            mean_boutRate_per_fish_bin_b=accumarray(fish_combinations_br, bout_rates_in_bin(:, 1), [], @mean);

            % Append results into the final matrix, associating with the current bin
            for fish_idx = 1:length(mean_IBI_per_fish_bin_b)
                mean_IBI_per_fish_bin = [mean_IBI_per_fish_bin; mean_IBI_per_fish_bin_b(fish_idx), b];
            end

            for fish_idx_br = 1:length(mean_boutRate_per_fish_bin_b)
                mean_boutrate_per_fish_bin = [mean_boutrate_per_fish_bin; mean_boutRate_per_fish_bin_b(fish_idx_br), b];
            end
            %paste the mean bins per fish into a matrix that keeps track of
            %the bins as well (2nd column)
            %mean_IBI_per_fish_bin=[mean_IBI_per_fish_bin; mean_IBI_per_fish,b_column];
            
            % Calculate the mean IBI and bout rate for this fish inthis bin
            %bin_combinations=
        %else % skip empty IBI and bout rate matrices
            %mean_IBI_per_fish{i, b} = NaN;
            %bin_means_bout_rate_per_fish{i, b} = NaN;
        end

        
    end
    % Remove negative values from the IBI_matrix
     mean_IBI_per_fish_bin(mean_IBI_per_fish_bin(:, 1) <= 0, :) = [];  % Remove rows where IBI is negative
     mean_boutrate_per_fish_bin(mean_boutrate_per_fish_bin(:, 1) <= 0, :) = [];

     % Compute the mean IBI per fish for the current bin
    bin_combinations = mean_IBI_per_fish_bin(:, 2);  % Extract fish indices (column 2)
    mean_IBI_per_bin = accumarray(bin_combinations, mean_IBI_per_fish_bin(:, 1), [], @mean);

    bin_combinations = mean_boutrate_per_fish_bin(:, 2);  % Extract fish indices (column 2)
    mean_boutrate_per_bin = accumarray(bin_combinations, mean_boutrate_per_fish_bin(:, 1), [], @mean);
end