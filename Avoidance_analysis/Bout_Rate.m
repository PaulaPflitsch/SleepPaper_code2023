% bout rate for left and right sides
function [total_boutrate, count_l, count_r, bpm_left, bpm_right,total_bout_left,total_bout_right,timeinterbout_left,timeinterbout_right,timebout_left,timebout_right,sum_timebout_left,sum_timeinterbout_left,sum_timeinterbout_right,sum_timebout_right,time_left,time_right,AllFish2_filtered,boutrate_left,boutrate_right,boutrateDiff,boutrate_left_for_fish_i] = Bout_Rate(AllFish2_filtered,d_thresh)
    noMove = [];
    for i = 1:size(AllFish2_filtered,2)
        dist = (sqrt(sum(diff(AllFish2_filtered(i).data(:,[1,2])).^2,2)));
        dist = smooth(dist,20);
        noMove(i) = sum(dist == 0)/length(dist);
    end
    AllFish2_filtered(noMove > 0.7) = [];%testing instead of .45

    %just for 0mM,1mM remove fish that never swim within 2cm on the left side
    % 1mM control: [11,26,27,41,44]
    % 1mM sleep: [2,15,19]
    % 0mM control: [13,16,17,25,26,28]
    %AllFish2_filtered([2,15,19]) = [];

    for i = 1:size(AllFish2_filtered,2)
       count_l=0;
       count_r=0;
       bpm_right=0,
       bpm_left=0;
       
        for n = 1:size(AllFish2_filtered(i).postBout_x,1)
            %disp(n)
            total_boutrate(i)=n;             % calculate the total bouts per fish
            bpm(i)=total_boutrate(i)/30;     % calculate the bouts per minute (30 min trials?!) --> check with Kristian
           
            %convert pixels to cm to distinguish between left and right
            %9.5 cm/570 pixels, 60 pixels/cm
            %450 or less by agarose, 780 or more by cadaverine (2 cm)
            if AllFish2_filtered(i).preBout_x(n)<(330+60*d_thresh) % d_thresh is the distance in centimeters we want to look at
                count_l=count_l+1;
                total_bout_left(i)=count_l; %bout rate left 
                bpm_left= total_bout_left/30;
                timebout_left(i,n) = AllFish2_filtered(i).postBout_t(n)-AllFish2_filtered(i).preBout_t(n);    % bout length
                for m=2:size(AllFish2_filtered(i).postBout_t,1)   %substract the second value from the first
                    timeinterbout_left(i,n)=AllFish2_filtered(i).preBout_t(m)-AllFish2_filtered(i).postBout_t(n);
                    %disp(timeinterbout_left)
                end
                sum_timebout_left=sum(timebout_left,2);
                sum_timeinterbout_left=sum(timeinterbout_left,2);
            end
          
            
            if AllFish2_filtered(i).preBout_x(n)>(900 - 60*d_thresh)
                count_r=count_r+1
                total_bout_right(i)=count_r
                bpm_right= total_bout_right/30
                timebout_right(i,n) = AllFish2_filtered(i).postBout_t(n)-AllFish2_filtered(i).preBout_t(n);    % bout length
                for m=2:size(AllFish2_filtered(i).postBout_t,1)   %substract the second value from the firsts
                    timeinterbout_right(i,n)=AllFish2_filtered(i).preBout_t(m)-AllFish2_filtered(i).postBout_t(n);
                    %disp(timeinterbout_right)
                end
                sum_timebout_right=sum(timebout_right,2);
                sum_timeinterbout_right=sum(timeinterbout_right,2);
                
            end
        end
        % NOT use this part
        %%both lines give the same result for time left
        %time_left(i)=sum_timebout_left(i)+sum_timeinterbout_left(i);
        %time_left=nonzeros(sum_timebout_left)+nonzeros(sum_timeinterbout_left);
        %time_left=time_left/60000; % from ms to minutes

        %time_right=nonzeros(sum_timebout_right)+nonzeros(sum_timeinterbout_right);
        %time_right=time_right/60000; % from ms to minutes

        %total_bout_left=nonzeros(total_bout_left);
        %total_bout_right=nonzeros(total_bout_right);

        %boutrate_left=total_bout_left./time_left;
        %boutrate_left = boutrate_left(~isinf(boutrate_left));
        %boutrate_right=total_bout_right./time_right;
        %boutrate_right = boutrate_right(~isinf(boutrate_right));


    end 
        % use this part
        time_left=sum_timebout_left+sum_timeinterbout_left;
        time_left=time_left/60000;
        
        time_right=sum_timebout_right+sum_timeinterbout_right;
        time_right=time_right/60000;

        total_bout_left=total_bout_left(:); % turn row into column
        total_bout_right=total_bout_right(:);

        boutrate_left=total_bout_left./time_left;
        boutrate_right=total_bout_right./time_right;
        
        %replace nans with 0
        boutrate_left( isnan(boutrate_left) ) = 0;
        boutrate_right( isnan(boutrate_right) ) = 0; 
        disp(boutrate_left);
        disp(boutrate_right);
        
        boutrate_left_for_fish_i = boutrate_left(i);
        boutrate_right_for_fish_i = boutrate_right(i);
      

    % calculate the bout rate Difference within 2cm of each side
        boutrateDiff=minus(boutrate_left,boutrate_right);
    
    % exclude bout rates that are higher than 200 because they are
    % unrealistic
    rowToDelete_left=boutrate_left>200;
    boutrate_left(rowToDelete_left)=[];
    rowsToDelete_right=boutrate_right>200;
    boutrate_right(rowsToDelete_right)=[];
    rowsToDelete = boutrateDiff>200 | boutrateDiff<-200;
    boutrateDiff(rowsToDelete) = [];
    
       
       
      
        
end





