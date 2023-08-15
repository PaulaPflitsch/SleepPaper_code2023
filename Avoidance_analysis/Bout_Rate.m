% bout rate for left and right sides
function [total_boutrate, count_l, count_r, bpm_left, bpm_right,boutrateDiff,total_bout_right,total_bout_left,timeinterbout_left,timeinterbout_right,timebout_left,timebout_right,sum_timebout_left,sum_timeinterbout_left,sum_timeinterbout_right,sum_timebout_right,time_left,time_right,boutrate_left,boutrate_right] = Bout_Rate(AllFish2,d_thresh)
    noMove = [];
    for i = 1:size(AllFish2,2)
        dist = (sqrt(sum(diff(AllFish2(i).data(:,[1,2])).^2,2)));
        dist = smooth(dist,20);
        noMove(i) = sum(dist == 0)/length(dist);
    end
    AllFish2(noMove > 0.7) = [];%testing instead of .45

    for i = 1:size(AllFish2,2)
       count_l=0;
       count_r=0;
       bpm_right=0,
       bpm_left=0;
       
        for n = 1:size(AllFish2(i).postBout_x,1)
            %disp(n)
            total_boutrate(i)=n;             % calculate the total boutrate per fish
            bpm(i)=total_boutrate(i)/30;     % calculate the bouts per minute (30 min trials?!) --> check with Kristian
           
            %convert pixels to cm to distinguish between left and right
            %9.5 cm/570 pixels, 60 pixels/cm
            %450 or less by agarose, 780 or more by cadaverine (2 cm)
            if AllFish2(i).preBout_x(n)<(330+60*d_thresh) % d_thresh is the distance in centimeters we want to look at
                timebout_left(i,n) = AllFish2(i).postBout_t(n)-AllFish2(i).preBout_t(n);    % bout length
                count_l=count_l+1;
                total_bout_left(i)=count_l; %bout rate left 
                bpm_left= total_bout_left/30;
                for m=2:size(AllFish2(i).postBout_t,1)   %substract the second value from the first
                    timeinterbout_left(i,n)=AllFish2(i).preBout_t(m)-AllFish2(i).postBout_t(n);
                    %disp(timeinterbout_left)
                end
                sum_timebout_left=sum(timebout_left,2);
                
                sum_timeinterbout_left=sum(timeinterbout_left,2);
            end
            
            if AllFish2(i).preBout_x(n)>(900 - 60*d_thresh)
                timebout_right(i,n) = AllFish2(i).postBout_t(n)-AllFish2(i).preBout_t(n);    % bout length
                count_r=count_r+1
                total_bout_right(i)=count_r
                bpm_right= total_bout_right/30
                for m=2:size(AllFish2(i).postBout_t,1)   %substract the second value from the firsts
                    timeinterbout_right(i,n)=AllFish2(i).preBout_t(m)-AllFish2(i).postBout_t(n);
                    %disp(timeinterbout_right)
                end
                sum_timebout_right=sum(timebout_right,2);
                
                sum_timeinterbout_right=sum(timeinterbout_right,2);
                
            end
        end
        %%both lines give the same result for time left
        %time_left(i)=sum_timebout_left(i)+sum_timeinterbout_left(i);
        time_left=nonzeros(sum_timebout_left)+nonzeros(sum_timeinterbout_left);
        time_left=time_left/60000; % from ms to minutes

        time_right=nonzeros(sum_timebout_right)+nonzeros(sum_timeinterbout_right);
        time_right=time_right/60000; % from ms to minutes

        total_bout_left=nonzeros(total_bout_left);
        total_bout_right=nonzeros(total_bout_right);

        boutrate_left=total_bout_left./time_left;
        boutrate_left = boutrate_left(~isinf(boutrate_left));
        boutrate_right=total_bout_right./time_right;
        boutrate_right = boutrate_right(~isinf(boutrate_right));
       
    end

    for f= 1:size(bpm_left,2)
        % calculate the bout rate difference between left and right
        boutrateDiff(f) = bpm_right(f) - bpm_left(f);

    end

    % exclude bout rates that are obviously unrealistic
    %for u=1:size boutrate_right

       
       
      
        
end





