function [PI,PI_mean,PI_SE,PI_binned,BadFish] = findPIs(AllFish2_filtered2,X)

    % value for pixel that do not change over time  
    %noMove = [];
        %for i = 1:size(AllFish2,2)
            %dist = (sqrt(sum(diff(AllFish2(i).data(:,[1,2])).^2,2)));
            %dist = smooth(dist,20);
            %noMove(i) = sum(dist == 0)/length(dist);

        %end
        %AllFish2(noMove > 0.7) = [];%can change, returned to original value but might want to use higher to make sure we are using fish that were succesfully trackedfor the majority of the time
        %BadFish = noMove > 1;%changed because does not matter for thresholding (not used for graphs)
        HalfwayPoint = X;
        for i = 1:size(AllFish2_filtered2,2)
            currFish = AllFish2_filtered2(i).data(:,1);
            SideOn = currFish > X;%changed from less than
            PI(i) = (sum(SideOn) - sum(~SideOn))/length(SideOn);
            totalT = size(currFish,1);
            ts = 1:totalT;
            bins = discretize(ts,linspace(1,totalT,10));
            for t = unique(bins)
                use = bins == t;
                SideOn = currFish(use) > X;
                PI_binned(i,t) = (sum(SideOn) - sum(~SideOn))/length(SideOn);
            end
                
                
        end

    PI_mean = mean(PI);
    %BadFish = noMove;
    PI_SE = std(PI)/sqrt(length(PI));
    PI;
end