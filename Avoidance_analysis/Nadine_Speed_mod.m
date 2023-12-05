function [lSE, rSE, fishLeft, fishRight, speedFishDiff, meanSpeed, meanLeft, meanRight, totalRight, count, totalSpeed, AllFish2_filtered, speedArray] = Nadine_Speed_mod(AllFish2, d_thresh)
    noMove = [];
    for i = 1:size(AllFish2, 2)
        dist = sqrt(sum(diff(AllFish2(i).data(:, [1,2])).^2, 2));
        dist = smooth(dist, 20);
        noMove(i) = sum(dist == 0) / length(dist);
    end
    AllFish2(noMove > 0.7) = [];

    totalSpeed = 0;
    l = 0;
    r = 0;
    count = 0;
    totalLeft = 0;
    totalRight = 0;
    speedLC = 0;
    numFish = size(AllFish2, 2);
    speedArray = zeros(numFish, 1);
    AllFish2_filtered = AllFish2;  % Initialize the filtered variable

    for i = 1:size(AllFish2, 2)
        fishSpeedR = 0;
        fishSpeedL = 0;
        lc = 0;
        rc = 0;
        for n = 1:size(AllFish2(i).postBout_x, 1)
            distanceCM = sqrt((AllFish2(i).postBout_x(n) - AllFish2(i).preBout_x(n)).^2 + (AllFish2(i).postBout_y(n) - AllFish2(i).preBout_y(n)).^2) / 60.0;
            timeCM = AllFish2(i).postBout_t(n) - AllFish2(i).preBout_t(n);
            speed = distanceCM / timeCM;
            speedArray(i)=speed; % save speed in an array

            % Exclude rows where speed > 6
            if speed < 6
                %450 or less by agarose, 780 or more by cadaverine (2 cm)
                %d_thresh is 2
                %510, 720 (3mc)
                if AllFish2(i).preBout_x(n) < (330 + 60 * d_thresh)
                    left(l + 1) = speed;
                    totalLeft = totalLeft + speed;
                    l = l + 1;
                    lc = lc + 1;
                    fishSpeedL = fishSpeedL + speed;
                end
                if AllFish2(i).preBout_x(n) > (900 - 60 * d_thresh)
                    totalRight = totalRight + speed;
                    right(r + 1) = speed;
                    r = r + 1;
                    rc = rc + 1;
                    fishSpeedR = fishSpeedR + speed;
                end
                % Store the filtered data
                AllFish2_filtered(i).speed(n) = speed;
            end
        end
        fishLeft(i) = fishSpeedL / lc;
        fishRight(i) = fishSpeedR / rc;
        speedFishDiff(i) = fishRight(i) - fishLeft(i);
    end

    meanSpeed = totalSpeed / count;
    meanLeft = totalLeft / l;
    meanRight = totalRight / r;
    lSE = std(fishLeft) / sqrt(length(fishLeft));
    rSE = std(fishRight) / sqrt(length(fishRight));
end