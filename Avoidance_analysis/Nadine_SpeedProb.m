function [speedProb,towardsCad,awayCad, locProb] = Nadine_SpeedProb(AllFish2)
    noMove = [];
        for i = 1:size(AllFish2,2)
            dist = (sqrt(sum(diff(AllFish2(i).data(:,[1,2])).^2,2)));
            dist = smooth(dist,20);
            noMove(i) = sum(dist == 0)/length(dist);
        end
        AllFish2(noMove > 0.6) = [];
        %9.5 cm/570 pixels, 60 pixels/cm distance = 330+i*60
        for i = 1:9
            totSpeed = 0;
            count = 0;
            speedCad = 0;
            C = 0;
            speedAway = 0;
            A = 0;
            for j = 1:size(AllFish2,2)
                for n = 1:size(AllFish2(j).postBout_x,1)
                    if AllFish2(j).preBout_x(n) > (330 + (i-1)*60) && AllFish2(j).preBout_x(n) < (330 + i*60)
                        distanceCM = (AllFish2(j).postBout_x(n)-AllFish2(j).preBout_x(n))/60.0;
                        timeCM = AllFish2(j).postBout_t(n)-AllFish2(j).preBout_t(n);
                        speed = distanceCM/timeCM;
                        if speed >0
                            speedCad = speedCad + speed;
                            C = C+1;
                        end
                        if speed < 0
                            speedAway = speedAway + abs(speed);
                            A = A+1;
                        end
                        
                        totSpeed = totSpeed + abs(speed);
                        count = count + 1;
                    end
                end
            end
            speedProb(i) = totSpeed/count;
            locProb(i) = count;
            towardsCad(i) = speedCad/C;
            awayCad(i) = speedAway/A;
        end
end
