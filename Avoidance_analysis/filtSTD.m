function [filtered] = filtSTD(signal,window)
    for t = 1:length(signal)
        if (t < (window+1)) || (t >(length(signal)-window))
            filtVar(t) = 0;
        else
            filtvar(t) = std(signal((t-window):t+window));
        end
    end
    filtered = filtvar;
end