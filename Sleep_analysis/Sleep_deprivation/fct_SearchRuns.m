function [RestBoutNumber, RestBoutLength ] = fct_SearchRuns(a)
%%Found this method to search for number runs online:
%http://stackoverflow.com/questions/18282458/how-to-count-length-of-runs-of-values-in-matrix-columns 


a_log = logical(a);
a_log_inverse = ~a_log;

nrows = size(a_log_inverse,1);
a_log_inverse_neg_cell = num2cell(~a_log_inverse,[nrows 1]);
zeros_a_log_inverse = cellfun(@find, a_log_inverse_neg_cell, 'UniformOutput', 0);
find_runs = @(v) nonzeros( diff([0; v; nrows+1])-1 ).';

rest_bouts = cellfun(find_runs, zeros_a_log_inverse, 'UniformOutput', 0);

%rest_bouts = cell array of matrixes. Number elements in each matrix = #rest
%bouts for that fish,


RestBoutNumber = cellfun(@numel,rest_bouts); %Number of rest bouts for each fish


RestBoutLength = cellfun(@mean, rest_bouts); %Mean length for each fish;

end

