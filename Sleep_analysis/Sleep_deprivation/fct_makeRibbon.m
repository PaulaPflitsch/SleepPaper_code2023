function [outx,outy] = fct_makeRibbon(x,bottom,top)
% use this function to make "ribbon" style polygons
% ex. a ribbon where the top and bottom are given by the mean +/- SD
% expects inputs to be a single dimension. 
% use fill(outx,outy,'color') to set the color
%
% no secret here, just trial and error to figure out how MATLAB's
% fill/patch commands want their inputs formatted to get sensible output.
%
% example: [x y] = makeRibbon(timestamps,mean-sd,mean+sd);
%          fill(x,y,'r','facecolor','r','edgecolor','r')
%          will yield a red ribbon at points timestamps that goes between
%          the mean-sd and mean+sd
x = x(:);
bottom = bottom(:);
top = top(:);

outx = [x;flipud(x)];
outy = [bottom;flipud(top)];