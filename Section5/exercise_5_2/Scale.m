function [dataout] = Scale(datain,minval,maxval)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
dataout = datain - min(datain(:));
dataout = (dataout/range(dataout(:)))*(maxval-minval);
dataout = dataout + minval;

end

