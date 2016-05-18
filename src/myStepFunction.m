function [ output ] = myStepFunction( input )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

output = nan(length(input),1);
output(input >= 0) = 1;  % the point zero has measure zero anyway
output(input < 0) = 0;

end

