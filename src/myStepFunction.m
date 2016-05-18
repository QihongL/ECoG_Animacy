function [ output ] = myStepFunction( input )
% step function: R^n -> R^n 
% f(x_i) =  1, if x_i >= 0 
%           0, if x_i < 0 
% for all i 

output = nan(length(input),1);
output(input >= 0) = 1;  % the point zero has measure zero anyway
output(input < 0) = 0;

end

