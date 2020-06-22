function y = base_cal( x )
%EVL_FUC Summary of this function goes here
%   Detailed explanation goes here

pi = 3.1415926;
if (x < 0)
    error('error! input should be bigger than 0!');
% elseif (x == 0)
%     y = 1;
elseif (x < 10)
    y = exp(-0.4527 * x^0.86 + 0.0218);
else
    y = sqrt(pi/x) * (1 - 10/(7*x)) * exp(-x/4);
end
