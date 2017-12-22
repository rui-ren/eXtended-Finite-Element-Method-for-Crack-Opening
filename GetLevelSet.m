function [ff] = GetLevelSet(a, b, c, xx, yy)
%Get the level-set function 
%Define discontinuity
% b*y = a*x +c
ff = (a * xx + b * yy + c) / (sqrt(a^2 + b^2));
end

