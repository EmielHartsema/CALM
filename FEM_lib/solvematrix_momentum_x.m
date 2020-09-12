function [c] = solvematrix_momentum_x(~,S,f)
%SOLVEMATRIX Summary of this function goes here
%   Detailed explanation goes here

%identify fixed value boundaries

c = S\f;
end

