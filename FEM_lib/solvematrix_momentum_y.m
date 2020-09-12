function [c] = solvematrix_momentum_y(~,S,f)
%SOLVEMATRIX Summary of this function goes here
%   Detailed explanation goes here

% solve the system for the non-fixed values
c = S\f;
end

