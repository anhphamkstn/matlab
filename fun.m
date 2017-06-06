function [ z ] = fun( x )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    H = [1 -1; -1 2]; 
    f = [-2; -6];
    z = (1/2)*transpose(x)*H*x + transpose(f)*x;

end

