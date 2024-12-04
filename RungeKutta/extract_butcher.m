function [A,b,c,q] = extract_butcher(B)
% usage: [A,b,c,q] = extract_butcher(B)
%
% Utility routine to extract the Butcher table components from a Butcher table returned by butcher.m.
%
% Input:    B -- Butcher table
% Outputs:  A -- main coefficient matrix
%           b -- solution weights
%           c -- abscissae
%           q -- method order
%
% Daniel R. Reynolds
% SMU Mathematics
% December 2024

[rows,cols] = size(B);
s = cols-1;
A = B(1:s,2:s+1);
b = B(s+1,2:s+1);
c = B(1:s,1);
q = B(s+1,1);
