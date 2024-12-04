function [A,b,c,d,q,p] = extract_butcher_embedded(B)
% usage: [A,b,c,d,q,p] = extract_butcher_embedded(B)
%
% Utility routine to extract the embedded Butcher table
% components from a Butcher table returned by butcher.m.
%
% Input:    B -- Butcher table
% Outputs:  A -- main coefficient matrix
%           b -- solution weights
%           c -- abscissae
%           d -- embedding weights
%           q -- method order
%           p -- embedding order
%
% Daniel R. Reynolds
% SMU Mathematics
% December 2024

[rows,cols] = size(B);
if (rows <= cols)
  error('Input Butcher table does not contain an embedding; call extract_butcher instead')
end
s = cols-1;
A = B(1:s,2:s+1);
b = B(s+1,2:s+1);
c = B(1:s,1);
d = B(s+2,2:s+1);
q = B(s+1,1);
p = B(s+2,1);
