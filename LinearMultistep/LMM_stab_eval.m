function [x,y] = LMM_stab_eval(theta,a,b)
% usage: [x,y] = LMM_stab_eval(theta,a,b)
%
% This function evaluates the boundary of the stability region for the LMM
% defined by the arrays a and b, for a given input angle theta.
%
% D.R. Reynolds
% Math 6321 @ SMU
% Fall 2016

% compute r
r = exp(i*theta);

% compute the number of steps for this LMM
n = length(a)-1;
m = length(b)-2;
p = max([n,m]);

% compute the stability region numerator
num = r^(p+1);
for j=0:n
   num = num - a(j+1)*r^(p-j);
end

% compute the stability region denominator
den = 0;
for j=-1:m
   den = den + b(j+2)*r^(p-j);
end

% compute the stability region boundary
z = num/den;

% extract the real and imaginary parts
x = real(z);
y = imag(z);

% end of function