function [Rmat,Rsym] = ark_stab_function(AE,AI,bE,bI)
% usage: [Rmat,Rsym] = ark_stab_function(AE,AI,bE,bI)
%
% Computes the stability function for an ARK method applied to the
% scalar, additive test problem,
%  
%    y'(t) = lI*y + lE*y,
%  
% using time step h, where Re(lI) < 0 and Re(lE) < 0.
% For an s-stage ARK method applied to this problem, then defining 
% the complex-valued inputs zI = h*lI and zE = h*lE, the function
% is given by
%
%   R(zE,zI) = 1 + (zE*bE + zI*bI)*((eye(s)-zE*AE-zI*AI)\e)
%
% where e is a vector of all ones in \R^s.
%
% Inputs:
%    AE, bE -- ERK Butcher table
%    AI, bI -- DIRK Butcher table
%
% Outputs:
%    Rmat -- Matlab function handle for stability function
%    Rsym -- symbolic expression for stability function
%
% Daniel R. Reynolds
% SMU Mathematics
% March 2020

% check Butcher tables for compatibility
[rowsE,colsE] = size(AE);
[rowsI,colsI] = size(AI);
if (length(bE) ~= length(bI))
    error('ark_stab_function: bI and bI are incompatible')
end
s = length(bE);
if ((s ~= rowsE) || (s ~= colsE))
    error('ark_stab_function: incompatible explicit Butcher table inputs')
end
if ((s ~= rowsI) || (s ~= colsI))
    error('ark_stab_function: incompatible implicit Butcher table inputs')
end

% ensure that bE and bI are row vectors
bE = reshape(bE,1,s);
bI = reshape(bI,1,s);

% create e, I
e = ones(s,1);
I = eye(s);

% create double-precision stability function handle
Rmat = @(zE,zI) 1 + (zE*double(bE) + zI*double(bI))*linsolve(I-zE*double(AE)-zI*double(AI),e);

% create symbolic stability function
syms zE zI;
Rsym = 1 + (zE*bE + zI*bI)*linsolve(I-zE*AE-zI*AI,e);

end
