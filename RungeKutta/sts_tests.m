%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2024, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------
% driver to check analytical properties of the RKC and RKL
% method coefficients in Butcher form

clear; close all

% test parameters
plotregions = true;      % create stability region plots
runRKC = false;           % check RKC methods
runRKL = true;           % check RKL methods

% utility routines to determine rho for a given number of stages
rhoRKC = @(s) -0.653*s^2;
rhoRKL = @(s) - ((2*s + 1)^2 - 9)/8;

% RKC tests
if (runRKC)
  % numbers of stages to test, and plot region boxes to use
  tests = {
           {2,   [1.1*rhoRKC(2),   0.5,   -2,   2]},
           {4,   [1.1*rhoRKC(4),   0.5,   -4,   4]},
           {8,   [1.1*rhoRKC(8),   0.5,   -8,   8]},
           {16,  [1.1*rhoRKC(16),  0.5,  -16,  16]},
           {32,  [1.1*rhoRKC(32),  0.5,  -32,  32]},
           {64,  [1.1*rhoRKC(64),  0.5,  -84,  84]},
           {128, [1.1*rhoRKC(128), 0.5, -144, 144]},
     };

  fprintf('Runge--Kutta--Chebyshev methods\n');
  fprintf('        rho   |   s  |  q | lq\n');
  fprintf('   ----------------------------\n');
  for i = 1:length(tests)
    s = tests{i}{1};
    box = tests{i}{2};
    [A,b,c] = rkc_butcher_coeffs(rhoRKC(s), false);
    B = [c, A; 0, b];
    name = sprintf('RKC-%i', s);
    [q,p,qs,lq,lp,tol] = check_rk(B,0,plotregions,box,name,name);
    fprintf('    %9.2e | %4i | %2i | %2i\n', rhoRKC(s),s,q,lq);
  end
  fprintf('   ----------------------------\n');
  fprintf('\n');
end

% RKL tests
if (runRKL)
  % numbers of stages to test, and plot region boxes to use
  tests = {
           {2,   [1.1*rhoRKC(2),   0.5,   -2,   2]},
           {4,   [1.1*rhoRKC(4),   0.5,   -4,   4]},
           {8,   [1.1*rhoRKC(8),   0.5,   -8,   8]},
           {16,  [1.1*rhoRKC(16),  0.5,  -16,  16]},
           {32,  [1.1*rhoRKC(32),  0.5,  -32,  32]},
           {64,  [1.1*rhoRKC(64),  0.5,  -84,  84]},
           {128, [1.1*rhoRKC(128), 0.5, -144, 144]},
     };

  fprintf('Runge--Kutta--Legendre methods\n');
  fprintf('        rho   |   s  |  q | lq\n');
  fprintf('   ----------------------------\n');
  for i = 1:length(tests)
    s = tests{i}{1};
    box = tests{i}{2};
    [A,b,c] = rkl_butcher_coeffs(rhoRKL(s), false);
    B = [c, A; 0, b];
    name = sprintf('RKL-%i', s);
    [q,p,qs,lq,lp,tol] = check_rk(B,0,plotregions,box,name,name);
    fprintf('    %9.2e | %4i | %2i | %2i\n', rhoRKL(s),s,q,lq);
  end
  fprintf('   ----------------------------\n');
  fprintf('\n');
end

% end of script
