function [C] = ssp_coefficient(A,b)
  % Usage: [C] = ssp_coefficient(A,b)
  %
  % Computes the SSP coefficient for a given RK method.  This uses the
  % approach outlined in https://doi.org/10.1007/s10915-017-0560-2,
  % Section 2.
  %
  % Inputs:
  %    A -- Butcher table matrix
  %    b -- Butcher table solution coefficients
  %
  % Outputs:
  %    C -- SSP coefficient
  %
  %------------------------------------------------------------
  % Programmer(s):  Daniel R. Reynolds @ SMU
  %------------------------------------------------------------
  % Copyright (c) 2024, Southern Methodist University.
  % All rights reserved.
  % For details, see the LICENSE file.
  %------------------------------------------------------------

  % disable warnings in this function (due to tests with singular matrices)
  warning('off', 'all')

  % ensure that A and b are compatible dimensions
  stages = length(b);
  if (size(A,1) ~= stages)
    error('the number of rows in A does not match the length of b');
  end
  if (size(A,2) ~= stages)
    error('the number of columns in A does not match the length of b');
  end

  % construct square matrix S that contains both A and b
  S = zeros(stages+1, stages+1);
  S(1:stages,1:stages) = A;
  S(stages+1,1:stages) = reshape(b, 1, stages);

  % construct identity matrix I, and vector of ones
  I = eye(stages+1);
  e = ones(stages+1, 1);

  % create utility functions
  Re = @(r) (I + r*S)\e;
  P = @(r) r*((I + r*S)\S);
  ReMin = @(r) min(Re(r));
  PMin = @(r) min(min(P(r)));
  Test = @(r) ((ReMin(r) >= 0) && (PMin(r) >= 0));

  % find bounds on SSP coefficient
  Cvals = logspace(-15,200);
  Cfound = false;
  for i = 1:length(Cvals)-1
    Clo = Cvals(i);
    Chi = Cvals(i+1);
    if Test(Clo) && ~Test(Chi)
      Cfound = true;
      break
    end
  end
  if ~Cfound
    error('cannot find lower/upper bounds on SSP coefficient')
  end

  % perform bisection method to hone in on SSP coefficient
  tol = 1e-12;
  while (Chi - Clo) / 2 > tol
    C = (Chi + Clo) / 2;
    if Test(C)
      Clo = C;
    else
      Chi = C;
    end
  end

end

% re-enable warnings
warning('on', 'all')

% end of function
