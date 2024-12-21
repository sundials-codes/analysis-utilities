function [A,b,c] = rkl_butcher_coeffs(rho,symbolic)
  % This function returns a Butcher table of Runge--Kutta coefficients corresponding
  % with a single step of the RKL2 method that would be stable for the spectral radius
  % rho.
  %
  % If "symbolic" is true, then the coefficients are given in symbolic variables
  % instead of double precision.

  % determine minimum number of stages to evaluate
  s = max(2,ceil((sqrt(9 + 8*abs(rho)) - 1)/2));

  if (symbolic)
    v = @(x) sym(x);
  else
    v = @(x) double(x);
  end

  % allocate Butcher tables (needs padding for initial state)
  A = zeros(s+1,s+1);
  b = zeros(1,s+1);

  % various constants
  one = v(1); two = v(2); three = v(3); four = v(4);
  w1 = four/((s + two)*(s - one));
  bjm2 = one/three;
  bjm1 = bjm2;

  % first stage is the old solution

  % second stage
  mus = w1*bjm1;
  cjm1 = mus;
  % tempv2 = yn;
  % tempv1 = yn + h*mus*fn;
  A(2,1) = mus;

  % remaining stages
  for j=2:s
    temj = (j + two)*(j - one);
    bj   = temj/(two*j*(j + one));
    ajm1 = one - bjm1;
    mu   = (two*j - one)/j*(bj/bjm1);
    nu   = -(j - one)/j*(bj/bjm2);
    mus  = w1*mu;
    cj   = temj*w1/four;

    %ycur = mu*tempv1 + nu*tempv2 + (one - mu - nu)*yn + h*mus*f(tempv1) - h*mus*ajm1*fn;
    A(j+1,1:j+1) = mu*A(j,1:j+1) + nu*A(j-1,1:j+1);
    A(j+1,j) = A(j+1,j) + mus;
    A(j+1,1) = A(j+1,1) - mus*ajm1;

    if(j < s)
      %tempv2 = tempv1;
      %tempv1 = ycur;
      cjm1 = cj;
      bjm2 = bjm1;
      bjm1 = bj;
    end
  end

  b = A(end,:);
  c = sum(A,2);

end
