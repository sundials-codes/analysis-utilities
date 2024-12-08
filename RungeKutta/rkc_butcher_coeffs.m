function [A,b,c] = rkc_butcher_coeffs(rho,symbolic)
  % This function returns a Butcher table of Runge--Kutta coefficients corresponding
  % with a single step of the RKC2 method that would be stable for the spectral radius
  % rho.
  %
  % If "symbolic" is true, then the coefficients are given in symbolic variables
  % instead of double precision.

  % determine minimum number of stages to evaluate
  s = ceil(sqrt(abs(rho)/0.653));

  if (symbolic)
    v = @(x) sym(x);
  else
    v = @(x) double(x);
  end

  % allocate Butcher tables (needs padding for initial state)
  A = zeros(s+1,s+1);
  b = zeros(1,s+1);

  % various constants
  one = v(1); two = v(2); four = v(4); c13 = v(13); zero = v(0);
  w0 = one + two/(c13*s^2);
  temp1 = w0^2 - one;
  temp2 = sqrt(temp1);
  arg = s*log(w0 + temp2);
  w1 = sinh(arg)*temp1 / (cosh(arg)*s*temp2 - w0*sinh(arg));
  bjm1 = one/(two*w0)^2;
  bjm2 = bjm1;
  mus = w1*bjm1;
  zjm1   = w0;
  zjm2   = one;
  dzjm1  = one;
  dzjm2  = zero;
  d2zjm1 = zero;
  d2zjm2 = zero;

  % first stage is the old solution

  % second stage
  %tempv2 = yn;
  %tempv1 = yn + h*mus*fn;
  A(2,1) = mus;

  % stages 3:s+1
  for j=3:s+1
    zj   =   two*w0*zjm1 - zjm2;
    dzj  =   two*w0*dzjm1 - dzjm2 + two*zjm1;
    d2zj =   two*w0*d2zjm1 - d2zjm2 + four*dzjm1;
    bj   =   d2zj/dzj^2;
    ajm1 =   one - zjm1*bjm1;
    mu   =   two*w0*bj/bjm1;
    nu   = - bj/bjm2;
    mus  =   mu*w1/w0;

    %ycur = mu*tempv1 + nu*tempv2 + (one - mu - nu)*yn + h*mus*f(tempv1) - h*mus*ajm1*fn;
    A(j,1:j) = mu*A(j-1,1:j) + nu*A(j-2,1:j);
    A(j,j-1) = A(j,j-1) + mus;
    A(j,1) = A(j,1) - mus*ajm1;

    if(j < s+1)
      %tempv2 = tempv1;
      %tempv1 = ycur;
      bjm2   = bjm1;
      bjm1   = bj;
      zjm2   = zjm1;
      zjm1   = zj;
      dzjm2  = dzjm1;
      dzjm1  = dzj;
      d2zjm2 = d2zjm1;
      d2zjm1 = d2zj;
    end
  end

  b = A(end,:);
  c = sum(A,2);

end
