function Rval = extsts_stability_function(Ai, Ae, DMethod, etaR, etaA, etaD)
  % Evaluates the ExtSTS stability function for a given Ai, Ae, DMethod, etaR, etaA, and etaD

  % handle the case where one of Ai or Ae is zero
  if ((length(Ai) > 1) && (length(Ae) > 1))
    s = size(Ae,1);
    c = Ai*ones(s,1);
    ce = Ae*ones(s,1);
    if (c ~= ce)
      error('Error: Ai and Ai are not internally consistent')
    end
  end
  if (length(Ai) == 1)
    Ai = zeros(size(Ae));
    s = size(Ae,1);
    c = Ae*ones(s,1);
  end
  if (length(Ae) == 1)
    Ae = zeros(size(Ai));
    s = size(Ai,1);
    c = Ai*ones(s,1);
  end

  % set the Delta c and e1, es arrays
  dc = [0; c(2:end)-c(1:end-1)];
  e1 = [1; zeros(s-1,1)];
  es = [zeros(s-1,1); 1];

  % generate arrays of the fast stability function components for each stage
  Rrk = zeros(s,1);
  Rf = zeros(s,1);
  for i=2:s
    [Rrk(i), Rf(i)] = fast_stability(DMethod, dc(i)*etaD);
  end

  % generate fast propagation matrices
  P = diag(Rrk(2:s), -1);
  Q = diag(Rf) - diag(Rf(2:s),-1);

  Rval = es' * ((eye(s) - P - Q*(etaR*Ai + etaA*Ae))\e1);

end

%%% internal utility routine %%%

function [Rrk, Rf] = fast_stability(DMethod,etaD)
  % This function evaluates the second-order super-time-stepping propagation
  % function factors over a scaled step of size etaD,
  %    v_i = Rrk(etaD)*v_{i-1} + Rf(etaD)*g_i
  % returning both Rrk(etaD) and Rf(etaD) as outputs.
  %
  % Here, DMethod must be either 'RKC' or 'RKL' to specify the type of
  % super-time-stepping method.

  addpath('../RungeKutta')

  % get the Butcher table for the STS method
  if (DMethod == 'RKC')
    [A,b,c] = rkc_butcher_coeffs(etaD,false);
  end
  if (DMethod == 'RKL')
    [A,b,c] = rkl_butcher_coeffs(etaD,false);
  end

  % set vector of ones and inner matrix
  one = ones(size(c));
  M = eye(length(c)) - etaD*A;

  % evaluate propagation function components
  Rrk = 1 + etaD*b*(M\one);
  Rf = 1 + etaD*b*(M\c);

end
