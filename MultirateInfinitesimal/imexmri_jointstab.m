function [xgrid,ygrid,Rmax] = imexmri_jointstab(G,W,dc,maxalpha,maxtheta,numRay,numAngle,box,numGrid,plottype)
  % Figure out the combined stability given A(theta) implicit stability region
  % and fast A(alpha) stability region
  % G cell array with gamma matrices
  % W cell array with omega matrices
  % dc differences between c values
  % maxalpha is the fast angle
  % box = [xl xr yb yt] defines x and y limits of box
  %       assume yt = -yb symmetry about the real axis
  % numRay is the number of points taken along each ray
  % numGrid is the number of points taken along each direction (x or y) in
  % the grid
  % numAngle is the number of angles sampled
  %
  % Rujeko Chinomona
  % Southern Methodist University
  % Department of Mathematics
  % Spring 2020
  %
  % Updated by Daniel Reynolds, SMU, Spring 2025

  coder.extrinsic('tic')
  coder.extrinsic('toc')

  % Get expressions for zf dependent functions
  [phivals,Mzf,Nzf] = zf_dependent_fns(G,W,dc,plottype);

  % Grid for explicit region sampling
  Nx = numGrid + 1;
  Ny = numGrid + 1;
  xl = box(1); xr = box(2); yb = box(3); yt = box(4);
  if abs(yb) ~= abs(yt)
    fprintf('Warning: yb should be -yt, code assumes symmetry about the real axis\n');
  end
  xgrid = linspace(xl,xr,Nx);
  ygrid = linspace(yb,yt,Ny);

  % Grid for fast region sampling
  Nrho =  numRay;
  Nalpha = numAngle;
  maxrho = 1;
  rho = -linspace(1e-6,maxrho,Nrho);
  alpha = linspace(-maxalpha*pi/180,maxalpha*pi/180,Nalpha);

  % Grid for implicit region sampling
  Nray = numRay;
  maxray = 5;
  Ntheta = numAngle;
  ray = -fliplr(logspace(-maxray,maxray,Nray));
  % ray = -linspace(1e-6,1,Nray);
  theta = linspace(-maxtheta*pi/180,maxtheta*pi/180,Ntheta);

  % define z functions
  switch plottype
  case 'explicit'
    ze_fun = @(x,y) xgrid(x) + 1i*ygrid(y);
    zi_fun = @(x,y) 0;
    zf_fun = @(x,y) 0;
  case 'implicit'
    ze_fun = @(x,y) 0;
    zi_fun = @(x,y) xgrid(x) + 1i*ygrid(y);
    zf_fun = @(x,y) 0;
  case 'imex'
    ze_fun = @(x,y) xgrid(x) + 1i*ygrid(y);
    zi_fun = @(x,y) ray(x)*(cos(theta(y)) + 1i*sin(theta(y)));
    zf_fun = @(x,y) 0;
  case 'explicit_mrigark'
    ze_fun = @(x,y) xgrid(x) + 1i*ygrid(y);
    zi_fun = @(x,y) 0;
    zf_fun = @(x,y) rho(x)*(cos(alpha(y)) + 1i*sin(alpha(y)));
  case 'implicit_mrigark'
    ze_fun = @(x,y) 0;
    zi_fun = @(x,y) xgrid(x) + 1i*ygrid(y);
    zf_fun = @(x,y) rho(x)*(cos(alpha(y)) + 1i*sin(alpha(y)));
  case 'imexmri'
    ze_fun = @(x,y) xgrid(x) + 1i*ygrid(y);
    zi_fun = @(x,y) ray(x)*(cos(theta(y)) + 1i*sin(theta(y)));
    zf_fun = @(x,y) rho(x)*(cos(alpha(y)) + 1i*sin(alpha(y)));
  otherwise
    warning('Unexpected plot type.')
  end

  totalcount = Nx*(Ny+1)/2;                   % Totalcount to gauge percent done
  Rthresh = 1.1;

  Rmax = zeros(Nx,Ny);
  count = 0;
  fprintf('Running %s with alpha = %4.2f, theta = %4.2f degrees.\n',plottype,maxalpha,maxtheta);
  tic
  for ix = 1:Nx
    for iy = 1:(Ny+1)/2
      ze = ze_fun(ix,iy);
      if isequal(plottype,'implicit') || isequal(plottype,'implicit_mrigark')
        zi = zi_fun(ix,iy);
      end
      for irho = 1:Nrho
        for ialpha = 1:Nalpha-1
          zf = zf_fun(irho,ialpha);
          for imray = 1:Nray
            for itheta = 1:Ntheta-1
              if ~isequal(plottype,'implicit') && ~isequal(plottype,'implicit_mrigark')
                zi = zi_fun(imray,itheta);
              end
              % keyboard
              Rval = stability_function(phivals,Mzf,Nzf,zf,zi,ze);
              Rmax(ix,iy) = max(Rmax(ix,iy),abs(Rval));
              Rmax(ix,Ny+1-iy) = Rmax(ix,iy);
            end
          end
        end
        if Rmax(ix,iy) > Rthresh
          break
        end
      end
      count = count + 1;
      if(mod(count,floor(totalcount/10))==0)
        percentdone = floor(count*100/(totalcount));
        fprintf('completed = %g, %g percent \n',count,percentdone);
        toc
      end
    end

  end
  toc

end



%%% internal utility routines %%%

function [p,M,N] = zf_dependent_fns(G,W,dc,plottype)
  % This function returns matlab anonymous functions
  % corresponding to:
  % p - matrix of phi values at each delta_c value:
  %     [phi_0(dc*zf), phi_1(dc*zf), ..., phi_K(dc*zf)]
  %     dimensions: s-by-(K+1) where s= num stages, K = num of coefficient matrices
  % M - matrix corresponding to sum_{k>=0} diag(phi_{k+1}(dc*zf))*G{k}
  % N - matrix corresponding to sum_{k>=0} diag(phi_{k+1}(dc*zf))*W{k}
  %
  % We compute these separately first to speed up the code
  % and avoid redundant computations

  % define symbolic variables
  syms zf t

  % definition of phi functions
  phi_0 = @(z) exp(z);
  phi_k = @(k,z)  int(exp(z.*(1-t)).*t.^(k-1),t,0,1);

  nummatrices = length(G);                   % number of coefficient matrices for G/W
  s = length(dc);                            % number of stages

  % initialize storage
  phivals = sym(zeros(s,nummatrices+1));
  Mzf = sym(zeros(s,s));
  Nzf = sym(zeros(s,s));

  switch plottype
  case {'implicit_mrigark','explicit_mrigark','imexmri'}
    % compute [phi_0(dc*zf), phi_1(dc*zf), ..., phi_K(dc*zf)]
    for i = 1:s
      phivals(i,1) = phi_0(dc(i)*zf);
      for j = 2:nummatrices + 1
        phivals(i,j) = phi_k(j-1,dc(i)*zf);
      end
    end

    % compute sum_{k>=0} diag(phi_{k+1}(dc*zf))*G{k} & sum_{k>=0} diag(phi_{k+1}(dc*zf))*W{k}
    for k = 1:nummatrices
      Nzf = Nzf + diag(phivals(:,k+1))*W{k};
      Mzf = Mzf + diag(phivals(:,k+1))*G{k};
    end

    % Convert to matlab anonymous functions
    M = matlabFunction(Mzf);
    N = matlabFunction(Nzf);
    p = matlabFunction(phivals);
    if (nargin(M) == 0)
      M = @(z) Mzf;
    end
    if (nargin(N) == 0)
      N = @(z) Nzf;
    end
    if (nargin(p) == 0)
      p = @(z) phivals;
    end

  case {'implicit','explicit','imex'}
    % No fast dependency, matrices are constant
    for i = 1:s
      phivals(i,1) = phi_0(0);
      for j = 2:nummatrices + 1
        phivals(i,j) = phi_k(j-1,0);
      end
    end

    for k = 1:nummatrices
      Nzf = Nzf + diag(phivals(:,k+1))*W{k};
      Mzf = Mzf + diag(phivals(:,k+1))*G{k};
    end

    M = @(z) double(Mzf);
    N = @(z) double(Nzf);
    p = @(z) double(phivals);
  otherwise
    warning('Unexpected plot type.')
  end
end



function Rval = stability_function(p,M,N,zf,zi,ze)
  % Computes R = (e_s^s)^T (I - diag(phi_0(dc*zf))*L - ze*N(zf) - zi*M(zf))^{-1}e_1

  % Compute zf dependencies first
  Mzf = M(zf);
  Nzf = N(zf);
  phivals = p(zf);
  [s,S] = size(Mzf);

  % (Faster) first way of computing
  Stage = zeros(s,1);
  Stage(1) = 1;
  for istage = 2:s
    implicit_prev = 0;
    explicit_prev = 0;
    for iprevious = 1:istage-1
      implicit_prev = implicit_prev + Mzf(istage,iprevious)*Stage(iprevious);
      explicit_prev = explicit_prev + Nzf(istage,iprevious)*Stage(iprevious);
    end
    st  = (phivals(istage,1)*Stage(istage-1) + zi*implicit_prev + ze*explicit_prev)/(1 - zi*Mzf(istage,istage));
    Stage(istage) = st;
  end
  Rval = st;

  % Second way of computing
  % I = eye(s);
  % L = diag(ones(s-1,1),-1);
  %
  % A = I - diag(phivals(:,1))*L - ze*Nzf - zi*Mzf;
  % b = phivals(1,1)*I(:,1);
  % Y = linsolve(A,b);
  % Rval = Y(end);

end
