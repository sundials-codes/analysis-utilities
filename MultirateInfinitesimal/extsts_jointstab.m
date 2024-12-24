function [xgrid,ygrid,Rmax] = extsts_jointstab(Ai,Ae,DMethod,maxtheta,maxrho,numRay,numAngle,box,numGrid)
  % Figure out the combined ExtendedSTS stability given A(theta) reaction stability
  % region and A(0) diffusion stability region
  % Ai -- slow implicit Butcher table (use 0 for explicit slow method)
  % Ae -- slow explicit Butcher table (use 0 for implicit slow method)
  % DMethod -- one of {'RKC', 'RKL'}, indicating the choice of super-time-stepping method
  % maxtheta -- maximum theta angle to test
  % maxrho -- maximum radius to test the A(theta) and A(0) regions
  % numRay - number of points along each ray emanating from the origin
  % numAngle -- number of angles sampled in in the [-maxtheta,maxtheta] sector
  % box = [xl xr yb yt] defines x and y limits of the joint stability plotting region
  %       assumes yt = -yb symmetry about the real axis
  % numGrid -- number of points taken along each direction (x or y) in the grid
  %
  % Daniel R. Reynolds
  % Southern Methodist University
  % Department of Mathematics
  % December 2024

  coder.extrinsic('tic')
  coder.extrinsic('toc')

  % Grid for explicit region sampling
  Nx = numGrid + 1;
  Ny = numGrid + 1;
  xl = box(1); xr = box(2); yb = box(3); yt = box(4);
  if abs(yb) ~= abs(yt)
    fprintf('Warning: yb should be -yt, code assumes symmetry about the real axis\n');
  end
  xgrid = linspace(xl,xr,Nx);
  ygrid = linspace(yb,yt,Ny);

  % Grid for diffusion region sampling
  Nrho =  numRay;
  Nalpha = 1;
  rho = -linspace(1e-6,maxrho,Nrho);
  alpha = zeros(1,Nalpha);

  % Grid for reaction region sampling
  Nray = numRay;
  maxray = 5;
  Ntheta = numAngle;
  ray = -fliplr(logspace(-maxray,maxray,Nray));
  % ray = -linspace(1e-6,1,Nray);
  theta = linspace(-maxtheta*pi/180,maxtheta*pi/180,Ntheta);

  % set plotting type based on Ai, Ae inputs
  if ((length(Ai) > 1) && (length(Ae) > 1))
    plottype = 'imex';
  elseif (length(Ai) == 1)
    plottype = 'explicit';
  else
    plottype = 'implicit';
  end

  % define z functions
  switch plottype
  case 'explicit'
    etaA_fun = @(x,y) xgrid(x) + 1i*ygrid(y);
    etaR_fun = @(x,y) 0;
    etaD_fun = @(x,y) rho(x)*(cos(alpha(y)) + 1i*sin(alpha(y)));
  case 'implicit'
    etaA_fun = @(x,y) 0;
    etaR_fun = @(x,y) xgrid(x) + 1i*ygrid(y);
    etaD_fun = @(x,y) rho(x)*(cos(alpha(y)) + 1i*sin(alpha(y)));
  case 'imex'
    etaA_fun = @(x,y) xgrid(x) + 1i*ygrid(y);
    etaR_fun = @(x,y) ray(x)*(cos(theta(y)) + 1i*sin(theta(y)));
    etaD_fun = @(x,y) rho(x)*(cos(alpha(y)) + 1i*sin(alpha(y)));
  end

  totalcount = Nx*(Ny+1)/2;                   % Totalcount to gauge percent done
  Rthresh = 1.05;

  Rmax = zeros(Nx,Ny);
  count = 0;
  fprintf('Running %s for %s with theta = %4.2f degrees.  Nx = %i, Ny = %i, Nrho = %i, Nalpha = %i\n', plottype, DMethod, maxtheta, Nx, Ny, Nrho, Nalpha);
  tic
  for ix = 1:Nx
    for iy = 1:(Ny+1)/2
      etaA = etaA_fun(ix,iy);
      if isequal(plottype,'implicit')
        etaR = etaR_fun(ix,iy);
      end
      for irho = 1:Nrho
        for ialpha = 1:Nalpha
          etaD = etaD_fun(irho,ialpha);
          for imray = 1:Nray
            for itheta = 1:Ntheta-1
              if ~isequal(plottype,'implicit')
                etaR = etaR_fun(imray,itheta);
              end
              Rval = extsts_stability_function(Ai, Ae, DMethod, etaR, etaA, etaD);
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
        fprintf('  completed %g of %g (%g%%). ', count, totalcount, percentdone);
        toc
      end
    end

  end

  % output some statistics regarding Rmax
  fprintf("Rmax: max = %g, min = %g, avg = %g. ", max(max(Rmax)), min(min(Rmax)), mean(Rmax, 'all'));
  toc; fprintf("\n");

end
