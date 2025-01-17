function [xgrid,ygrid,Rmax] = extsts_jointstab(Ai,Ae,DMethod,RxParams,DiffParams,box,numGrid)
  % Figure out the combined ExtendedSTS stability given A(theta) reaction stability
  % region and A(0) diffusion stability region
  % Ai -- slow implicit Butcher table (use 0 for explicit slow method)
  % Ae -- slow explicit Butcher table (use 0 for implicit slow method)
  % DMethod -- one of {'RKC', 'RKL'}, indicating the choice of super-time-stepping method
  % RxParams -- reaction sector parameters [max angle to test, num angles to test, max radius to test, num samples per ray]
  % DiffParams -- diffusion sampling parameters [max radius to test, num samples per ray]
  % box = [xl xr yb yt] defines x and y limits of the joint stability plotting region
  %       assumes yt = -yb symmetry about the real axis
  % numGrid -- number of points taken along each direction (x or y) in the grid (should be an even number)
  %
  % Daniel R. Reynolds
  % Southern Methodist University
  % Department of Mathematics
  % December 2024

  coder.extrinsic('tic')
  coder.extrinsic('toc')

  % ensure that numGrid is even (otherwise, increment by one)
  if (mod(numGrid,2) == 1)
    fprintf('Warning: numGrid should be even, incrementing by one\n');
    numGrid = numGrid + 1;
  end

  % Sampling grid to use in generating the plot
  Nx = numGrid + 1;
  Ny = numGrid + 1;
  xl = box(1); xr = box(2); yb = box(3); yt = box(4);
  if abs(yb) ~= abs(yt)
    fprintf('Warning: yb should be -yt, code assumes symmetry about the real axis\n');
    yb = -max(abs(yb), abs(yt));
    yt =  max(abs(yb), abs(yt));
  end
  xgrid = linspace(xl, xr, Nx);
  ygrid = linspace(yb, yt, Ny);

  % Diffusion region sampling values
  if (abs(DiffParams(1)) <= 0)
    fprintf('Warning: max diffusion radius must be positive, resetting to default');
    DiffParams(1) = 100;
  end
  if (abs(DiffParams(2)) < 2)
    fprintf('Warning: num diffusion samples must be strictly greater than 1, resetting to default');
    DiffParams(2) = 2;
  end
  DRadii = -linspace(1e-6, DiffParams(1), DiffParams(2));

  % Reaction region sampling values
  if (abs(RxParams(1)) < 0)
    fprintf('Warning: max reaction angle must be non-negative, resetting to default');
    RxParams(1) = 0;
  end
  if (abs(RxParams(2)) < 1)
    fprintf('Warning: num reaction angles must be at least 1, resetting to default');
    RxParams(2) = 1;
  end
  if (abs(RxParams(3)) <= 0)
    fprintf('Warning: max reaction radius must be strictly greater than 0, resetting to default');
    RxParams(3) = 1;
  end
  if (abs(RxParams(4)) < 1)
    fprintf('Warning: num reaction samples per ray must be at least 1, resetting to default');
    RxParams(4) = 1;
  end
  RAngles = linspace(0, RxParams(1)*pi/180, RxParams(2));
  RRadii = -linspace(1e-6, RxParams(3), RxParams(4));

  % set plotting type based on Ai, Ae inputs, and set sampling functions
  if ((length(Ai) > 1) && (length(Ae) > 1))
    plottype = 'imex';
    etaA_fun = @(ix,iy) xgrid(ix) + 1i*ygrid(iy);
    etaR_fun = @(ix,iy) RRadii(ix)*(cos(RAngles(iy)) + 1i*sin(RAngles(iy)));
  elseif (length(Ai) == 1)
    plottype = 'explicit';
    etaA_fun = @(ix,iy) xgrid(ix) + 1i*ygrid(iy);
    etaR_fun = @(ix,iy) 0;
    RAngles = [0];
    RRadii = [0];
  else
    plottype = 'implicit';
    etaA_fun = @(ix,iy) 0;
    etaR_fun = @(ix,iy) xgrid(ix) + 1i*ygrid(iy);
  end


  totalcount = Nx*(Ny+1)/2;                   % Totalcount to gauge percent done
  Rthresh = 1.05;

  Rmax = zeros(Nx,Ny);
  count = 0;
  fprintf('Running %s for %s with max(RAngle) = %4.2f degrees and max(DiffRho) = %4.2f.  Nx = %i, Ny = %i, NDiff = %i, NRx = %i\n', plottype, DMethod, RAngles(end)*180/pi, DiffParams(1), Nx, Ny, length(DRadii), length(RAngles)*length(RRadii));
  tic
  for ix = 1:Nx  % loop over mesh
    for iy = 1:(Ny+1)/2

      % get sampling coordinate
      etaA = etaA_fun(ix,iy);
      if isequal(plottype,'implicit')
        etaR = etaR_fun(ix,iy);
      end

      % loop over diffusion region
      for idiff = 1:length(DRadii)

        % get diffusion coordinate
        etaD = DRadii(idiff);

        % loop over reaction region, evaluating stability function, and storing values symmetrically above/below real axis
        for irxradius = 1:length(RRadii)
          for irxangle = 1:length(RAngles)
            if ~isequal(plottype,'implicit')
              etaR = etaR_fun(irxradius,irxangle);
            end
            Rval = extsts_stability_function(Ai, Ae, DMethod, etaR, etaA, etaD);
            Rmax(ix,iy) = max(Rmax(ix,iy),abs(Rval));
            Rmax(ix,Ny+1-iy) = Rmax(ix,iy);
          end
        end

        % break out of loop if we've violated stability
        if Rmax(ix,iy) > Rthresh
          break
        end
      end

      % periodically output overall progress
      count = count + 1;
      if ((mod(count,ceil(totalcount/10))==0) && (count ~= totalcount))
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
