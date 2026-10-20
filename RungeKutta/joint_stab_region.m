function joint_stab_region(AE,AI,bE,bI,thetas,box,fig,rho)
% Usage: joint_stab_region(AE,AI,bE,bI,thetas,box,fig,rho)
%
% Computes the joint stability region for an ARK method applied
% to the scalar, additive test problem,
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
% For a given angle 0 <= theta <= pi/2 and radius rho > 0, we define
% the joint stability region as
%
%   J(theta,rho) = { zE \in \C : |R(zE,zI)| <= 1 forall zI in S(theta,rho) },
%
% where S(theta,rho) is the sector
%
%   S(theta,rho) = { zI = -a+i*b : |zI| <= rho, a >= 0, and atan(|b|/a) <= theta }.
%
% Since AE is strictly lower triangular and AI is lower triangular with
% nonnegative diagonal, R(zE,zI) is analytic in zI throughout S(theta,rho),
% so by the maximum modulus principle we need only test the boundary of
% S(theta,rho): the two rays zI = r*exp(i*(pi +/- theta)), 0 <= r <= rho,
% and the arc zI = rho*exp(i*phi), pi-theta <= phi <= pi+theta.  We sample
% each ray at NI points with logarithmically-scaled distance from the
% origin (along with zI=0), and the arc at NA uniformly-spaced angles.
% We test an NE x NE uniform mesh of points zE within the pre-defined
% "box" in the complex plane.  Since all coefficients are real,
% J(theta,rho) is symmetric across the real axis, so we only evaluate
% zE in the upper half-plane, and mirror the results.
%
% The input 'thetas' is array-valued -- we plot the joint stability
% region J(theta,rho) for each value in this array, and overlay these plots.
%
%
% Inputs:
%    AE, bE -- ERK Butcher table
%    AI, bI -- DIRK Butcher table
%    thetas -- array of sector angles (in degrees) to use in creating
%              overlaid plots
%    box    -- [xl, xr, yl, yr] is the bounding box for the sub-region
%              of the complex plane in which to perform the test.  We
%              assume that yl=-yr.
%    fig    -- figure handle to use
%    rho    -- (optional) sector radius [default 1e8]
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2024, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% set general parameters
if (nargin < 8)
  rho = 1e8;
end
NE = 201;  % must be odd
NI = 100;
NA = 50;
Rthresh = 1.25;
CM = lines(length(thetas));

% retrieve the general ARK stability function
[R,~] = ark_stab_function(AE,AI,bE,bI);

% set mesh of ERK sample points
xl = box(1);
xr = box(2);
yl = box(3);
yr = box(4);
x = linspace(xl,xr,NE);
y = linspace(yl,yr,NE);

% create new figure window
xlim = box(1:2);  ylim = box(3:4);
xax = plot(linspace(xlim(1),xlim(2),10),zeros(1,10),'k:'); hold on
yax = plot(zeros(1,10),linspace(ylim(1),ylim(2),10),'k:');

% reusable value for progress output
Ntot = NE*(NE-1)/2;

% loop over theta values, creating contour plot data for each
for itheta = 1:length(thetas)
  theta = thetas(itheta)*pi/180;  % convert to radians
  %fprintf('\n joint_stab_region: theta = %3.0f (%i of %i):', theta*180/pi, itheta, length(thetas));

  % initialize max|R| over box
  Rmax = zeros(NE,NE);

  % set array of DIRK sample points: both sector rays (log-spaced in |zI|),
  % and the outer arc of the sector
  r = logspace(-4,log10(rho),NI);
  phi = linspace(pi-theta,pi+theta,NA);
  zI = [0, r*exp(1i*(pi-theta)), r*exp(1i*(pi+theta)), rho*exp(1i*phi)];

  % loop over zE mesh
  iTot = 0;
  for j=0:(NE-1)/2
    for i=1:NE

      % output intermittent progress
      iTot = iTot+1;
      %if (mod(iTot,floor(Ntot/10)) == 0)
      %  fprintf(' %i%%', floor(100*iTot/Ntot));
      %end

      % set zE value
      zE = x(i) + y((NE+1)/2+j)*sqrt(-1);

      % loop over zI values, breaking the moment |R(zE,zI)| > Rthresh
      for k=1:length(zI)
        Rval = abs(R(zE,zI(k)));
        j1 = (NE+1)/2+j;
        j2 = (NE+1)/2-j;
        Rmax(j1,i) = max(Rmax(j1,i),Rval);
        Rmax(j2,i) = max(Rmax(j2,i),Rval);
        if (Rval > Rthresh)
          break;
        end
      end
    end
  end

  % create contours and add to figure
  %c = contourc(x,y,Rmax,[1 1]);
  c = contourc(x,y,Rmax,[1+eps 1+eps]);

  % assemble plot
  figure(fig);
  hold on
  idx = 1;
  while idx < size(c,2)
    cols = idx+1:idx+c(2,idx);
    X = c(1,cols);
    Y = c(2,cols);
    lstring = ['$\theta =\;$', sprintf('%g', thetas(itheta))];
    plot(X, Y, 'color', CM(itheta,:), 'DisplayName', lstring, 'LineWidth', 2)
    idx = idx + c(2,idx) + 1;
  end
  hold off

end

% finish up figure
set(get(get(xax,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
set(get(get(yax,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
axis(box);
xlabel('Re(zE)');
ylabel('Im(zE)');
%legend(arrayfun(@(th) sprintf('theta = %g', th), thetas, 'UniformOutput', false),'Location', 'northwest')
legend('Location', 'northwest', 'Interpreter', 'latex');
hold off


% end of function
