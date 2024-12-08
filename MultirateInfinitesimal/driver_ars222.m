function driver_ars222(maxalpha)

  addpath('../RungeKutta')

  % create base Butcher tables (including padding and stiff accuracy)
  zed = 0;  % zed = sym(0);
  one = 1;  % one = sym(1);
  two = 2;  % two = sym(2);
  three = 3;  % three = sym(3);
  gamma = (two-sqrt(two))/two;
  delta = one-one/(two*gamma);
  c = [zed; gamma; gamma; one; one];
  be = [delta, zed, one-delta, zed, zed];
  de = [zed, zed, two/three, zed, one/three];
  Ae = [zed, zed, zed, zed, zed;
        gamma, zed, zed, zed, zed;
        gamma, zed, zed, zed, zed;
        delta, zed, one-delta, zed, zed;
        delta, zed, one-delta, zed, zed];
  Be = [c, Ae; 2, be; 1, de];
  bi = [zed, zed, one-gamma, zed, gamma];
  di = de;
  Ai = [zed, zed, zed, zed, zed;
        gamma, zed, zed, zed, zed;
        zed, zed, gamma, zed, zed;
        zed, zed, one, zed, zed;
        zed, zed, one-gamma, zed, gamma];
  Bi = [c, Ai; 2, bi; 1, di];

  % verify properties of ARK table
  disp('checking ARK method properties')
  check_ark_embedded(c,c,Ae,Ai,be,bi,de,di,1e-11,1,false,0,0,0);

  % convert Butcher tables to MRI "Gamma" and "Omega" matrices
  [cmri, Wmri] = mis_to_mri(Be);
  if (cmri ~= c)
    error('cmri does not match c for explicit MIS table')
  end
  [cmri, Gmri] = mis_to_mri(Bi);
  if (cmri ~= c)
    error('cmri does not match c for implicit MIS table')
  end

  % pad W and G with an initial row of zeros to match expected table structure
  W{1} = [zeros(1,5); Wmri{1}];
  G{1} = [zeros(1,5); Gmri{1}];

  % set "dc" increment array (pad with initial 0)
  dc = [0; c(2:end)-c(1:end-1)];

  % test parameters
  thetavals = [10,30,45,60,80,90];
  numRay = 25;
  numGrid = 200;
  numAngle = 2;
  box = [-3,0.5,-3,3];
  header = {};
  plottype = 'imexmri';
  plotcolors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4940 0.1840 0.5560], ...
                [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
  plotlinestyle = {'-','--','-.',':','-','--'};

  filename = ['imexmri_ars222_alpha_',num2str(maxalpha),'.mat'];
  q = matfile(filename,'Writable',true);
  q.box = box;
  q.maxalpha = maxalpha;
  q.thetavals = thetavals;
  q.numRay = numRay;
  q.numAngle = numAngle;
  q.numGrid = numGrid;
  figure
  hold on

  for itheta = 1:length(thetavals)
    maxtheta = thetavals(itheta);
    [xgrid,ygrid,Rmax] = imexmri_jointstab(G,W,dc,maxalpha,maxtheta,numRay,numAngle,box,numGrid,plottype);
    R{itheta} = Rmax;
    contour(xgrid, ygrid, Rmax', [1+eps,1+eps], 'color', plotcolors{itheta}, 'LineStyle', ...
            plotlinestyle{itheta}, 'LineWidth', 2);
    header{itheta} = [num2str(maxtheta),char(176)];
    hold on
  end

  q.R = R;
  q.xgrid = xgrid;
  q.ygrid = ygrid;

  plottype = 'imex';
  maxtheta = 90;
  numRay = 50;
  numGrid = 50;
  [xgrid,ygrid,Rmax] = imexmri_jointstab(G,W,dc,maxalpha,maxtheta,numRay,numAngle,box,numGrid,plottype);
  contour(xgrid,ygrid,Rmax',[1+eps,1+eps],'k:','LineWidth',2);
  header{itheta+1} = ['Base (',num2str(maxtheta),char(176),')'];
  q.Rbase = Rmax;
  q.xbase = xgrid;
  q.ybase = ygrid;


  hold off
  title(['\alpha = ',num2str(maxalpha), char(176)]);
  lgd = legend(header);
  lgd.Location = 'best';
  lgd.Title.String = '\theta values';
  axis square

  plotname = ['imexmri_ars222_alpha_',num2str(maxalpha)];
  print('-dpng',plotname);
  % savefig(plotname);
