function driver_ssp_dirk(maxAlpha,plotRK,plotMRI,plotExtSTS)

  addpath('../RungeKutta')

  % shared plotting information
  box = [-3,0.5,-3,3];
  mname = 'SSP-DIRK2';
  fname = 'SSPDIRK2';

  % create base Butcher tables (including padding and stiff accuracy)
  zed = 0;  % zed = sym(0);
  one = 1;  % one = sym(1);
  two = 2;  % two = sym(2);
  five = 5;  % five = sym(5);
  seven = 7;  % seven = sym(7);
  twelve = 12;  % twelve = sym(12);
  gamma = one - one / sqrt(two);
  c = [zed; gama; gamma; one-gamma; one-gamma; one];
  Ai = [zed, zed, zed, zed, zed, zed;
        gamma, zed, zed, zed, zed, zed;
        zed, zed, gamma, zed, zed, zed;
        zed, zed, one-gamma, zed, zed, zed;
        zed, zed, one-two*gamma, zed, gamma, zed;
        zed, zed, one/two, zed, one/two, zed];
  bi = [zed, zed, one/two, zed, one/two, zed];
  di = [zed, zed, five/twelve, zed, seven/twelve, zed];
  Bi = [c, Ai; 2, bi; 1, di];
  Ae = 0;
  be = 0;
  de = 0;
  Be = 0;

  % verify properties and generate plots of ImEx-ARK method
  if (plotRK)
    disp('checking DIRK method properties')
    check_rk(c,Ai,bi,di,1,true,box,mname,fname);
  end

  % generate joint stability plot for this as an ImEx-MRI-GARK method
  if (plotMRI)

    % convert Butcher table to MRI "Omega" matrix
    [cmri, Gmri] = mis_to_mri(Bi);
    if (cmri ~= c)
      error('cmri does not match c for implicit MIS table')
    end

    % pad G with an initial row of zeros to match expected table structure,
    % and remove embedding row
    Gm = Gmri{1};
    G{1} = [zeros(1,6); Gm(1:end-1,:)];
    W = {};

    % set "dc" increment array (pad with initial 0)
    dc = [0; c(2:end)-c(1:end-1)];

    % test parameters
    thetavals = [10,30,45,60,80,90];
    numRay = 25;
    numGrid = 200;
    numAngle = 2;
    header = {};
    plottype = 'implicit_mrigark';
    plotcolors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4940 0.1840 0.5560], ...
                  [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    plotlinestyle = {'-','--','-.',':','-','--'};

    filename = ['mri_',fname,'_alpha_',num2str(maxAlpha),'.mat'];
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

    plottype = 'implicit';
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

    plotname = ['mri_',fname,'_alpha_',num2str(maxAlpha)];
    print('-dpng',plotname);
    savefig(plotname);

  end

  % generate joint stability plot for this as an ExtSTS method
  if (plotExtSTS)

    % test parameters
    box = [-1.5,0.5,-2,2];
    %thetavals = [10,30,45,60,80,90];
    thetavals = [0,30,45,60,90];
    %numRay = 25;
    numRay = 10;
    maxRho = 1e2;
    %numGrid = 200;
    numGrid = 60;
    numAngle = 2;
    header = {};
    plotcolors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4940 0.1840 0.5560], ...
                  [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    plotlinestyle = {'-','--','-.',':','-','--'};

    % RKC
    filename = ['extsts_',fname,'_rkc.mat'];
    q = matfile(filename,'Writable',true);
    q.box = box;
    q.thetavals = thetavals;
    q.numRay = numRay;
    q.maxRho = maxRho;
    q.numAngle = numAngle;
    q.numGrid = numGrid;
    figure
    hold on

    for itheta = 1:length(thetavals)
      maxTheta = thetavals(itheta);
      [xgrid,ygrid,Rmax] = extsts_jointstab(Ai, Ae, 'RKC', maxTheta, maxRho, numRay, numAngle, box, numGrid);
      R{itheta} = Rmax;
      contour(xgrid, ygrid, Rmax', [1+eps,1+eps], 'color', plotcolors{itheta}, 'LineStyle', ...
              plotlinestyle{itheta}, 'LineWidth', 2);
      header{itheta} = [num2str(maxTheta),char(176)];
      hold on
    end

    q.R = R;
    q.xgrid = xgrid;
    q.ygrid = ygrid;

    xl = box(1:2);  yl = box(3:4);
    xax = plot( linspace(xl(1),xl(2),10), zeros(1,10), 'k:');
    yax = plot( zeros(1,10), linspace(yl(1),yl(2),10), 'k:');
    hold off
    tstring = [mname,' ExtSTS method with RKC'];
    title(tstring);
    lgd = legend(header);
    lgd.Location = 'best';
    lgd.Title.String = '\theta values';

    plotfile = ['extsts_',fname,'_rkc'];
    print('-dpng',plotfile);
    savefig(plotfile);


    % RKL
    filename = ['extsts_',fname,'_rkl.mat'];
    q = matfile(filename,'Writable',true);
    q.box = box;
    q.thetavals = thetavals;
    q.numRay = numRay;
    q.maxRho = maxRho;
    q.numAngle = numAngle;
    q.numGrid = numGrid;
    figure
    hold on

    for itheta = 1:length(thetavals)
      maxTheta = thetavals(itheta);
      [xgrid,ygrid,Rmax] = extsts_jointstab(Ai, Ae, 'RKL', maxTheta, maxRho, numRay, numAngle, box, numGrid);
      R{itheta} = Rmax;
      contour(xgrid, ygrid, Rmax', [1+eps,1+eps], 'color', plotcolors{itheta}, 'LineStyle', ...
              plotlinestyle{itheta}, 'LineWidth', 2);
      header{itheta} = [num2str(maxTheta),char(176)];
      hold on
    end

    q.R = R;
    q.xgrid = xgrid;
    q.ygrid = ygrid;

    xl = box(1:2);  yl = box(3:4);
    xax = plot( linspace(xl(1),xl(2),10), zeros(1,10), 'k:');
    yax = plot( zeros(1,10), linspace(yl(1),yl(2),10), 'k:');
    hold off
    tstring = [mname,' ExtSTS method with RKL'];
    title(tstring);
    lgd = legend(header);
    lgd.Location = 'best';
    lgd.Title.String = '\theta values';

    plotfile = ['extsts_',fname,'_rkl'];
    print('-dpng',plotfile);
    savefig(plotfile);

  end
