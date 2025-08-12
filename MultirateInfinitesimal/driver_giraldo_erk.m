function driver_giraldo_erk(maxAlpha,plotRK,plotMRI,plotExtSTS)

  addpath('../RungeKutta')

  % shared plotting information
  box = [-3.5,0.5,-3,3];
  mname = 'Giraldo-ERK2';
  fname = 'GiraldoERK2';

  % create base Butcher tables (including padding and stiff accuracy)
  zed = 0;  % zed = sym(0);
  one = 1;  % one = sym(1);
  two = 2;  % two = sym(2);
  three = 3;  % three = sym(3);
  four = 4;  % four = sym(4);
  six = 6;  % six = sym(6);
  eight = 8;  % eight = sym(8);
  c = [zed; two-sqrt(two); two-sqrt(two); one; one; one];
  be = [one/(two*sqrt(two)), zed, one/(two*sqrt(two)), zed, (sqrt(two)-one)/sqrt(two), zed];
  de = [(four-sqrt(two))/eight, zed, (four-sqrt(two))/eight, zed, one/(two*sqrt(two)), zed];
  Ae = [zed, zed, zed, zed, zed, zed;
        two-sqrt(two), zed, zed, zed, zed, zed;
        two-sqrt(two), zed, zed, zed, zed, zed;
        (three-two*sqrt(two))/six, zed, (three+two*sqrt(two))/six, zed, zed, zed;
        (three-two*sqrt(two))/six, zed, (three+two*sqrt(two))/six, zed, zed, zed;
        one/(two*sqrt(two)), zed, one/(two*sqrt(two)), zed, (sqrt(two)-one)/sqrt(two), zed];
  Be = [c, Ae; 2, be; 1, de];
  bi = 0;
  di = 0;
  Ai = 0;
  Bi = 0;

  % verify properties and generate plots of RK method
  if (plotRK)
    fprintf('\nChecking ERK method properties for %s method\n', mname)
    check_rk(Be,1,true,box,mname,fname);
  end

  % generate joint stability plot for this as an ExtSTS method
  if (plotExtSTS)
    fprintf('\nPlotting ExtSTS joint stability region for %s method\n', mname)

    % test parameters
    box = [-3.5,0.5,-3,3];
    thetavals = [0];  % maxRxAngle values
    numDiff = 3;
    maxDiff = 1e2;
    numRxRadii = 3;
    numRxAngle = 2;
    maxRxRadius = 1;
    numGrid = 60;
    plotcolors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4940 0.1840 0.5560], ...
                  [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    plotlinestyle = {'-','--','-.',':','-','--'};

    % RKC
    filename = ['extsts_',fname,'_rkc.mat'];
    q = matfile(filename,'Writable',true);
    q.box = box;
    q.thetavals = thetavals;
    q.numDiff = numDiff;
    q.maxDiff = maxDiff;
    q.numRxRadii = numRxRadii;
    q.numRxAngle = numRxAngle;
    q.maxRxRadius = maxRxRadius;
    q.numGrid = numGrid;
    fig = figure;
    stab_region(double(Ae),double(be),box,fig,'k--','base');  % base explicit method stability region
    hold on

    for itheta = 1:length(thetavals)
      maxTheta = thetavals(itheta);
      RxParams = [maxTheta, numRxAngle, maxRxRadius, numRxRadii];
      DiffParams = [maxDiff, numDiff];
      [xgrid,ygrid,Rmax] = extsts_jointstab(Ai, Ae, 'RKC', RxParams, DiffParams, box, numGrid);
      R{itheta} = Rmax;
      contour(xgrid, ygrid, Rmax', [1+eps,1+eps], 'color', plotcolors{itheta}, 'LineStyle', ...
              plotlinestyle{itheta}, 'LineWidth', 2);
    end

    q.R = R;
    q.xgrid = xgrid;
    q.ygrid = ygrid;

    xl = box(1:2);  yl = box(3:4);
    xax = plot( linspace(xl(1),xl(2),10), zeros(1,10), 'k:');
    yax = plot( zeros(1,10), linspace(yl(1),yl(2),10), 'k:');
    hold off
    tstring = ['ExtSTS joint stability -- ', mname,' + RKC'];
    title(tstring);
    lgd = legend('Base','ExtSTS');
    lgd.Location = 'best';

    plotfile = ['extsts_',fname,'_rkc'];
    print('-dpng',plotfile);
    savefig(plotfile);


    % RKL
    filename = ['extsts_',fname,'_rkl.mat'];
    q = matfile(filename,'Writable',true);
    q.box = box;
    q.thetavals = thetavals;
    q.numDiff = numDiff;
    q.maxDiff = maxDiff;
    q.numRxRadii = numRxRadii;
    q.numRxAngle = numRxAngle;
    q.maxRxRadius = maxRxRadius;
    q.numGrid = numGrid;
    fig = figure;
    stab_region(double(Ae),double(be),box,fig,'k--','base');  % base method stability region
    hold on

    for itheta = 1:length(thetavals)
      maxTheta = thetavals(itheta);
      RxParams = [maxTheta, numRxAngle, maxRxRadius, numRxRadii];
      DiffParams = [maxDiff, numDiff];
      [xgrid,ygrid,Rmax] = extsts_jointstab(Ai, Ae, 'RKL', RxParams, DiffParams, box, numGrid);
      R{itheta} = Rmax;
      contour(xgrid, ygrid, Rmax', [1+eps,1+eps], 'color', plotcolors{itheta}, 'LineStyle', ...
              plotlinestyle{itheta}, 'LineWidth', 2);
    end

    q.R = R;
    q.xgrid = xgrid;
    q.ygrid = ygrid;

    xl = box(1:2);  yl = box(3:4);
    xax = plot( linspace(xl(1),xl(2),10), zeros(1,10), 'k:');
    yax = plot( zeros(1,10), linspace(yl(1),yl(2),10), 'k:');
    hold off
    tstring = ['ExtSTS joint stability -- ', mname,' + RKL'];
    title(tstring);
    lgd = legend('Base','ExtSTS');
    lgd.Location = 'best';

    plotfile = ['extsts_',fname,'_rkl'];
    print('-dpng',plotfile);
    savefig(plotfile);

  end

  % generate joint stability plot for this as an ImEx-MRI-GARK method
  if (plotMRI)
    fprintf('\nPlotting MRI joint stability region for %s method\n', mname)

    % convert Butcher tables to MRI "Gamma" and "Omega" matrices
    [cmri, Wmri] = mis_to_mri(Be);
    if (cmri ~= c)
      error('cmri does not match c for explicit MIS table')
    end

    % pad W with an initial row of zeros to match expected table structure,
    % and remove embedding row
    Wm = Wmri{1};
    W{1} = [zeros(1,6); Wm(1:end-1,:)];
    G{1} = 0*W{1};;

    % set "dc" increment array (pad with initial 0)
    dc = [0; c(2:end)-c(1:end-1)];

    % test parameters
    thetavals = [10,30,45,60,80,90];
    numRay = 25;
    numGrid = 200;
    numAngle = 2;
    header = {};
    plottype = 'explicit_mrigark';
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

    plottype = 'explicit';
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
