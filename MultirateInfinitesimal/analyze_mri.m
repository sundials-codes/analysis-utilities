function analyze_mri(maxAlpha, embedding, mname, fname, box, c, Be, Bi, ...
                     thetavals, jointPlotType, basePlotType, filePrefix)

  if embedding
    fprintf('\nPlotting MRI joint stability region for %s embedding\n', mname)
  else
    fprintf('\nPlotting MRI joint stability region for %s method\n', mname)
  end

  nstages = length(c);
  [c, G, W] = mri_coupling_tables(Be, Bi, nstages);
  dc = [0; c(2:end)-c(1:end-1)];

  numRay = 25;
  numGrid = 200;
  numAngle = 2;
  header = {};
  plotcolors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4940 0.1840 0.5560], ...
                [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
  plotlinestyle = {'-','--','-.',':','-','--'};

  filebase = [filePrefix, '_', fname, '_alpha_', num2str(maxAlpha)];
  if embedding
    filebase = [filebase, '_embedding'];
  end

  q = matfile([filebase, '.mat'], 'Writable', true);
  q.box = box;
  q.maxAlpha = maxAlpha;
  q.thetavals = thetavals;
  q.numRay = numRay;
  q.numAngle = numAngle;
  q.numGrid = numGrid;

  R = cell(1, length(thetavals));
  figure
  hold on

  for itheta = 1:length(thetavals)
    maxTheta = thetavals(itheta);
    [xgrid, ygrid, Rmax] = imexmri_jointstab(G, W, dc, maxAlpha, maxTheta, ...
                                             numRay, numAngle, box, numGrid, jointPlotType);
    R{itheta} = Rmax;
    contour(xgrid, ygrid, Rmax', [1+eps,1+eps], 'color', plotcolors{itheta}, ...
            'LineStyle', plotlinestyle{itheta}, 'LineWidth', 2);
    header{itheta} = [num2str(maxTheta), char(176)];
  end

  q.R = R;
  q.xgrid = xgrid;
  q.ygrid = ygrid;

  maxTheta = 90;
  numRay = 50;
  numGrid = 50;
  [xgrid, ygrid, Rmax] = imexmri_jointstab(G, W, dc, maxAlpha, maxTheta, ...
                                           numRay, numAngle, box, numGrid, basePlotType);
  contour(xgrid, ygrid, Rmax', [1+eps,1+eps], 'k:', 'LineWidth', 2);
  header{end+1} = ['Base (', num2str(maxTheta), char(176), ')'];
  q.Rbase = Rmax;
  q.xbase = xgrid;
  q.ybase = ygrid;

  hold off
  title(['\alpha = ', num2str(maxAlpha), char(176)]);
  lgd = legend(header);
  lgd.Location = 'best';
  lgd.Title.String = '\theta values';

  print('-dpng', filebase);
  savefig(filebase);

end

function [cmri, G, W] = mri_coupling_tables(Be, Bi, nstages)

  if ~isscalar(Be) || Be ~= 0
    [cmri, Wmri] = mis_to_mri(Be);
    c = cmri;
    nstages = length(c);
    Wm = Wmri{1};
    W{1} = [zeros(1,nstages); Wm(1:end-1,:)];
  else
    W{1} = [];
    c = [];
  end

  if ~isscalar(Bi) || Bi ~= 0
    [cmri, Gmri] = mis_to_mri(Bi);
    nstages = length(cmri);
    if (~isempty(c)) && any(cmri ~= c)
      error('cmri does not match c for implicit MIS table')
    end
    Gm = Gmri{1};
    G{1} = [zeros(1,nstages); Gm(1:end-1,:)];
  else
    G{1} = [];
  end

  if isempty(W{1}) && ~isempty(G{1})
    W{1} = 0*G{1};
  elseif isempty(G{1}) && ~isempty(W{1})
    G{1} = 0*W{1};
  end

end
