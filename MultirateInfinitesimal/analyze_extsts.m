function analyze_extsts(embedding, mname, fname, Ae, be, Ai, bi, box, ...
                        sweepMode, sweepVals, numDiff, maxDiff, ...
                        numRxRadii, numRxAngle, maxRxAngle, maxRxRadius, ...
                        numGrid, baseMode, legendLocation, legendTitle, ...
                        fileSuffix, titlePrefixes)

  if embedding
    fprintf('\nPlotting ExtSTS joint stability region for %s embedding\n', mname)
  else
    fprintf('\nPlotting ExtSTS joint stability region for %s method\n', mname)
  end

  plotcolors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4940 0.1840 0.5560], ...
                [0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
  plotlinestyle = {'-','--','-.',':','-','--'};
  integrators = {'RKC', 'RKL'};

  for iint = 1:length(integrators)
    integrator = integrators{iint};
    filebase = extsts_filebase(fname, embedding, lower(integrator), fileSuffix);
    q = matfile([filebase, '.mat'], 'Writable', true);
    q.box = box;
    q.numDiff = numDiff;
    q.numGrid = numGrid;

    if strcmp(sweepMode, 'theta')
      q.thetavals = sweepVals;
      q.maxDiff = maxDiff;
      q.numRxRadii = numRxRadii;
      q.numRxAngle = numRxAngle;
      q.maxRxRadius = maxRxRadius;
    elseif strcmp(sweepMode, 'rho')
      q.maxDiff = sweepVals;
      q.numRxRadii = numRxRadii;
      q.numRxAngle = numRxAngle;
      q.maxRxAngle = maxRxAngle;
      q.maxRxRadius = maxRxRadius;
    else
      error('Unknown ExtSTS sweep mode "%s"', sweepMode)
    end

    fig = figure;
    stab_region(double(extsts_base_matrix(baseMode, Ae, Ai)), ...
                double(extsts_base_weights(baseMode, be, bi)), ...
                box, fig, 'k--', 'base');
    hold on

    [header, useHeader] = extsts_initial_legend(sweepMode, sweepVals, legendTitle);
    R = cell(1, length(sweepVals));

    for isweep = 1:length(sweepVals)
      if strcmp(sweepMode, 'theta')
        RxParams = [sweepVals(isweep), numRxAngle, maxRxRadius, numRxRadii];
        DiffParams = [maxDiff, numDiff];
      else
        RxParams = [maxRxAngle, numRxAngle, maxRxRadius, numRxRadii];
        DiffParams = [sweepVals(isweep), numDiff];
      end

      [xgrid, ygrid, Rmax] = extsts_jointstab(Ai, Ae, integrator, RxParams, DiffParams, box, numGrid);
      R{isweep} = Rmax;
      contour(xgrid, ygrid, Rmax', [1+eps,1+eps], 'color', plotcolors{isweep}, ...
              'LineStyle', plotlinestyle{isweep}, 'LineWidth', 2);

      if useHeader
        if strcmp(sweepMode, 'theta')
          header{isweep+1} = [num2str(sweepVals(isweep)), char(176)];
        else
          header{isweep+1} = ['\rho = ', num2str(sweepVals(isweep))];
        end
      end
    end

    q.R = R;
    q.xgrid = xgrid;
    q.ygrid = ygrid;

    xl = box(1:2);
    yl = box(3:4);
    plot(linspace(xl(1), xl(2), 10), zeros(1,10), 'k:');
    plot(zeros(1,10), linspace(yl(1), yl(2), 10), 'k:');
    hold off

    title(extsts_title(titlePrefixes{iint}, mname, embedding, integrator));
    if useHeader
      lgd = legend(header);
      lgd.Location = legendLocation;
      if ~isempty(legendTitle)
        lgd.Title.String = legendTitle;
      end
    else
      lgd = legend('Base', 'ExtSTS');
      lgd.Location = legendLocation;
    end

    print('-dpng', filebase);
    savefig(filebase);
  end

end

function filebase = extsts_filebase(fname, embedding, integrator, fileSuffix)

  filebase = ['extsts_', fname, '_'];
  if embedding
    filebase = [filebase, 'embedding_'];
  end
  filebase = [filebase, integrator];
  if ~isempty(fileSuffix)
    filebase = [filebase, '-', fileSuffix];
  end

end

function A = extsts_base_matrix(baseMode, Ae, Ai)

  if strcmp(baseMode, 'explicit')
    A = Ae;
  elseif strcmp(baseMode, 'implicit')
    A = Ai;
  else
    error('Unknown baseMode "%s"', baseMode)
  end

end

function b = extsts_base_weights(baseMode, be, bi)

  if strcmp(baseMode, 'explicit')
    b = be;
  elseif strcmp(baseMode, 'implicit')
    b = bi;
  else
    error('Unknown baseMode "%s"', baseMode)
  end

end

function [header, useHeader] = extsts_initial_legend(sweepMode, sweepVals, legendTitle)

  useHeader = strcmp(sweepMode, 'rho') || length(sweepVals) > 1 || ~isempty(legendTitle);
  if useHeader
    header = {'Base'};
  else
    header = {};
  end

end

function tstring = extsts_title(titlePrefix, mname, embedding, integrator)

  if embedding
    tstring = [titlePrefix, mname, ' embedding + ', integrator];
  else
    tstring = [titlePrefix, mname, ' + ', integrator];
  end

end
