function driver_ssp43(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS)

  addpath('../RungeKutta')

  box = [-8,1,-5,5];
  mname = 'SSP(4,3)';
  fname = 'SSP43';

  Be = butcher('SSP(4,3)-ERK');
  s = size(Be,2)-1;
  c = Be(1:s,1);
  Ae = Be(1:s,2:s+1);
  be = Be(s+1,2:s+1);
  de = Be(s+2,2:s+1);

  if embedding
    btmp = be;
    be = de;
    de = btmp;
  end
  Be = [c, Ae; 3, be; 2, de];
  Ai = 0;
  bi = 0;
  di = 0;

  if plotRK
    fprintf('\nChecking ERK method properties for %s method\n', mname)
    check_rk(Be,1,true,box,mname,fname);
  end

  if plotMRI
    analyze_mri(maxAlpha, embedding, mname, fname, box, c, Be, 0, ...
                [10,30,45,60,80,90], 'explicit_mrigark', 'explicit', 'mri');
  end

  if plotExtSTS
    analyze_extsts(embedding, mname, fname, Ae, be, Ai, bi, box, ...
                   'theta', [0], 3, 1e6, 1, 1, 0, 1, 60, ...
                   'explicit', 'best', '', '', ...
                   {'ExtSTS joint stability -- ','ExtSTS joint stability -- '});
  end

end
