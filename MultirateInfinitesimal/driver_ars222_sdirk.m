function newdriver_ars222_sdirk(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS)

  addpath('../RungeKutta')

  box = [-5,110,-55,55];
  mname = 'Ascher(2,2,2)-SDIRK';
  fname = 'ARS222-SDIRK';

  zed = 0;
  one = 1;
  two = 2;
  gamma = (two-sqrt(two))/two;
  c = [zed; gamma; gamma; one; one];
  di = [zed, zed, 0.6, zed, 0.4];
  bi = [zed, zed, one-gamma, zed, gamma];
  Ai = [zed, zed, zed, zed, zed;
        gamma, zed, zed, zed, zed;
        zed, zed, gamma, zed, zed;
        zed, zed, one, zed, zed;
        zed, zed, one-gamma, zed, gamma];

  if embedding
    btmp = bi;
    bi = di;
    di = btmp;
  end

  Bi = [c, Ai; 2, bi; 1, di];
  Ae = 0;
  be = 0;

  if plotRK
    fprintf('\nChecking DIRK method properties for %s method\n', mname)
    check_rk(Bi,1,true,box,mname,fname);
  end

  if plotMRI
    analyze_mri(maxAlpha, embedding, mname, fname, box, c, 0, Bi, ...
                [10,30,45,60,80,90], 'implicit_mrigark', 'implicit', 'mri');
  end

  if plotExtSTS
    analyze_extsts(embedding, mname, fname, Ae, be, Ai, bi, [-5,110,-55,55], ...
                   'theta', [0], 3, 1e2, 3, 2, 0, 1, 60, ...
                   'explicit', 'best', '', '', ...
                   {'ExtSTS joint stability -- ','ExtSTS joint stability -- '});
  end

end
