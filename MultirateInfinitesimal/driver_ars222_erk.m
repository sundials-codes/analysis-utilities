function newdriver_ars222_erk(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS)

  addpath('../RungeKutta')

  box = [-3,0.5,-3,3];
  mname = 'Ascher(2,2,2)-ERK';
  fname = 'ARS222-ERK';

  zed = 0;
  one = 1;
  two = 2;
  gamma = (two-sqrt(two))/two;
  delta = one-one/(two*gamma);
  c = [zed; gamma; gamma; one; one];
  be = [delta, zed, one-delta, zed, zed];
  de = [zed, zed, 0.6, zed, 0.4];
  Ae = [zed, zed, zed, zed, zed;
        gamma, zed, zed, zed, zed;
        gamma, zed, zed, zed, zed;
        delta, zed, one-delta, zed, zed;
        delta, zed, one-delta, zed, zed];

  if embedding
    btmp = be;
    be = de;
    de = btmp;
  end

  Be = [c, Ae; 2, be; 1, de];
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
                   'theta', [0], 3, 1e2, 3, 2, 0, 1, 60, ...
                   'explicit', 'best', '', '', ...
                   {'ExtSTS joint stability -- ','ExtSTS joint stability -- '});
  end

end
