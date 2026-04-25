function newdriver_heuneuler(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS)

  addpath('../RungeKutta')

  box = [-2.5,0.5,-2.5,2.5];
  mname = 'Heun-Euler';
  fname = 'HeunEuler';

  zed = 0;
  one = 1;
  two = 2;
  half = one/two;
  c = [zed; one; one];
  be = [half, half, zed];
  de = [one, zed, zed];
  Ae = [zed, zed, zed;
        one, zed, zed;
        half, half, zed];

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
                   'theta', [0], 3, 1e6, 1, 1, 0, 1, 60, ...
                   'explicit', 'best', '', '', ...
                   {'ExtSTS joint stability -- ','ExtSTS joint stability -- '});
  end

end
