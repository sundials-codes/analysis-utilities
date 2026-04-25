function newdriver_erk22a(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS)

  addpath('../RungeKutta')

  box = [-5,0.5,-3,3];
  mname = 'ERK22a';
  fname = 'ERK22a';

  c = [0; 1];
  be = [0, 1];
  de = [1, 0];
  Ae = [0, 0;
        0.5, 0];

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
