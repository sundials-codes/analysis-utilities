function newdriver_irk21a(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS)

  addpath('../RungeKutta')

  box = [-5,25,-15,15];
  mname = 'IRK21a';
  fname = 'IRK21a';

  c = [0; 1; 1];
  Ai = [0, 0, 0;
        1, 0, 0;
        0.5, 0, 0.5];
  bi = [0.5, 0, 0.5];
  di = [0, 0, 1];

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
    analyze_extsts(embedding, mname, fname, Ae, be, Ai, bi, box, ...
                   'rho', [1, 1e2, 1e4, 1e6], 3, [1, 1e2, 1e4, 1e6], ...
                   1, 1, 0, 1, 60, 'implicit', 'northeast', '\rho values', '', ...
                   {'ExtSTS joint stability -- ','ExtSTS joint stability -- '});
  end

end
