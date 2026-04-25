function newdriver_giraldo_dirk(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS)

  addpath('../RungeKutta')

  box = [-5,25,-15,15];
  mname = 'Giraldo-DIRK2';
  fname = 'GiraldoDIRK2';

  zed = 0;
  one = 1;
  two = 2;
  four = 4;
  eight = 8;
  c = [zed; two-sqrt(two); two-sqrt(two); one; one; one];
  Ai = [zed, zed, zed, zed, zed, zed;
        two-sqrt(two), zed, zed, zed, zed, zed;
        (sqrt(two)-one)/sqrt(two), zed, (sqrt(two)-one)/sqrt(two), zed, zed, zed;
        zed, zed, one, zed, zed, zed;
        one/(two*sqrt(two)), zed, one/(two*sqrt(two)), zed, (sqrt(two)-one)/sqrt(two), zed;
        one/(two*sqrt(two)), zed, one/(two*sqrt(two)), zed, (sqrt(two)-one)/sqrt(two), zed];
  bi = [one/(two*sqrt(two)), zed, one/(two*sqrt(two)), zed, (sqrt(two)-one)/sqrt(two), zed];
  di = [(four-sqrt(two))/eight, zed, (four-sqrt(two))/eight, zed, one/(two*sqrt(two)), zed];

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
                   'rho', [1, 1e2, 1e4, 1e6], 3, [1, 1e2, 1e4, 1e6], ...
                   1, 1, 0, 1, 60, 'implicit', 'northeast', '', '', ...
                   {'ExtSTS joint stability -- ','ExtSTS joint stability -- '});

    analyze_extsts(embedding, mname, fname, Ae, be, Ai, bi, [-50,5000,-10000,10000], ...
                   'rho', [1, 1e2, 1e4, 1e6], 3, [1, 1e2, 1e4, 1e6], ...
                   1, 1, 0, 1, 60, 'implicit', 'northeast', '', 'zoomout', ...
                   {'ExtSTS joint stability -- ','ExtSTS joint stability -- '});
  end

end
