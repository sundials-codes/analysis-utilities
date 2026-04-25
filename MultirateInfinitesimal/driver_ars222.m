function newdriver_ars222(maxAlpha,embedding,plotImEx,plotImExMRI,plotExtSTS)

  addpath('../RungeKutta')

  box = [-3,0.5,-3,3];
  mname = 'Ascher(2,2,2)';
  fname = 'ARS222';

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
  bi = [zed, zed, one-gamma, zed, gamma];
  di = de;
  Ai = [zed, zed, zed, zed, zed;
        gamma, zed, zed, zed, zed;
        zed, zed, gamma, zed, zed;
        zed, zed, one, zed, zed;
        zed, zed, one-gamma, zed, gamma];

  if embedding
    btmp = be;
    be = de;
    de = btmp;
  end
  Be = [c, Ae; 2, be; 1, de];

  if embedding
    btmp = bi;
    bi = di;
    di = btmp;
  end
  Bi = [c, Ai; 2, bi; 1, di];

  if plotImEx
    fprintf('\nChecking ImEx-ARK method properties for %s method\n', mname)
    check_ark_embedded(c,c,Ae,Ai,be,bi,de,di,1e-11,1,true,box,mname,fname);
  end

  if plotImExMRI
    analyze_mri(maxAlpha, embedding, mname, fname, box, c, Be, Bi, ...
                [20,40,60,80], 'imexmri', 'imex', 'imexmri');
  end

  if plotExtSTS
    analyze_extsts(embedding, mname, fname, Ae, be, Ai, bi, box, ...
                   'theta', [20,40,60,80], 3, 1e6, 3, 2, 0, 1e2, 60, ...
                   'explicit', 'best', '\theta values', '', ...
                   {'ExtSTS joint stability -- ','ExtSTS joint stability -- '});
  end

end
