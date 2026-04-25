function newdriver_giraldo(maxAlpha,embedding,plotImEx,plotImExMRI,plotExtSTS)

  addpath('../RungeKutta')

  box = [-3,0.5,-3,3];
  mname = 'Giraldo-ARK2';
  fname = 'GiraldoARK2';

  zed = 0;
  one = 1;
  two = 2;
  three = 3;
  four = 4;
  six = 6;
  eight = 8;
  c = [zed; two-sqrt(two); two-sqrt(two); one; one; one];
  be = [one/(two*sqrt(two)), zed, one/(two*sqrt(two)), zed, (sqrt(two)-one)/sqrt(two), zed];
  de = [(four-sqrt(two))/eight, zed, (four-sqrt(two))/eight, zed, one/(two*sqrt(two)), zed];
  Ae = [zed, zed, zed, zed, zed, zed;
        two-sqrt(two), zed, zed, zed, zed, zed;
        two-sqrt(two), zed, zed, zed, zed, zed;
        (three-two*sqrt(two))/six, zed, (three+two*sqrt(two))/six, zed, zed, zed;
        (three-two*sqrt(two))/six, zed, (three+two*sqrt(two))/six, zed, zed, zed;
        one/(two*sqrt(two)), zed, one/(two*sqrt(two)), zed, (sqrt(two)-one)/sqrt(two), zed];
  bi = be;
  di = de;
  Ai = [zed, zed, zed, zed, zed, zed;
        two-sqrt(two), zed, zed, zed, zed, zed;
        (sqrt(two)-one)/sqrt(two), zed, (sqrt(two)-one)/sqrt(two), zed, zed, zed;
        zed, zed, one, zed, zed, zed;
        one/(two*sqrt(two)), zed, one/(two*sqrt(two)), zed, (sqrt(two)-one)/sqrt(two), zed;
        one/(two*sqrt(two)), zed, one/(two*sqrt(two)), zed, (sqrt(two)-one)/sqrt(two), zed];

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
