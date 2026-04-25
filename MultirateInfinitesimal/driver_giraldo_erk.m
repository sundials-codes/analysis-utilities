function newdriver_giraldo_erk(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS)

  addpath('../RungeKutta')

  box = [-3.5,0.5,-3,3];
  mname = 'Giraldo-ERK2';
  fname = 'GiraldoERK2';

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
