%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2018, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------
% driver to check analytical properties of various RK methods
% in butcher.m

clear; close all
!\rm new_tests.txt
diary new_tests.txt
fontsize(14, "points")

% test parameters
plotdir = '.';        % folder to use for stability region plots
plotregions = true;   % create stability region plots
use_symbolic = false; % use symbolic storage for tables
reportL = 1;          % level of reporting output from check_* routines
do_explicit = true;
do_implicit = true;
do_imex = true;

% explicit methods
if (do_explicit)
  %         {bname,                        mname,                      filename,              stability region box}
  tests = {
           {'SSP(2,2)-ERK',          'SSP-ERK-2-1-2',        'ARKODE_SSP_ERK_2_1_2',       [-4,0.5,-3,3]},
           {'SSP(3,2)-ERK',          'SSP-ERK-3-1-2',        'ARKODE_SSP_ERK_3_1_2',       [-6,0.5,-4,4]},
           {'SSP(4,2)-ERK',          'SSP-ERK-4-1-2',        'ARKODE_SSP_ERK_4_1_2',       [-8,0.5,-5,5]},
           {'SSP(10,2)-ERK',         'SSP-ERK-10-1-2',       'ARKODE_SSP_ERK_10_1_2',      [-20,2,-15,15]},
           {'SSP(4,3)-ERK',          'SSP-ERK-4-2-3',        'ARKODE_SSP_ERK_4_2_3',       [-8,1,-5,5]},
           {'SSP(9,3)-ERK',          'SSP-ERK-9-2-3',        'ARKODE_SSP_ERK_9_2_3',       [-20,2,-10,10]},
           {'SSP(10,4)-ERK',         'SSP-ERK-10-3-4',       'ARKODE_SSP_ERK_10_3_4',      [-20,2,-10,10]},
           {'SSP2(3,3,2)-lspum-ERK', 'SSP-LSPUM-ERK-3-1-2', 'ARKODE_SSP_LSPUM_ERK_3_1_2', [-3.5,0.5,-3.5,3.5]},
           {'Giraldo-ARK2-ERK',      'ARK2-ERK-3-1-2',       'ARKODE_ARK2_ERK_3_1_2',      [-3.5,0.5,-3.5,3.5]},
           {'Ascher(2,2,2)-ERK',     'ASCHER-ERK-3-1-2',     'ARKODE_ASCHER_ERK_3_1_2',    [-3.5,0.5,-3.5,3.5]},
     };
  fprintf('                             |    |  Method |  Embedding  |\n');
  fprintf('            Name             |  s |  q  lq  |  p  lp      | tol\n');
  fprintf('  -----------------------------------------------------------------\n');
  for i = 1:length(tests)
    bname = tests{i}{1};
    mname = tests{i}{2};
    fname = [ plotdir, '/', tests{i}{3} ];
    box = tests{i}{4};
    B = butcher(bname,use_symbolic);
    s = size(B,2)-1;
    [q,p,qs,lq,lp,tol] = check_rk(B,reportL,plotregions,box,mname,fname);
    fprintf(' %26s  | %2i | %2i  %2i  | %2i  %2i      | %.0e\n', mname,s,q,lq,p,lp,tol);
  end
  fprintf('  -----------------------------------------------------------------\n');
  fprintf('\n');
end

% implicit methods
if (do_implicit)
  %         {bname,                          name,                        filename,                  stability region box}
  tests = {
           {'SSP(2,2)-SDIRK',          'SSP-SDIRK-2-1-2',        'ARKODE_SSP_SDIRK_2_1_2',       [-10,20,-15,15]},
           {'SSP(3,2)-DIRK',           'SSP-DIRK-3-1-2',         'ARKODE_SSP_DIRK_3_1_2',        [-10,20,-15,15]},
           {'SSP2(3,3,2)-lspum-SDIRK', 'SSP-LSPUM-SDIRK-3-1-2', 'ARKODE_SSP_LSPUM_SDIRK_3_1_2', [-10,20,-15,15]},
           {'SSP(4,3)-ESDIRK',         'SSP-ESDIRK-4-2-3',       'ARKODE_SSP_ESDIRK_4_2_3',      [-10,20,-15,15]},
           {'Giraldo-ARK2-ESDIRK',     'ARK2-EDIRK-3-1-2',       'ARKODE_ARK2_EDIRK_3_1_2',      [-10,20,-15,15]},
           {'Ascher(2,2,2)-SDIRK',     'ASCHER-SDIRK-3-1-2',     'ARKODE_ASCHER_SDIRK_3_1_2',    [-10,20,-15,15]},
     };

  fprintf('                             |    |       Method        |      Embedding      |\n');
  fprintf('            Name             |  s |  q  lq   A   B   L  |  p  lp   A   B   L  |  qs  tol\n');
  fprintf('  ----------------------------------------------------------------------------------------\n');
  for i = 1:length(tests)
    bname = tests{i}{1};
    mname = tests{i}{2};
    fname = [ plotdir, '/', tests{i}{3} ];
    box = tests{i}{4};
    B = butcher(bname,use_symbolic);
    s = size(B,2)-1;
    [q,p,qs,lq,lp,tol,Bs,As,Ls,BsE,AsE,LsE] = check_rk(B,reportL,plotregions,box,mname,fname);
    Bs_ = ' ';
    As_ = ' ';
    Ls_ = ' ';
    BsE_ = ' ';
    AsE_ = ' ';
    LsE_ = ' ';
    if (Bs==1)
      Bs_ = 'Y';
    end
    if (As==1)
      As_ = 'Y';
    end
    if (Ls==1)
      Ls_ = 'Y';
    end
    if (BsE==1)
      BsE_ = 'Y';
    end
    if (AsE==1)
      AsE_ = 'Y';
    end
    if (LsE==1)
      LsE_ = 'Y';
    end
    fprintf(' %26s  | %2i | %2i  %2i   %s   %s   %s  | %2i  %2i   %s   %s   %s  | %2i  %.0e\n',...
            mname,s,q,lq,As_,Bs_,Ls_,p,lp,AsE_,BsE_,LsE_,qs,tol);
  end
  fprintf('  ----------------------------------------------------------------------------------------\n');
  fprintf('\n');
end


% imex methods
if (do_imex)
  %         {expname,                   impname,                 method name              filename,              stab region box}
  tests = {
           {'SSP(2,2)-ERK',          'SSP(2,2)-SDIRK',          'SSP-ARK-2-1-2',        'ARKODE_SSP_ARK_2_1_2',       [-2.5,0.5,-2.5,2.5]},
           {'Ascher(2,2,2)-ERK',     'Ascher(2,2,2)-SDIRK',     'ASCHER-ARK-3-1-2',     'ARKODE_ASCHER_ARK_3_1_2',    [-2.5,0.5,-2.5,2.5]},
           {'SSP(3,2)-ERK',          'SSP(3,2)-DIRK',           'SSP-ARK-3-1-2',        'ARKODE_SSP_ARK_3_1_2',       [-5,0.5,-3.5,3.5]},
           {'SSP2(3,3,2)-lspum-ERK', 'SSP2(3,3,2)-lspum-SDIRK', 'SSP-LSPUM-ARK-3-1-2',  'ARKODE_SSP_LSPUM_ARK_3_1_2', [-3.5,0.5,-3.5,3.5]},
           {'Giraldo-ARK2-ERK',      'Giraldo-ARK2-ESDIRK',     'ARK2-3-1-2',           'ARKODE_ARK2_3_1_2',          [-2,0.5,-2,2]},
           {'SSP(4,3)-ERK',          'SSP(4,3)-ESDIRK',         'SSP-ARK-4-2-3',        'ARKODE_SSP_ARK_4_2_3',       [-3.5,0.5,-2.5,2.5]},
           {'ARK3(2)4L[2]SA-ERK',    'ARK3(2)4L[2]SA-ESDIRK',   'ARK324L2SA-4-2-3',     'ARKODE_ARK324L2SA_4_2_3',    [-4,0.5,-4,4]},
           {'ARK4(3)6L[2]SA-ERK',    'ARK4(3)6L[2]SA-ESDIRK',   'ARK436L2SA-6-3-4',     'ARKODE_ARK436L2SA_6_3_4',    [-5,1,-5,5]},
           {'ARK4(3)7L[2]SA-ERK',    'ARK4(3)7L[2]SA-ESDIRK',   'ARK437L2SA-7-3-4',     'ARKODE_ARK437L2SA_7_3_4',    [-4,0.5,-3,3]},
           {'ARK5(4)8L[2]SA-ERK',    'ARK5(4)8L[2]SA-ESDIRK',   'ARK548L2SA-8-4-5',     'ARKODE_ARK548L2SA_8_4_5',    [-5,1,-4,4]},
           {'ARK5(4)8L[2]SAb-ERK',   'ARK5(4)8L[2]SAb-ESDIRK',  'ARK548L2SAb-8-4-5',    'ARKODE_ARK548L2SAb_8_4_5',   [-5,1,-4,4]},
     };

  fprintf('                             |     | Method | Embedding |\n');
  fprintf('            Name             |  s  |   q    |    p      | qs\n');
  fprintf('  ---------------------------------------------------------------\n');
  for i = 1:length(tests)
    ename = tests{i}{1};
    iname = tests{i}{2};
    mname = tests{i}{3};
    fname = [ plotdir, '/', tests{i}{4} ];
    box = tests{i}{5};
    Be = butcher(ename,use_symbolic);
    Bi = butcher(iname,use_symbolic);
    s = size(Be,2)-1;
    [Ae,be,ce,de,~,~] = extract_butcher_embedded(Be);
    [Ai,bi,ci,di,~,~] = extract_butcher_embedded(Bi);
    [qE,qI,q,pE,pI,p,qsE,qsI,qsA] = check_ark_embedded(ce,ci,Ae,Ai,be,bi,de,di,1e-11,reportL,plotregions,box,mname,fname);
    fprintf(' %26s  | %2i  |  %2i    |   %2i      | %2i\n', mname, s, q, p, qsA);
  end
  fprintf('  -----------------------------------------------------------------\n');
  fprintf('\n');
end

diary off
% end of script
