function [qE,qI,q,pE,pI,p,qsE,qsI,qsA] = check_ark_embedded(cE,cI,AE,AI,bE,bI,dE,dI,tol,reportL,doPlot,box,mname,fname)
% Usage: [qE,qI,q,pE,pI,p,qsE,qsI,qsA] = check_ark_embedded(cE,cI,AE,AI,bE,bI,dE,dI,tol,reportL,doPlot,box,mname,fname)
%
% Checks the pair of Butcher tables given by
%
%    cE | AE     cI | AI
%    --------    --------
%       | bE        | bI
%       | dE        | dI
%
% to determine the analytical order of accuracy for each method
% alone (qE and qI), their embeddings (pE and pI), as well as the
% combined ARK method order (q) and combined ARK embedding order
% (p).  This function uses the order condition equations from the
% article Sandu & Gunther, SINUM 53, 2015.
%
% We also compute the 'stage order' of the ERK, DIRK and ARK methods,
% respectively, through checking the Butcher simplifying assumption C(q)
% as prescribed in the article
% J.C. Butcher, "Implicit Runge-Kutta processes", Math Comp 18 (1964), pp 50-64.
%
% Inputs:
%   cE, AE, bE, dE -- ERK Butcher table (dE contains embedding)
%   cI, AI, bI, dI -- second Butcher table (dI contains embedding)
%   tol -- tolerance for checking order conditions (use negative
%          or zero for default value)
%   reportL -- integer flag denoting print level:
%              >1 -> all results
%              =1 -> final results
%            else -> no results
%   doPlot -- boolean flag denoting whether to create/save plot
%   box -- [xl,xr,yl,yr] sub-region of the complex plane to
%          use in plot (ignored if doPlot is false)
%   mname -- string containing method name to insert into plot title
%   fname -- string containing method name to insert into plot filenames
%
% Outputs:
%   qE -- order of accuracy for ERK method
%   qI -- order of accuracy for DIRK method
%   q -- order of accuracy for combined ARK method
%   pE -- order of accuracy for ERK embedding
%   pI -- order of accuracy for DIRK embedding
%   p -- order of accuracy for combined ARK embedding
%   qsE -- stage order of accuracy for ERK method
%   qsI -- stage order of accuracy for DIRK method
%   qsA -- stage order of accuracy for combined ARK method
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2016, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% set tolerance on 'equality'
if (tol <= 0)
   tol = 1e-8;
end

% check individual Butcher tables using single-table routine
BE = [cE, AE; 0, bE; 0, dE];
[qE,pE] = check_butcher(BE);
BI = [cI, AI; 0, bI; 0, dI];
[qI,pI,BsImp,AsImp,LsImp,BsEmb,AsEmb,LsEmb] = check_butcher(BI);

% check stage order of component and overall methods
[qsE, qsI, qsA] = check_ark_stage_order(cE,cI,AE,AI,tol);

% determine whether either method satisfies the 'stiffly accurate' property
expSA = true;
impSA = true;
for i=1:length(bE)
   if (abs(bE(i) - AE(end,i)) > tol)
      expSA = false;
      break;
   end
end
for i=1:length(bI)
   if (abs(bI(i) - AI(end,i)) > tol)
      impSA = false;
      break;
   end
end

% initialize failure flags, order of accuracy
failed = false;
failE = false;
failI = false;
failC = false;
q = -1;

% convert bE and bI to column vectors for these tests
bE = reshape(bE,length(bE),1);
bI = reshape(bI,length(bI),1);

% create cell arrays to streamline pairwise testing
bs = {bE,bI};
As = {AE,AI};
cs = {cE,cI};

% ensure row sum condition c(i) = sum(A(i,:)) is satisfied for each method
[s,tmp] = size(AE);
for i=1:s
   tst = cE(i) - sum(AE(i,:));
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails row sum condition (i = %i, tst = %g)\n', i, tst);
      end
   end
   tst = cI(i) - sum(AI(i,:));
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails row sum condition (i = %i, tst = %g)\n', i, tst);
      end
   end
end
if (~failed)
   q = 0;
end
if (~failE)
   qE = max(qE,0);
end
if (~failI)
   qI = max(qI,0);
end

% check for first order
if (~failed || ~failE || ~failI)
   tst = sum(bE)-1;
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 1st order condition (tst = %g)\n', tst);
      end
   end

   tst = sum(bI)-1;
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 1st order condition (tst = %g)\n', tst);
      end
   end
      if (reportL>1)
         if (~failE), fprintf('  ERK method passes order 1 conditions\n'); end
         if (~failI), fprintf('  DIRK method passes order 1 conditions\n'); end
      end
      if (~failed)
         q = 1;
      end
      if (~failE)
         qE = max(qE,1);
      end
      if (~failI)
         qI = max(qI,1);
      end
end

% check for second order
if (~failed || ~failE || ~failI)

   tst = bE'*cE - 0.5;
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 2nd order condition (tst = %g)\n', tst);
      end
   end

   tst = bI'*cI - 0.5;
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 2nd order condition (tst = %g)\n', tst);
      end
   end

   tst = bE'*cI - 0.5;
   if (abs(tst) > tol)
      failed = true;
      failC = true;
      if (reportL>1)
         fprintf('    ARK method fails 2nd order coupling condition A1 (tst = %g)\n', tst);
      end
   end

   tst = bI'*cE - 0.5;
   if (abs(tst) > tol)
      failed = true;
      failC = true;
      if (reportL>1)
         fprintf('    ARK method fails 2nd order coupling condition A2 (tst = %g)\n', tst);
      end
   end

   if (reportL>1)
      if (~failE),  fprintf('  ERK method passes order 2 conditions\n'); end
      if (~failI),  fprintf('  DIRK method passes order 2 conditions\n'); end
      if (~failC),  fprintf('  ARK method pass order 2 coupling conditions\n'); end
   end
   if (~failed)
      q = 2;
   end
   if (~failE)
      qE = max(qE,2);
   end
   if (~failI)
      qI = max(qI,2);
   end
end


% check for third order
if (~failed || ~failE || ~failI)

   tst = bE'*(cE.*cE) - (1/3);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 3nd order condition A (tst = %g)\n', tst);
      end
   end

   tst = bI'*(cI.*cI) - (1/3);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 3rd order condition A (tst = %g)\n', tst);
      end
   end

   tst = bE'*(AE*cE) - (1/6);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 3nd order condition B (tst = %g)\n', tst);
      end
   end

   tst = bI'*(AI*cI) - (1/6);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 3nd order condition B (tst = %g)\n', tst);
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for ic2 = 1:2
            c2 = cs{ic2};
            tst = b'*(c1.*c2) - (1/3);
            counter = counter+1;
            if (abs(tst) > tol)
               failed = true;
               failC = true;
               if (reportL>1)
                  fprintf('    ARK method fails 3rd order coupling condition A%i (tst = %g)\n', counter, tst);
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA = 1:2
         A = As{iA};
         for ic = 1:2
            c = cs{ic};
            tst = b'*(A*c) - (1/6);
            counter = counter+1;
            if (abs(tst) > tol)
               failed = true;
               failC = true;
               if (reportL>1)
                  fprintf('    ARK method fails 3rd order coupling condition B%i (tst = %g)\n', counter, tst);
               end
            end
         end
      end
   end

   if (reportL>1)
      if (~failE),  fprintf('  ERK method passes order 3 conditions\n'); end
      if (~failI),  fprintf('  DIRK method passes order 3 conditions\n'); end
      if (~failC),  fprintf('  ARK method passes order 3 coupling conditions\n'); end
   end
   if (~failed)
      q = 3;
   end
   if (~failE)
      qE = max(qE,3);
   end
   if (~failI)
      qI = max(qI,3);
   end
end




% check for fourth order
if (~failed || ~failE || ~failI)

   tst = bE'*(cE.*cE.*cE) - (1/4);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 4th order condition A (tst = %g)\n', tst);
      end
   end

   tst = bI'*(cI.*cI.*cI) - (1/4);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 4th order condition A (tst = %g)\n', tst);
      end
   end

   tst = (bE.*cE)'*(AE*cE) - (1/8);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 4th order condition B (tst = %g)\n', tst);
      end
   end

   tst = (bI.*cI)'*(AI*cI) - (1/8);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 4th order condition B (tst = %g)\n', tst);
      end
   end

   tst = bE'*AE*(cE.*cE) - (1/12);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 4th order condition C (tst = %g)\n', tst);
      end
   end

   tst = bI'*AI*(cI.*cI) - (1/12);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 4th order condition C (tst = %g)\n', tst);
      end
   end

   tst = bE'*AE*AE*cE - (1/24);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 4th order condition D (tst = %g)\n', tst);
      end
   end

   tst = bI'*AI*AI*cI - (1/24);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 4th order condition D (tst = %g)\n', tst);
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for ic2 = 1:2
            c2 = cs{ic2};
            for ic3 = 1:2
               c3 = cs{ic3};
               tst = b'*(c1.*c2.*c3) - (1/4);
               counter = counter+1;
               if (abs(tst) > tol)
                  failed = true;
                  failC = true;
                  if (reportL>1)
                     fprintf('    ARK method fails 4th order coupling condition A%i (tst = %g)\n', counter, tst);
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for iA = 1:2
            A = As{iA};
            for ic2 = 1:2
               c2 = cs{ic2};
               tst = (b.*c1)'*(A*c2) - (1/8);
               counter = counter+1;
               if (abs(tst) > tol)
                  failed = true;
                  failC = true;
                  if (reportL>1)
                     fprintf('    ARK method fails 4th order coupling condition B%i (tst = %g)\n', counter, tst);
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA = 1:2
         A = As{iA};
         for ic1 = 1:2
            c1 = cs{ic1};
            for ic2 = 1:2
               c2 = cs{ic2};
               tst = b'*A*(c1.*c2) - (1/12);
               counter = counter+1;
               if (abs(tst) > tol)
                  failed = true;
                  failC = true;
                  if (reportL>1)
                     fprintf('    ARK method fails 4th order coupling condition C%i (tst = %g)\n', counter, tst);
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for iA2 = 1:2
            A2 = As{iA2};
            for ic = 1:2
               c = cs{ic};
               tst = b'*A1*A2*c - (1/24);
               counter = counter+1;
               if (abs(tst) > tol)
                  failed = true;
                  failC = true;
                  if (reportL>1)
                     fprintf('    ARK method fails 4th order coupling condition D%i (tst = %g)\n', counter, tst);
                  end
               end
            end
         end
      end
   end

   if (reportL>1)
      if (~failE),  fprintf('  ERK method passes order 4 conditions\n'); end
      if (~failI),  fprintf('  DIRK method passes order 4 conditions\n'); end
      if (~failC),  fprintf('  ARK method passes order 4 coupling conditions\n'); end
   end
   if (~failed)
      q = 4;
   end
   if (~failE)
      qE = max(qE,4);
   end
   if (~failI)
      qI = max(qI,4);
   end
end




% check for fifth order
if (~failed || ~failE || ~failI)

   tst = bE'*(cE.*cE.*cE.*cE) - (1/5);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 5th order condition A (tst = %g)\n', tst);
      end
   end

   tst = bI'*(cI.*cI.*cI.*cI) - (1/5);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 5th order condition A (tst = %g)\n', tst);
      end
   end

   tst = (bE.*cE.*cE)'*(AE*cE) - (1/10);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 5th order condition B (tst = %g)\n', tst);
      end
   end

   tst = (bI.*cI.*cI)'*(AI*cI) - (1/10);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 5th order condition B (tst = %g)\n', tst);
      end
   end

   tst = bE'*((AE*cE).*(AE*cE)) - (1/20);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 5th order condition C (tst = %g)\n', tst);
      end
   end

   tst = bI'*((AI*cI).*(AI*cI)) - (1/20);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    DIRK method fails 5th order condition C (tst = %g)\n', tst);
      end
   end

   tst = (bE.*cE)'*AE*(cE.*cE) - (1/15);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 5th order condition D (tst = %g)\n', tst);
      end
   end

   tst = (bI.*cI)'*AI*(cI.*cI) - (1/15);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 5th order condition D (tst = %g)\n', tst);
      end
   end

   tst = bE'*AE*(cE.*cE.*cE) - (1/20);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 5th order condition E (tst = %g)\n', tst);
      end
   end

   tst = bI'*AI*(cI.*cI.*cI) - (1/20);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    DIRK method fails 5th order condition E (tst = %g)\n', tst);
      end
   end

   tst = (bE.*cE)'*AE*AE*cE - (1/30);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 5th order condition F (tst = %g)\n', tst);
      end
   end

   tst = (bI.*cI)'*AI*AI*cI - (1/30);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    DIRK method fails 5th order condition F (tst = %g)\n', tst);
      end
   end

   tst = -(1/40);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + bE(i)*AE(i,j)*cE(j)*AE(j,k)*cE(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 5th order condition G (tst = %g)\n', tst);
      end
   end

   tst = -(1/40);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + bI(i)*AI(i,j)*cI(j)*AI(j,k)*cI(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    DIRK method fails 5th order condition G (tst = %g)\n', tst);
      end
   end

   tst = bE'*AE*AE*(cE.*cE) - (1/60);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 5th order condition H (tst = %g)\n', tst);
      end
   end

   tst = bI'*AI*AI*(cI.*cI) - (1/60);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    DIRK method fails 5th order condition H (tst = %g)\n', tst);
      end
   end

   tst = bE'*AE*AE*AE*cE - (1/120);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 5th order condition I (tst = %g)\n', tst);
      end
   end

   tst = bI'*AI*AI*AI*cI - (1/120);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    DIRK method fails 5th order condition I (tst = %g)\n', tst);
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for ic2 = 1:2
            c2 = cs{ic2};
            for ic3 = 1:2
               c3 = cs{ic3};
               for ic4 = 1:2
                  c4 = cs{ic4};
                  tst = b'*(c1.*c2.*c3.*c4) - (1/5);
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK method fails 5th order coupling condition A%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for ic2 = 1:2
            c2 = cs{ic2};
            for iA = 1:2
               A = As{iA};
               for ic3 = 1:2
                  c3 = cs{ic3};
                  tst = (b.*c1.*c2)'*(A*c3) - 1/10;
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK method fails 5th order coupling condition B%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for ic1 = 1:2
            c1 = cs{ic1};
            for iA2 = 1:2
               A2 = As{iA2};
               for ic2 = 1:2
                  c2 = cs{ic2};
                  tst = b'*((A1*c1).*(A2*c2)) - (1/20);
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK method fails 5th order coupling condition C%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for iA = 1:2
            A = As{iA};
            for ic2 = 1:2
               c2 = cs{ic2};
               for ic3 = 1:2
                  c3 = cs{ic3};
                  tst = (b.*c1)'*A*(c2.*c3) - (1/15);
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK method fails 5th order coupling condition D%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA = 1:2
         A = As{iA};
         for ic1 = 1:2
            c1 = cs{ic1};
            for ic2 = 1:2
               c2 = cs{ic2};
               for ic3 = 1:2
                  c3 = cs{ic3};
                  tst = b'*A*(c1.*c2.*c3) - (1/20);
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK method fails 5th order coupling condition E%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for iA1 = 1:2
            A1 = As{iA1};
            for iA2 = 1:2
               A2 = As{iA2};
               for ic2 = 1:2
                  c2 = cs{ic2};
                  tst = (b.*c1)'*A1*A2*c2 - (1/30);
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK method fails 5th order coupling condition F%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for ic1 = 1:2
            c1 = cs{ic1};
            for iA2 = 1:2
               A2 = As{iA2};
               for ic2 = 1:2
                  c2 = cs{ic2};
                  tst = -(1/40);
                  for i=1:s
                     for j=1:s
                        for k=1:s
                           tst = tst + b(i)*A1(i,j)*c1(j)*A2(j,k)*c2(k);
                        end
                     end
                  end
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK method fails 5th order coupling condition G%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for iA2 = 1:2
            A2 = As{iA2};
            for ic1 = 1:2
               c1 = cs{ic1};
               for ic2 = 1:2
                  c2 = cs{ic2};
                  tst = b'*A1*A2*(c1.*c2) - (1/60);
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK method fails 5th order coupling condition H%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for iA2 = 1:2
            A2 = As{iA2};
            for iA3 = 1:2
               A3 = As{iA3};
               for ic = 1:2
                  c = cs{ic};
                  tst = b'*A1*A2*A3*c - (1/120);
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK method fails 5th order coupling condition I%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end


   if (reportL>1)
      if (~failE),  fprintf('  ERK method passes order 5 conditions\n'); end
      if (~failI),  fprintf('  DIRK method passes order 5 conditions\n'); end
      if (~failC),  fprintf('  ARK method passes order 5 coupling conditions\n'); end
   end
   if (~failed)
      q = 5;
   end
   if (~failE)
      qE = max(qE,5);
   end
   if (~failI)
      qI = max(qI,5);
   end
end




% check for sixth order
if (~failed || ~failE || ~failI)

   tst = bE'*(cE.*cE.*cE.*cE.*cE) - (1/6);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition A (tst = %g)\n', tst);
      end
   end

   tst = bI'*(cI.*cI.*cI.*cI.*cI) - (1/6);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition A (tst = %g)\n', tst);
      end
   end

   tst = (bE.*cE.*cE.*cE)'*AE*cE - (1/12);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition B (tst = %g)\n', tst);
      end
   end

   tst = (bI.*cI.*cI.*cI)'*AI*cI - (1/12);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition B (tst = %g)\n', tst);
      end
   end

   tst = -(1/24);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + bE(i)*cE(i)*AE(i,j)*cE(j)*AE(j,k)*cE(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition C (tst = %g)\n', tst);
      end
   end

   tst = -(1/24);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + bI(i)*cI(i)*AI(i,j)*cI(j)*AI(j,k)*cI(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition C (tst = %g)\n', tst);
      end
   end

   tst = (bE.*cE.*cE)'*AE*(cE.*cE) - (1/18);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition D (tst = %g)\n', tst);
      end
   end

   tst = (bI.*cI.*cI)'*AI*(cI.*cI) - (1/18);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition D (tst = %g)\n', tst);
      end
   end

   tst = -(1/36);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + bE(i)*AE(i,j)*cE(j)*cE(j)*AE(i,k)*cE(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition E (tst = %g)\n', tst);
      end
   end

   tst = -(1/36);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + bI(i)*AI(i,j)*cI(j)*cI(j)*AI(i,k)*cI(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition E (tst = %g)\n', tst);
      end
   end

   tst = (bE.*cE)'*AE*(cE.*cE.*cE) - (1/24);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition F (tst = %g)\n', tst);
      end
   end

   tst = (bI.*cI)'*AI*(cI.*cI.*cI) - (1/24);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition F (tst = %g)\n', tst);
      end
   end

   tst = bE'*AE*(cE.*cE.*cE.*cE) - (1/30);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition G (tst = %g)\n', tst);
      end
   end

   tst = bI'*AI*(cI.*cI.*cI.*cI) - (1/30);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition G (tst = %g)\n', tst);
      end
   end

   tst = (bE.*cE.*cE)'*AE*AE*cE - (1/36);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition H (tst = %g)\n', tst);
      end
   end

   tst = (bI.*cI.*cI)'*AI*AI*cI - (1/36);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition H (tst = %g)\n', tst);
      end
   end

   tst = -(1/72);
   for i=1:s
      for j=1:s
         for k=1:s
            for l=1:s
               tst = tst + bE(i)*AE(i,j)*AE(j,k)*cE(k)*AE(i,l)*cE(l);
            end
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition I (tst = %g)\n', tst);
      end
   end

   tst = -(1/72);
   for i=1:s
      for j=1:s
         for k=1:s
            for l=1:s
               tst = tst + bI(i)*AI(i,j)*AI(j,k)*cI(k)*AI(i,l)*cI(l);
            end
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition I (tst = %g)\n', tst);
      end
   end

   tst = -(1/48);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + bE(i)*cE(i)*AE(i,j)*cE(j)*AE(j,k)*cE(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition J (tst = %g)\n', tst);
      end
   end

   tst = -(1/48);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + bI(i)*cI(i)*AI(i,j)*cI(j)*AI(j,k)*cI(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition J (tst = %g)\n', tst);
      end
   end

   tst = -(1/60);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + bE(i)*AE(i,j)*cE(j)*cE(j)*AE(j,k)*cE(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition K (tst = %g)\n', tst);
      end
   end

   tst = -(1/60);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + bI(i)*AI(i,j)*cI(j)*cI(j)*AI(j,k)*cI(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition K (tst = %g)\n', tst);
      end
   end

   tst = -(1/120);
   for i=1:s
      for j=1:s
         for k=1:s
            for l=1:s
               tst = tst + bE(i)*AE(i,j)*AE(j,k)*cE(k)*AE(j,l)*cE(l);
            end
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition L (tst = %g)\n', tst);
      end
   end

   tst = -(1/120);
   for i=1:s
      for j=1:s
         for k=1:s
            for l=1:s
               tst = tst + bI(i)*AI(i,j)*AI(j,k)*cI(k)*AI(j,l)*cI(l);
            end
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition L (tst = %g)\n', tst);
      end
   end

   tst = (bE.*cE)'*AE*AE*(cE.*cE) - (1/72);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition M (tst = %g)\n', tst);
      end
   end

   tst = (bI.*cI)'*AI*AI*(cI.*cI) - (1/72);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition M (tst = %g)\n', tst);
      end
   end

   tst = -(1/90);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + bE(i)*AE(i,j)*cE(j)*AE(j,k)*cE(k)*cE(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition N (tst = %g)\n', tst);
      end
   end

   tst = -(1/90);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + bI(i)*AI(i,j)*cI(j)*AI(j,k)*cI(k)*cI(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition N (tst = %g)\n', tst);
      end
   end

   tst = bE'*AE*AE*(cE.*cE.*cE) - (1/120);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition O (tst = %g)\n', tst);
      end
   end

   tst = bI'*AI*AI*(cI.*cI.*cI) - (1/120);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition O (tst = %g)\n', tst);
      end
   end

   tst = (bE.*cE)'*AE*AE*AE*cE - (1/144);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition P (tst = %g)\n', tst);
      end
   end

   tst = (bI.*cI)'*AI*AI*AI*cI - (1/144);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition P (tst = %g)\n', tst);
      end
   end

   tst = -(1/180);
   for i=1:s
      for j=1:s
         for k=1:s
            for l=1:s
               tst = tst + bE(i)*AE(i,j)*cE(j)*AE(j,k)*AE(k,l)*cE(l);
            end
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition Q (tst = %g)\n', tst);
      end
   end

   tst = -(1/180);
   for i=1:s
      for j=1:s
         for k=1:s
            for l=1:s
               tst = tst + bI(i)*AI(i,j)*cI(j)*AI(j,k)*AI(k,l)*cI(l);
            end
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition Q (tst = %g)\n', tst);
      end
   end

   tst = -(1/240);
   for i=1:s
      for j=1:s
         for k=1:s
            for l=1:s
               tst = tst + bE(i)*AE(i,j)*AE(j,k)*cE(k)*AE(k,l)*cE(l);
            end
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition R (tst = %g)\n', tst);
      end
   end

   tst = -(1/240);
   for i=1:s
      for j=1:s
         for k=1:s
            for l=1:s
               tst = tst + bI(i)*AI(i,j)*AI(j,k)*cI(k)*AI(k,l)*cI(l);
            end
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition R (tst = %g)\n', tst);
      end
   end

   tst = bE'*AE*AE*AE*(cE.*cE) - (1/360);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition S (tst = %g)\n', tst);
      end
   end

   tst = bI'*AI*AI*AI*(cI.*cI) - (1/360);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition S (tst = %g)\n', tst);
      end
   end

   tst = bE'*AE*AE*AE*AE*cE - (1/720);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK method fails 6th order condition T (tst = %g)\n', tst);
      end
   end

   tst = bI'*AI*AI*AI*AI*cI - (1/720);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK method fails 6th order condition T (tst = %g)\n', tst);
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for ic2 = 1:2
            c2 = cs{ic2};
            for ic3 = 1:2
               c3 = cs{ic3};
               for ic4 = 1:2
                  c4 = cs{ic4};
                  for ic5 = 1:2
                     c5 = cs{ic5};
                     tst = b'*(c1.*c2.*c3.*c4.*c5) - (1/5);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition A%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for ic2 = 1:2
            c2 = cs{ic2};
            for ic3 = 1:2
               c3 = cs{ic3};
               for iA = 1:2
                  A = As{iA};
                  for ic4 = 1:2
                     c4 = cs{ic4};
                     tst = (b.*c1.*c2.*c3)'*A*c4 - (1/12);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition B%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for iA1 = 1:2
            A1 = As{iA1};
            for ic2 = 1:2
               c2 = cs{ic2};
               for iA2 = 1:2
                  A2 = As{iA2};
                  for ic3 = 1:2
                     c3 = cs{ic3};
                     tst = -(1/24);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              tst = tst + b(i)*c1(i)*A1(i,j)*c2(j)*A2(j,k)*c3(k);
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition C%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for ic2 = 1:2
            c2 = cs{ic2};
            for iA = 1:2
               A = As{iA};
               for ic3 = 1:2
                  c3 = cs{ic3};
                  for ic4 = 1:2
                     c4 = cs{ic4};
                     tst = (b.*c1.*c2)'*A*(c3.*c4) - (1/18);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition D%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end


   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for ic1 = 1:2
            c1 = cs{ic1};
            for ic2 = 1:2
               c2 = cs{ic2};
               for iA2 = 1:2
                  A2 = As{iA2};
                  for ic3 = 1:2
                     c3 = cs{ic3};
                     tst = -(1/36);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              tst = tst + b(i)*A1(i,j)*c1(j)*c2(j)*A2(i,k)*c3(k);
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition E%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for iA = 1:2
            A = As{iA};
            for ic2 = 1:2
               c2 = cs{ic2};
               for ic3 = 1:2
                  c3 = cs{ic3};
                  for ic4 = 1:2
                     c4 = cs{ic4};
                     tst = (b.*c1)'*A*(c2.*c3.*c4) - (1/24);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition F%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA = 1:2
         A = As{iA};
         for ic1 = 1:2
            c1 = cs{ic1};
            for ic2 = 1:2
               c2 = cs{ic2};
               for ic3 = 1:2
                  c3 = cs{ic3};
                  for ic4 = 1:2
                     c4 = cs{ic4};
                     tst = b'*A*(c1.*c2.*c3.*c4) - (1/30);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition G%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for ic2 = 1:2
            c2 = cs{ic2};
            for iA1 = 1:2
               A1 = As{iA};
               for iA2 = 1:2
                  A2 = As{iA2};
                  for ic3 = 1:2
                     c3 = cs{ic3};
                     tst = (b.*c1.*c2)'*A1*A2*c3 - (1/36);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition H%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for iA2 = 1:2
            A2 = As{iA2};
            for ic1 = 1:2
               c1 = cs{ic1};
               for iA3 = 1:2
                  A3 = As{iA3};
                  for ic2 = 1:2
                     c2 = cs{ic2};
                     tst = -(1/72);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              for l=1:s
                                 tst = tst + b(i)*A1(i,j)*A2(j,k)*c1(k)*A3(i,l)*c2(l);
                              end
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition I%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for iA1 = 1:2
            A1 = As{iA1};
            for ic2 = 1:2
               c2 = cs{ic2};
               for iA2 = 1:2
                  A2 = As{iA2};
                  for ic3 = 1:2
                     c3 = cs{ic3};
                     tst = -(1/48);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              tst = tst + b(i)*c1(i)*A1(i,j)*c2(j)*A2(j,k)*c3(k);
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition J%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for ic1 = 1:2
            c1 = cs{ic1};
            for ic2 = 1:2
               c2 = cs{ic2};
               for iA2 = 1:2
                  A2 = As{iA2};
                  for ic3 = 1:2
                     c3 = cs{ic3};
                     tst = -(1/60);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              tst = tst + b(i)*A1(i,j)*c1(j)*c2(j)*A2(j,k)*c3(k);
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition K%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for iA2 = 1:2
            A2 = As{iA2};
            for ic1 = 1:2
               c1 = cs{ic1};
               for iA3 = 1:2
                  A3 = As{iA3};
                  for ic2 = 1:2
                     c2 = cs{ic2};
                     tst = -(1/120);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              for l=1:s
                                 tst = tst + b(i)*A1(i,j)*A2(j,k)*c1(k)*A3(j,l)*c2(l);
                              end
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition L%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for iA1 = 1:2
            A1 = As{iA};
            for iA2 = 1:2
               A2 = As{iA2};
               for ic2 = 1:2
                  c2 = cs{ic2};
                  for ic3 = 1:2
                     c3 = cs{ic3};
                     tst = (b.*c1)'*A1*A2*(c2.*c3) - (1/72);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition M%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for ic1 = 1:2
            c1 = cs{ic1};
            for iA2 = 1:2
               A2 = As{iA2};
               for ic2 = 1:2
                  c2 = cs{ic2};
                  for ic3 = 1:2
                     c3 = cs{ic3};
                     tst = -(1/90);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              tst = tst + b(i)*A1(i,j)*c1(j)*A2(j,k)*c2(k)*c3(k);
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition N%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA};
         for iA2 = 1:2
            A2 = As{iA2};
            for ic1 = 1:2
               c1 = cs{ic1};
               for ic2 = 1:2
                  c2 = cs{ic2};
                  for ic3 = 1:2
                     c3 = cs{ic3};
                     tst = b'*A1*A2*(c1.*c2.*c3) - (1/120);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition O%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for iA1 = 1:2
            A1 = As{iA};
            for iA2 = 1:2
               A2 = As{iA2};
               for iA3 = 1:2
                  A3 = As{iA3};
                  for ic2 = 1:2
                     c2 = cs{ic2};
                     tst = (b.*c1)'*A1*A2*A3*c2 - (1/144);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition P%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for ic1 = 1:2
            c1 = cs{ic1};
            for iA2 = 1:2
               A2 = As{iA2};
               for iA3 = 1:2
                  A3 = As{iA3};
                  for ic2 = 1:2
                     c2 = cs{ic2};
                     tst = -(1/180);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              for l=1:s
                                 tst = tst + b(i)*A1(i,j)*c1(j)*A2(j,k)*A3(k,l)*c2(l);
                              end
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition Q%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for ic1 = 1:2
            c1 = cs{ic1};
            for iA2 = 1:2
               A2 = As{iA2};
               for iA3 = 1:2
                  A3 = As{iA3};
                  for ic2 = 1:2
                     c2 = cs{ic2};
                     tst = -(1/240);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              for l=1:s
                                 tst = tst + b(i)*A1(i,j)*A2(j,k)*c1(k)*A3(k,l)*c2(l);
                              end
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition R%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA};
         for iA2 = 1:2
            A2 = As{iA2};
            for iA3 = 1:2
               A3 = As{iA3};
               for ic1 = 1:2
                  c1 = cs{ic1};
                  for ic2 = 1:2
                     c2 = cs{ic2};
                     tst = b'*A1*A2*A3*(c1.*c2) - (1/360);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition S%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA};
         for iA2 = 1:2
            A2 = As{iA2};
            for iA3 = 1:2
               A3 = As{iA3};
               for iA4 = 1:2
                  A4 = As{iA4};
                  for ic = 1:2
                     c = cs{ic};
                     tst = b'*A1*A2*A3*A4*c - (1/720);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK method fails 6th order coupling condition T%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end


   if (reportL>1)
      if (~failE),  fprintf('  ERK method passes order 6 conditions\n'); end
      if (~failI),  fprintf('  DIRK method passes order 6 conditions\n'); end
      if (~failC),  fprintf('  ARK method passes order 6 coupling conditions\n'); end
   end
   if (~failed)
      q = 6;
   end
   if (~failE)
      qE = max(qE,6);
   end
   if (~failI)
      qI = max(qE,6);
   end
end





% initialize failure flags, order of accuracy
failed = false;
failE = false;
failI = false;
failC = false;
p = -1;

% convert dE and dI to column vectors for these tests
dE = reshape(dE,length(dE),1);
dI = reshape(dI,length(dI),1);

% create cell arrays to streamline pairwise testing
bs = {dE,dI};
As = {AE,AI};
cs = {cE,cI};

% ensure row sum condition c(i) = sum(A(i,:)) is satisfied for each method
[s,tmp] = size(AE);
for i=1:s
   tst = cE(i) - sum(AE(i,:));
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails row sum condition (i = %i, tst = %g)\n', i, tst);
      end
   end
   tst = cI(i) - sum(AI(i,:));
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails row sum condition (i = %i, tst = %g)\n', i, tst);
      end
   end
end
if (~failed)
   p = 0;
end
if (~failE)
   pE = max(pE,0);
end
if (~failI)
   pI = max(pI,0);
end

% check for first order embedding
if (~failed || ~failE || ~failI)
   tst = sum(dE)-1;
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 1st order condition (tst = %g)\n', tst);
      end
   end

   tst = sum(dI)-1;
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 1st order condition (tst = %g)\n', tst);
      end
   end
      if (reportL>1)
         if (~failE), fprintf('  ERK embedding passes order 1 conditions\n'); end
         if (~failI), fprintf('  DIRK embedding passes order 1 conditions\n'); end
      end
      if (~failed)
         p = 1;
      end
      if (~failE)
         pE = max(pE,1);
      end
      if (~failI)
         pI = max(pI,1);
      end
end

% check for second order embedding
if (~failed || ~failE || ~failI)

   tst = dE'*cE - 0.5;
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 2nd order condition (tst = %g)\n', tst);
      end
   end

   tst = dI'*cI - 0.5;
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 2nd order condition (tst = %g)\n', tst);
      end
   end

   tst = dE'*cI - 0.5;
   if (abs(tst) > tol)
      failed = true;
      failC = true;
      if (reportL>1)
         fprintf('    ARK embedding fails 2nd order coupling condition A1 (tst = %g)\n', tst);
      end
   end

   tst = dI'*cE - 0.5;
   if (abs(tst) > tol)
      failed = true;
      failC = true;
      if (reportL>1)
         fprintf('    ARK embedding fails 2nd order coupling condition A2 (tst = %g)\n', tst);
      end
   end

   if (reportL>1)
      if (~failE),  fprintf('  ERK embedding passes order 2 conditions\n'); end
      if (~failI),  fprintf('  DIRK embedding passes order 2 conditions\n'); end
      if (~failC),  fprintf('  ARK embedding passes order 2 coupling conditions\n'); end
   end
   if (~failed)
      p = 2;
   end
   if (~failE)
      pE = max(pE,2);
   end
   if (~failI)
      pI = max(pI,2);
   end
end


% check for third order embedding
if (~failed || ~failE || ~failI)

   tst = dE'*(cE.*cE) - (1/3);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 3nd order condition A (tst = %g)\n', tst);
      end
   end

   tst = dI'*(cI.*cI) - (1/3);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 3rd order condition A (tst = %g)\n', tst);
      end
   end

   tst = dE'*(AE*cE) - (1/6);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 3nd order condition B (tst = %g)\n', tst);
      end
   end

   tst = dI'*(AI*cI) - (1/6);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 3nd order condition B (tst = %g)\n', tst);
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for ic2 = 1:2
            c2 = cs{ic2};
            tst = b'*(c1.*c2) - (1/3);
            counter = counter+1;
            if (abs(tst) > tol)
               failed = true;
               failC = true;
               if (reportL>1)
                  fprintf('    ARK embedding fails 3rd order coupling condition A%i (tst = %g)\n', counter, tst);
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA = 1:2
         A = As{iA};
         for ic = 1:2
            c = cs{ic};
            tst = b'*(A*c) - (1/6);
            counter = counter+1;
            if (abs(tst) > tol)
               failed = true;
               failC = true;
               if (reportL>1)
                  fprintf('    ARK embedding fails 3rd order coupling condition B%i (tst = %g)\n', counter, tst);
               end
            end
         end
      end
   end

   if (reportL>1)
      if (~failE),  fprintf('  ERK embedding passes order 3 conditions\n'); end
      if (~failI),  fprintf('  DIRK embedding passes order 3 conditions\n'); end
      if (~failC),  fprintf('  ARK embedding passes order 3 coupling conditions\n'); end
   end
   if (~failed)
      p = 3;
   end
   if (~failE)
      pE = max(pE,3);
   end
   if (~failI)
      pI = max(pI,3);
   end
end




% check for fourth order embedding
if (~failed || ~failE || ~failI)

   tst = dE'*(cE.*cE.*cE) - (1/4);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 4th order condition A (tst = %g)\n', tst);
      end
   end

   tst = dI'*(cI.*cI.*cI) - (1/4);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 4th order condition A (tst = %g)\n', tst);
      end
   end

   tst = (dE.*cE)'*(AE*cE) - (1/8);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 4th order condition B (tst = %g)\n', tst);
      end
   end

   tst = (dI.*cI)'*(AI*cI) - (1/8);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 4th order condition B (tst = %g)\n', tst);
      end
   end

   tst = dE'*AE*(cE.*cE) - (1/12);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 4th order condition C (tst = %g)\n', tst);
      end
   end

   tst = dI'*AI*(cI.*cI) - (1/12);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 4th order condition C (tst = %g)\n', tst);
      end
   end

   tst = dE'*AE*AE*cE - (1/24);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 4th order condition D (tst = %g)\n', tst);
      end
   end

   tst = dI'*AI*AI*cI - (1/24);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 4th order condition D (tst = %g)\n', tst);
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for ic2 = 1:2
            c2 = cs{ic2};
            for ic3 = 1:2
               c3 = cs{ic3};
               tst = b'*(c1.*c2.*c3) - (1/4);
               counter = counter+1;
               if (abs(tst) > tol)
                  failed = true;
                  failC = true;
                  if (reportL>1)
                     fprintf('    ARK embedding fails 4th order coupling condition A%i (tst = %g)\n', counter, tst);
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for iA = 1:2
            A = As{iA};
            for ic2 = 1:2
               c2 = cs{ic2};
               tst = (b.*c1)'*(A*c2) - (1/8);
               counter = counter+1;
               if (abs(tst) > tol)
                  failed = true;
                  failC = true;
                  if (reportL>1)
                     fprintf('    ARK embedding fails 4th order coupling condition B%i (tst = %g)\n', counter, tst);
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA = 1:2
         A = As{iA};
         for ic1 = 1:2
            c1 = cs{ic1};
            for ic2 = 1:2
               c2 = cs{ic2};
               tst = b'*A*(c1.*c2) - (1/12);
               counter = counter+1;
               if (abs(tst) > tol)
                  failed = true;
                  failC = true;
                  if (reportL>1)
                     fprintf('    ARK embedding fails 4th order coupling condition C%i (tst = %g)\n', counter, tst);
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for iA2 = 1:2
            A2 = As{iA2};
            for ic = 1:2
               c = cs{ic};
               tst = b'*A1*A2*c - (1/24);
               counter = counter+1;
               if (abs(tst) > tol)
                  failed = true;
                  failC = true;
                  if (reportL>1)
                     fprintf('    ARK embedding fails 4th order coupling condition D%i (tst = %g)\n', counter, tst);
                  end
               end
            end
         end
      end
   end

   if (reportL>1)
      if (~failE),  fprintf('  ERK embedding passes order 4 conditions\n'); end
      if (~failI),  fprintf('  DIRK embedding passes order 4 conditions\n'); end
      if (~failC),  fprintf('  ARK embedding passes order 4 coupling conditions\n'); end
   end
   if (~failed)
      p = 4;
   end
   if (~failE)
      pE = max(pE,4);
   end
   if (~failI)
      pI = max(pI,4);
   end
end




% check for fifth order embedding
if (~failed || ~failE || ~failI)

   tst = dE'*(cE.*cE.*cE.*cE) - (1/5);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 5th order condition A (tst = %g)\n', tst);
      end
   end

   tst = dI'*(cI.*cI.*cI.*cI) - (1/5);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 5th order condition A (tst = %g)\n', tst);
      end
   end

   tst = (dE.*cE.*cE)'*(AE*cE) - (1/10);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 5th order condition B (tst = %g)\n', tst);
      end
   end

   tst = (dI.*cI.*cI)'*(AI*cI) - (1/10);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 5th order condition B (tst = %g)\n', tst);
      end
   end

   tst = dE'*((AE*cE).*(AE*cE)) - (1/20);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 5th order condition C (tst = %g)\n', tst);
      end
   end

   tst = dI'*((AI*cI).*(AI*cI)) - (1/20);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 5th order condition C (tst = %g)\n', tst);
      end
   end

   tst = (dE.*cE)'*AE*(cE.*cE) - (1/15);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 5th order condition D (tst = %g)\n', tst);
      end
   end

   tst = (dI.*cI)'*AI*(cI.*cI) - (1/15);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 5th order condition D (tst = %g)\n', tst);
      end
   end

   tst = dE'*AE*(cE.*cE.*cE) - (1/20);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 5th order condition E (tst = %g)\n', tst);
      end
   end

   tst = dI'*AI*(cI.*cI.*cI) - (1/20);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 5th order condition E (tst = %g)\n', tst);
      end
   end

   tst = (dE.*cE)'*AE*AE*cE - (1/30);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 5th order condition F (tst = %g)\n', tst);
      end
   end

   tst = (dI.*cI)'*AI*AI*cI - (1/30);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 5th order condition F (tst = %g)\n', tst);
      end
   end

   tst = -(1/40);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + dE(i)*AE(i,j)*cE(j)*AE(j,k)*cE(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 5th order condition G (tst = %g)\n', tst);
      end
   end

   tst = -(1/40);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + dI(i)*AI(i,j)*cI(j)*AI(j,k)*cI(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 5th order condition G (tst = %g)\n', tst);
      end
   end

   tst = dE'*AE*AE*(cE.*cE) - (1/60);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 5th order condition H (tst = %g)\n', tst);
      end
   end

   tst = dI'*AI*AI*(cI.*cI) - (1/60);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 5th order condition H (tst = %g)\n', tst);
      end
   end

   tst = dE'*AE*AE*AE*cE - (1/120);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 5th order condition I (tst = %g)\n', tst);
      end
   end

   tst = dI'*AI*AI*AI*cI - (1/120);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 5th order condition I (tst = %g)\n', tst);
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for ic2 = 1:2
            c2 = cs{ic2};
            for ic3 = 1:2
               c3 = cs{ic3};
               for ic4 = 1:2
                  c4 = cs{ic4};
                  tst = b'*(c1.*c2.*c3.*c4) - (1/5);
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK embedding fails 5th order coupling condition A%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for ic2 = 1:2
            c2 = cs{ic2};
            for iA = 1:2
               A = As{iA};
               for ic3 = 1:2
                  c3 = cs{ic3};
                  tst = (b.*c1.*c2)'*(A*c3) - 1/10;
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK embedding fails 5th order coupling condition B%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for ic1 = 1:2
            c1 = cs{ic1};
            for iA2 = 1:2
               A2 = As{iA2};
               for ic2 = 1:2
                  c2 = cs{ic2};
                  tst = b'*((A1*c1).*(A2*c2)) - (1/20);
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK embedding fails 5th order coupling condition C%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for iA = 1:2
            A = As{iA};
            for ic2 = 1:2
               c2 = cs{ic2};
               for ic3 = 1:2
                  c3 = cs{ic3};
                  tst = (b.*c1)'*A*(c2.*c3) - (1/15);
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK embedding fails 5th order coupling condition D%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA = 1:2
         A = As{iA};
         for ic1 = 1:2
            c1 = cs{ic1};
            for ic2 = 1:2
               c2 = cs{ic2};
               for ic3 = 1:2
                  c3 = cs{ic3};
                  tst = b'*A*(c1.*c2.*c3) - (1/20);
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK embedding fails 5th order coupling condition E%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for iA1 = 1:2
            A1 = As{iA1};
            for iA2 = 1:2
               A2 = As{iA2};
               for ic2 = 1:2
                  c2 = cs{ic2};
                  tst = (b.*c1)'*A1*A2*c2 - (1/30);
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK embedding fails 5th order coupling condition F%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for ic1 = 1:2
            c1 = cs{ic1};
            for iA2 = 1:2
               A2 = As{iA2};
               for ic2 = 1:2
                  c2 = cs{ic2};
                  tst = -(1/40);
                  for i=1:s
                     for j=1:s
                        for k=1:s
                           tst = tst + b(i)*A1(i,j)*c1(j)*A2(j,k)*c2(k);
                        end
                     end
                  end
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK embedding fails 5th order coupling condition G%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for iA2 = 1:2
            A2 = As{iA2};
            for ic1 = 1:2
               c1 = cs{ic1};
               for ic2 = 1:2
                  c2 = cs{ic2};
                  tst = b'*A1*A2*(c1.*c2) - (1/60);
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK embedding fails 5th order coupling condition H%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for iA2 = 1:2
            A2 = As{iA2};
            for iA3 = 1:2
               A3 = As{iA3};
               for ic = 1:2
                  c = cs{ic};
                  tst = b'*A1*A2*A3*c - (1/120);
                  counter = counter+1;
                  if (abs(tst) > tol)
                     failed = true;
                     failC = true;
                     if (reportL>1)
                        fprintf('    ARK embedding fails 5th order coupling condition I%i (tst = %g)\n', counter, tst);
                     end
                  end
               end
            end
         end
      end
   end


   if (reportL>1)
      if (~failE),  fprintf('  ERK embedding passes order 5 conditions\n'); end
      if (~failI),  fprintf('  DIRK embedding passes order 5 conditions\n'); end
      if (~failC),  fprintf('  ARK embedding passes order 5 coupling conditions\n'); end
   end
   if (~failed)
      p = 5;
   end
   if (~failE)
      pE = max(pE,5);
   end
   if (~failI)
      pI = max(pI,5);
   end
end




% check for sixth order embedding
if (~failed || ~failE || ~failI)

   tst = dE'*(cE.*cE.*cE.*cE.*cE) - (1/6);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition A (tst = %g)\n', tst);
      end
   end

   tst = dI'*(cI.*cI.*cI.*cI.*cI) - (1/6);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition A (tst = %g)\n', tst);
      end
   end

   tst = (dE.*cE.*cE.*cE)'*AE*cE - (1/12);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition B (tst = %g)\n', tst);
      end
   end

   tst = (dI.*cI.*cI.*cI)'*AI*cI - (1/12);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition B (tst = %g)\n', tst);
      end
   end

   tst = -(1/24);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + dE(i)*cE(i)*AE(i,j)*cE(j)*AE(j,k)*cE(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition C (tst = %g)\n', tst);
      end
   end

   tst = -(1/24);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + dI(i)*cI(i)*AI(i,j)*cI(j)*AI(j,k)*cI(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition C (tst = %g)\n', tst);
      end
   end

   tst = (dE.*cE.*cE)'*AE*(cE.*cE) - (1/18);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition D (tst = %g)\n', tst);
      end
   end

   tst = (dI.*cI.*cI)'*AI*(cI.*cI) - (1/18);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition D (tst = %g)\n', tst);
      end
   end

   tst = -(1/36);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + dE(i)*AE(i,j)*cE(j)*cE(j)*AE(i,k)*cE(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition E (tst = %g)\n', tst);
      end
   end

   tst = -(1/36);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + dI(i)*AI(i,j)*cI(j)*cI(j)*AI(i,k)*cI(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition E (tst = %g)\n', tst);
      end
   end

   tst = (dE.*cE)'*AE*(cE.*cE.*cE) - (1/24);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition F (tst = %g)\n', tst);
      end
   end

   tst = (dI.*cI)'*AI*(cI.*cI.*cI) - (1/24);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition F (tst = %g)\n', tst);
      end
   end

   tst = dE'*AE*(cE.*cE.*cE.*cE) - (1/30);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition G (tst = %g)\n', tst);
      end
   end

   tst = dI'*AI*(cI.*cI.*cI.*cI) - (1/30);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition G (tst = %g)\n', tst);
      end
   end

   tst = (dE.*cE.*cE)'*AE*AE*cE - (1/36);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition H (tst = %g)\n', tst);
      end
   end

   tst = (dI.*cI.*cI)'*AI*AI*cI - (1/36);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition H (tst = %g)\n', tst);
      end
   end

   tst = -(1/72);
   for i=1:s
      for j=1:s
         for k=1:s
            for l=1:s
               tst = tst + dE(i)*AE(i,j)*AE(j,k)*cE(k)*AE(i,l)*cE(l);
            end
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition I (tst = %g)\n', tst);
      end
   end

   tst = -(1/72);
   for i=1:s
      for j=1:s
         for k=1:s
            for l=1:s
               tst = tst + dI(i)*AI(i,j)*AI(j,k)*cI(k)*AI(i,l)*cI(l);
            end
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition I (tst = %g)\n', tst);
      end
   end

   tst = -(1/48);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + dE(i)*cE(i)*AE(i,j)*cE(j)*AE(j,k)*cE(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition J (tst = %g)\n', tst);
      end
   end

   tst = -(1/48);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + dI(i)*cI(i)*AI(i,j)*cI(j)*AI(j,k)*cI(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition J (tst = %g)\n', tst);
      end
   end

   tst = -(1/60);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + dE(i)*AE(i,j)*cE(j)*cE(j)*AE(j,k)*cE(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition K (tst = %g)\n', tst);
      end
   end

   tst = -(1/60);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + dI(i)*AI(i,j)*cI(j)*cI(j)*AI(j,k)*cI(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition K (tst = %g)\n', tst);
      end
   end

   tst = -(1/120);
   for i=1:s
      for j=1:s
         for k=1:s
            for l=1:s
               tst = tst + dE(i)*AE(i,j)*AE(j,k)*cE(k)*AE(j,l)*cE(l);
            end
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition L (tst = %g)\n', tst);
      end
   end

   tst = -(1/120);
   for i=1:s
      for j=1:s
         for k=1:s
            for l=1:s
               tst = tst + dI(i)*AI(i,j)*AI(j,k)*cI(k)*AI(j,l)*cI(l);
            end
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition L (tst = %g)\n', tst);
      end
   end

   tst = (dE.*cE)'*AE*AE*(cE.*cE) - (1/72);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition M (tst = %g)\n', tst);
      end
   end

   tst = (dI.*cI)'*AI*AI*(cI.*cI) - (1/72);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition M (tst = %g)\n', tst);
      end
   end

   tst = -(1/90);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + dE(i)*AE(i,j)*cE(j)*AE(j,k)*cE(k)*cE(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition N (tst = %g)\n', tst);
      end
   end

   tst = -(1/90);
   for i=1:s
      for j=1:s
         for k=1:s
            tst = tst + dI(i)*AI(i,j)*cI(j)*AI(j,k)*cI(k)*cI(k);
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition N (tst = %g)\n', tst);
      end
   end

   tst = dE'*AE*AE*(cE.*cE.*cE) - (1/120);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition O (tst = %g)\n', tst);
      end
   end

   tst = dI'*AI*AI*(cI.*cI.*cI) - (1/120);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition O (tst = %g)\n', tst);
      end
   end

   tst = (dE.*cE)'*AE*AE*AE*cE - (1/144);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition P (tst = %g)\n', tst);
      end
   end

   tst = (dI.*cI)'*AI*AI*AI*cI - (1/144);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition P (tst = %g)\n', tst);
      end
   end

   tst = -(1/180);
   for i=1:s
      for j=1:s
         for k=1:s
            for l=1:s
               tst = tst + dE(i)*AE(i,j)*cE(j)*AE(j,k)*AE(k,l)*cE(l);
            end
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition Q (tst = %g)\n', tst);
      end
   end

   tst = -(1/180);
   for i=1:s
      for j=1:s
         for k=1:s
            for l=1:s
               tst = tst + dI(i)*AI(i,j)*cI(j)*AI(j,k)*AI(k,l)*cI(l);
            end
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition Q (tst = %g)\n', tst);
      end
   end

   tst = -(1/240);
   for i=1:s
      for j=1:s
         for k=1:s
            for l=1:s
               tst = tst + dE(i)*AE(i,j)*AE(j,k)*cE(k)*AE(k,l)*cE(l);
            end
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition R (tst = %g)\n', tst);
      end
   end

   tst = -(1/240);
   for i=1:s
      for j=1:s
         for k=1:s
            for l=1:s
               tst = tst + dI(i)*AI(i,j)*AI(j,k)*cI(k)*AI(k,l)*cI(l);
            end
         end
      end
   end
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition R (tst = %g)\n', tst);
      end
   end

   tst = dE'*AE*AE*AE*(cE.*cE) - (1/360);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition S (tst = %g)\n', tst);
      end
   end

   tst = dI'*AI*AI*AI*(cI.*cI) - (1/360);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition S (tst = %g)\n', tst);
      end
   end

   tst = dE'*AE*AE*AE*AE*cE - (1/720);
   if (abs(tst) > tol)
      failed = true;
      failE = true;
      if (reportL>1)
         fprintf('    ERK embedding fails 6th order condition T (tst = %g)\n', tst);
      end
   end

   tst = dI'*AI*AI*AI*AI*cI - (1/720);
   if (abs(tst) > tol)
      failed = true;
      failI = true;
      if (reportL>1)
         fprintf('    DIRK embedding fails 6th order condition T (tst = %g)\n', tst);
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for ic2 = 1:2
            c2 = cs{ic2};
            for ic3 = 1:2
               c3 = cs{ic3};
               for ic4 = 1:2
                  c4 = cs{ic4};
                  for ic5 = 1:2
                     c5 = cs{ic5};
                     tst = b'*(c1.*c2.*c3.*c4.*c5) - (1/5);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition A%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for ic2 = 1:2
            c2 = cs{ic2};
            for ic3 = 1:2
               c3 = cs{ic3};
               for iA = 1:2
                  A = As{iA};
                  for ic4 = 1:2
                     c4 = cs{ic4};
                     tst = (b.*c1.*c2.*c3)'*A*c4 - (1/12);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition B%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for iA1 = 1:2
            A1 = As{iA1};
            for ic2 = 1:2
               c2 = cs{ic2};
               for iA2 = 1:2
                  A2 = As{iA2};
                  for ic3 = 1:2
                     c3 = cs{ic3};
                     tst = -(1/24);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              tst = tst + b(i)*c1(i)*A1(i,j)*c2(j)*A2(j,k)*c3(k);
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition C%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for ic2 = 1:2
            c2 = cs{ic2};
            for iA = 1:2
               A = As{iA};
               for ic3 = 1:2
                  c3 = cs{ic3};
                  for ic4 = 1:2
                     c4 = cs{ic4};
                     tst = (b.*c1.*c2)'*A*(c3.*c4) - (1/18);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition D%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end


   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for ic1 = 1:2
            c1 = cs{ic1};
            for ic2 = 1:2
               c2 = cs{ic2};
               for iA2 = 1:2
                  A2 = As{iA2};
                  for ic3 = 1:2
                     c3 = cs{ic3};
                     tst = -(1/36);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              tst = tst + b(i)*A1(i,j)*c1(j)*c2(j)*A2(i,k)*c3(k);
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition E%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for iA = 1:2
            A = As{iA};
            for ic2 = 1:2
               c2 = cs{ic2};
               for ic3 = 1:2
                  c3 = cs{ic3};
                  for ic4 = 1:2
                     c4 = cs{ic4};
                     tst = (b.*c1)'*A*(c2.*c3.*c4) - (1/24);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition F%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA = 1:2
         A = As{iA};
         for ic1 = 1:2
            c1 = cs{ic1};
            for ic2 = 1:2
               c2 = cs{ic2};
               for ic3 = 1:2
                  c3 = cs{ic3};
                  for ic4 = 1:2
                     c4 = cs{ic4};
                     tst = b'*A*(c1.*c2.*c3.*c4) - (1/30);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition G%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for ic2 = 1:2
            c2 = cs{ic2};
            for iA1 = 1:2
               A1 = As{iA};
               for iA2 = 1:2
                  A2 = As{iA2};
                  for ic3 = 1:2
                     c3 = cs{ic3};
                     tst = (b.*c1.*c2)'*A1*A2*c3 - (1/36);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition H%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for iA2 = 1:2
            A2 = As{iA2};
            for ic1 = 1:2
               c1 = cs{ic1};
               for iA3 = 1:2
                  A3 = As{iA3};
                  for ic2 = 1:2
                     c2 = cs{ic2};
                     tst = -(1/72);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              for l=1:s
                                 tst = tst + b(i)*A1(i,j)*A2(j,k)*c1(k)*A3(i,l)*c2(l);
                              end
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition I%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for iA1 = 1:2
            A1 = As{iA1};
            for ic2 = 1:2
               c2 = cs{ic2};
               for iA2 = 1:2
                  A2 = As{iA2};
                  for ic3 = 1:2
                     c3 = cs{ic3};
                     tst = -(1/48);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              tst = tst + b(i)*c1(i)*A1(i,j)*c2(j)*A2(j,k)*c3(k);
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition J%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for ic1 = 1:2
            c1 = cs{ic1};
            for ic2 = 1:2
               c2 = cs{ic2};
               for iA2 = 1:2
                  A2 = As{iA2};
                  for ic3 = 1:2
                     c3 = cs{ic3};
                     tst = -(1/60);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              tst = tst + b(i)*A1(i,j)*c1(j)*c2(j)*A2(j,k)*c3(k);
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition K%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for iA2 = 1:2
            A2 = As{iA2};
            for ic1 = 1:2
               c1 = cs{ic1};
               for iA3 = 1:2
                  A3 = As{iA3};
                  for ic2 = 1:2
                     c2 = cs{ic2};
                     tst = -(1/120);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              for l=1:s
                                 tst = tst + b(i)*A1(i,j)*A2(j,k)*c1(k)*A3(j,l)*c2(l);
                              end
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition L%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for iA1 = 1:2
            A1 = As{iA};
            for iA2 = 1:2
               A2 = As{iA2};
               for ic2 = 1:2
                  c2 = cs{ic2};
                  for ic3 = 1:2
                     c3 = cs{ic3};
                     tst = (b.*c1)'*A1*A2*(c2.*c3) - (1/72);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition M%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for ic1 = 1:2
            c1 = cs{ic1};
            for iA2 = 1:2
               A2 = As{iA2};
               for ic2 = 1:2
                  c2 = cs{ic2};
                  for ic3 = 1:2
                     c3 = cs{ic3};
                     tst = -(1/90);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              tst = tst + b(i)*A1(i,j)*c1(j)*A2(j,k)*c2(k)*c3(k);
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition N%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA};
         for iA2 = 1:2
            A2 = As{iA2};
            for ic1 = 1:2
               c1 = cs{ic1};
               for ic2 = 1:2
                  c2 = cs{ic2};
                  for ic3 = 1:2
                     c3 = cs{ic3};
                     tst = b'*A1*A2*(c1.*c2.*c3) - (1/120);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition O%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for ic1 = 1:2
         c1 = cs{ic1};
         for iA1 = 1:2
            A1 = As{iA};
            for iA2 = 1:2
               A2 = As{iA2};
               for iA3 = 1:2
                  A3 = As{iA3};
                  for ic2 = 1:2
                     c2 = cs{ic2};
                     tst = (b.*c1)'*A1*A2*A3*c2 - (1/144);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition P%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for ic1 = 1:2
            c1 = cs{ic1};
            for iA2 = 1:2
               A2 = As{iA2};
               for iA3 = 1:2
                  A3 = As{iA3};
                  for ic2 = 1:2
                     c2 = cs{ic2};
                     tst = -(1/180);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              for l=1:s
                                 tst = tst + b(i)*A1(i,j)*c1(j)*A2(j,k)*A3(k,l)*c2(l);
                              end
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition Q%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA1};
         for ic1 = 1:2
            c1 = cs{ic1};
            for iA2 = 1:2
               A2 = As{iA2};
               for iA3 = 1:2
                  A3 = As{iA3};
                  for ic2 = 1:2
                     c2 = cs{ic2};
                     tst = -(1/240);
                     for i=1:s
                        for j=1:s
                           for k=1:s
                              for l=1:s
                                 tst = tst + b(i)*A1(i,j)*A2(j,k)*c1(k)*A3(k,l)*c2(l);
                              end
                           end
                        end
                     end
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition R%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA};
         for iA2 = 1:2
            A2 = As{iA2};
            for iA3 = 1:2
               A3 = As{iA3};
               for ic1 = 1:2
                  c1 = cs{ic1};
                  for ic2 = 1:2
                     c2 = cs{ic2};
                     tst = b'*A1*A2*A3*(c1.*c2) - (1/360);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition S%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end

   counter=0;
   for ib = 1:2
      b = bs{ib};
      for iA1 = 1:2
         A1 = As{iA};
         for iA2 = 1:2
            A2 = As{iA2};
            for iA3 = 1:2
               A3 = As{iA3};
               for iA4 = 1:2
                  A4 = As{iA4};
                  for ic = 1:2
                     c = cs{ic};
                     tst = b'*A1*A2*A3*A4*c - (1/720);
                     counter = counter+1;
                     if (abs(tst) > tol)
                        failed = true;
                        failC = true;
                        if (reportL>1)
                           fprintf('    ARK embedding fails 6th order coupling condition T%i (tst = %g)\n', counter, tst);
                        end
                     end
                  end
               end
            end
         end
      end
   end


   if (reportL>1)
      if (~failE),  fprintf('  ERK embedding passes order 6 conditions\n'); end
      if (~failI),  fprintf('  DIRK embedding passes order 6 conditions\n'); end
      if (~failC),  fprintf('  ARK embedding passes order 6 coupling conditions\n'); end
   end
   if (~failed)
      p = 6;
   end
   if (~failE)
      pE = max(pE,6);
   end
   if (~failI)
      pI = max(pI,6);
   end
end


% report
if (reportL>0)
   % determine SSP coefficients for component methods and embeddings
   try
      ERK_ssp = ssp_coefficient(double(AE), double(bE));
   catch
      ERK_ssp = -1;
   end
   try
      ERKembed_ssp = ssp_coefficient(double(AE), double(dE));
   catch
      ERKembed_ssp = -1;
   end
   try
      DIRK_ssp = ssp_coefficient(double(AI), double(bI));
   catch
      DIRK_ssp = -1;
   end
   try
      DIRKembed_ssp = ssp_coefficient(double(AI), double(dI));
   catch
      DIRKembed_ssp = -1;
   end

   fprintf('  Overall results:\n');
   fprintf('    ERK:   order = %i,  embed order = %i,  stage order = %i,  stiffly accurate = %i,  SSP coeff = %g,  embed SSP = %g\n', ...
           qE, pE, qsE, expSA, ERK_ssp, ERKembed_ssp);
   fprintf('    DIRK:  order = %i,  embed order = %i,  stage order = %i,  stiffly accurate = %i,  SSP coeff = %g,  embed SSP = %g\n', ...
           qI, pI, qsI, impSA, DIRK_ssp, DIRKembed_ssp);
   fprintf('           method    B/A/L stability = %i/%i/%i\n', BsImp, AsImp, LsImp);
   fprintf('           embedding B/A/L stability = %i/%i/%i\n', BsEmb, AsEmb, LsEmb);
   fprintf('    ARK:   order = %i,  embed order = %i,  stage order = %i\n', q, p, qsA);
end


% generate plots of stability regions
if (doPlot)

   % joint stability region
   fig = figure();
   xl = box(1:2);  yl = box(3:4);
   xax = plot( linspace(xl(1),xl(2),10), zeros(1,10), 'k:'); hold on
   yax = plot( zeros(1,10), linspace(yl(1),yl(2),10), 'k:');
%   stab_region(double(AE), double(bE), box, fig, 'r-', 'ERK');
%   stab_region(double(AE), double(dE), box, fig, 'g-', 'ERK embedding');
%   stab_region(double(AI), double(bI), box, fig, 'b-', 'DIRK');
%   stab_region(double(AI), double(dI), box, fig, 'm-', 'DIRK embedding');
   thetas = [0, 20, 40, 60, 80];
   joint_stab_region(double(AE),double(AI),double(bE),double(bI),thetas,box,fig);
   set(get(get(xax,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
   set(get(get(yax,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
%   legend('Location', 'northwest')
   axis(box)
   xlabel('Re(z)')
   ylabel('Im(z)')
%   title(sprintf('%s joint stability regions, ord %i, emb ord %i',mname,q,p))
%   print(sprintf('%s_stab_regions.png', fname), '-dpng');
%   print(sprintf('%s_stab_regions.eps', fname), '-depsc');
%   savefig(sprintf('%s_stab_regions.fig', fname));
   title(sprintf('ARK joint stability for %s method',mname))
   print(sprintf('%s_joint_stab_region.png', fname), '-dpng');
   print(sprintf('%s_joint_stab_region.eps', fname), '-depsc');
   savefig(sprintf('%s_joint_stab_region.fig', fname));

   % individual stability regions
   fig = figure();
   xl = box(1:2);  yl = box(3:4);
   xax = plot(linspace(xl(1),xl(2),10),zeros(1,10),'k:'); hold on
   yax = plot(zeros(1,10),linspace(yl(1),yl(2),10),'k:');
   stab_region(double(AE),double(bE),box,fig,'r-','ERK method');  % explicit method
   stab_region(double(AI),double(bI),box,fig,'b-','DIRK method');  % implicit method
   stab_region(double(AE),double(dE),box,fig,'g-','ERK embedding');  % explicit embedding
   stab_region(double(AI),double(dI),box,fig,'m-','DIRK embedding');  % implicit embedding
   set(get(get(xax,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
   set(get(get(yax,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
   axis(box)
   xlabel('Re(z)')
   ylabel('Im(z)')
   title(sprintf('%s linear stability regions',mname))
   lgd = legend;
   s = lgd.String;
   for i=5:length(s)-1
     s{i} = '';
   end
   legend(s);
   print(sprintf('%s_stab_regions.png', fname), '-dpng');
   print(sprintf('%s_stab_regions.eps', fname), '-depsc');
   savefig(sprintf('%s_stab_regions.fig', fname));
end


end

% end of function
