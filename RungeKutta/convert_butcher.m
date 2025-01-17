%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2018, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

diary 'Verner-decimal.txt'
output_table('Verner-6-5b-ERK')
output_table('Verner-7-6-ERK')
output_table('Verner-8-7-ERK')
diary off


function output_table(table)
  B = butcher(table);
  s = size(B,2)-1;
  c = B(1:s,1);
  A = B(1:s,2:s+1);
  b = B(s+1,2:s+1);
  d = B(s+2,2:s+1);
  disp(table)
  disp('  ')
  for i=1:s
    for j=1:s
      if (A(i,j) ~= 0)
        fprintf('B->A[%i][%i] = ', i-1, j-1)
        vpa(A(i,j),40)
      end
    end
  end
  disp('  ')
  for i=1:s
    if (b(i) ~= 0)
      fprintf('B->b[%i] = ', i-1)
      vpa(b(i),40)
    end
  end
  disp('  ')
  for i=1:s
    if (d(i) ~= 0)
      fprintf('B->d[%i] = ', i-1)
      vpa(d(i),40)
    end
  end
  disp('  ')
  for i=1:s
    if (c(i) ~= 0)
      fprintf('B->c[%i] = ', i-1)
      vpa(c(i),40)
    end
  end

end


% end of script
