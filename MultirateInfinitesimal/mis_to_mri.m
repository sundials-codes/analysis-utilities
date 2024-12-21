function [c,G] = mis_to_mri(B)
% usage: [c,G] = mis_to_mri(B)
%
% This routine takes as input a "slow" Butcher table (B) and
% generates the cell 'array' of MRI "Gamma" matrices, {G0}
% so that an MIS method with B is equivalent to the MRI-GARK method
% with G.  The slow table abcissae, c, that accompany G are also returned.
%
% This routine checks for validity of the abcissae in the B table:
%     0 = c(1), c(i) <= c(i+1), c(s) <= 1
% and that B is at most diagonally implicit.  If any of these
% conditions are violated this returns an error.
%
% This routine additionally handles the case where B requires padding to
% satisfy the MRI requirement that c(s)=1.
%
% This routine will also add the embedding row to G if an embedding
% is included in the original Butcher table, B.


% extract RK method information from B
s = size(B,2) - 1;    % number of stages
c = B(1:s,1);         % abcissae array
b = B(s+1,2:s+1);     % solution weight array (row)
A = B(1:s,2:s+1);     % coefficients
embedded = (size(B,1) > size(B,2));

% check that the table is valid
if ((abs(c(1)) > 100*eps) || (sum(abs(A(1,:))) > 100*eps))
   error('Invalid input table, first stage is not explicit')
end
if (c(end) > 1+100*eps)
   error('Invalid input table, final abcissa > 1')
end
if (min(c(2:end)-c(1:end-1)) < -100*eps)
   error('Invalid input table, abcissae are not sorted')
end
if (sum(abs(triu(A,1))) > 100*eps)
   error('Invalid input table, must be at most diagonally implicit')
end

% determine whether the table needs padding
if ((norm(A(s,:)-b) > 100*eps) || (abs(c(s)-1) > 100*eps))
   c = [c; 1];
   A = [A, zeros(s,1); b, 0];
end

% construct output G and return (c was already computed)
if embedded
   d = B(s+2,2:s+1);
   G{1} = [A(2:end,:) - A(1:end-1,:); d - A(end-1,:)];
else
   G{1} = A(2:end,:) - A(1:end-1,:);
end

% end of function
