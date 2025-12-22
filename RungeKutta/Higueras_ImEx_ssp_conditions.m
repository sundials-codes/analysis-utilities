function Higueras_ImEx_ssp_conditions(ERKtableName, DIRKtableName)

% This function takes an explicit Butcher table and its corresponding
% implicit Butcher table and determines if the IMEX scheme satifies the
% 'L', 'S', 'P', 'U', 'M' properties in Higueras paper:
% "In Optimized strong stability preserving IMEX Runge–Kutta methods by
% Higueras et.al (2014)"
% 'L': the implicit method is L-stable
% 'S': the stability region for the explicit part contains an interval on the imaginary axis
% 'P': the amplification factor g for the implicit method is always positive
% 'U': the IMEX RK method features uniform convergence
% 'M': the IMEX RK method has a nontrivial region of absolute monotonicity.
% Input: ERKtableName  - name of the explicit method in butcher.m
%        DIRKtableName - name of the implicit method in butcher.m
%------------------------------------------------------------
% Programmer(s): Sylvia Amihere @ UMBC
%------------------------------------------------------------
% Copyright (c) 2025, University of Maryland Baltimore County.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------


% ------ Butcher tables for explicit and implicit methods ------
B_explicit = butcher(ERKtableName,true);
B_implicit = butcher(DIRKtableName,true);

% ------ number of stages of methods ------
stages_ex = length(B_explicit(1, 2:end));
stages_im = length(B_implicit(1, 2:end));
if (stages_ex ~= stages_im)
    error('The number of stages for both explicit and implcit methods should be equal.\n');
else
    fprintf("Both method have %d stages\n", stages_ex);
end

% ------ determine if the methods have an embedding ------
[rows_ex, cols_ex] = size(B_explicit);
[rows_im, cols_im] = size(B_implicit);

if (stages_ex + 2 == rows_ex)
    new_rows_ex = rows_ex - 1; %remove the embedding
    fprintf('The explicit method has an embedding.\n');
elseif (stages_ex + 1 == rows_ex)
    new_rows_ex = rows_ex;
    fprintf('The explicit method has no embedding.\n');
end

if (stages_im + 2 == rows_im)
    new_rows_im = rows_im - 1; %remove the embedding
    fprintf('The implicit method has an embedding.\n');
elseif (stages_im + 1 == rows_im)
    new_rows_im = rows_im;
    fprintf('The implicit method has no embedding.\n');
end

% ---- extract Butcher table for explcit method without embedding -----
butcher_ex = B_explicit(1:new_rows_ex, 1:cols_ex);
A_ex       = butcher_ex(1:end-1, 2:end); %matrix A from explicit method
b_ex       = butcher_ex(new_rows_ex, 2:end); %vector b from explicit method
c_ex       = butcher_ex(1:end-1, 1); % vector c from explicit method

% ---- extract Butcher table for implcit method without embedding -----
butcher_im = B_implicit(1:new_rows_im, 1:cols_im);
Atilde_im  = butcher_im(1:end-1, 2:end);  %matrix A from implicit method
btilde_im  = butcher_im(new_rows_ex, 2:end); %vector b from explicit method
ctilde_im  = butcher_im(1:end-1, 1); % vector c from explicit method


% -------------------------------------------------------------------------
%                   Check the L S P U M properties
% -------------------------------------------------------------------------

% ================ Check 'L' property using stab_function.m ===============
[alpha, beta] = stab_function(Atilde_im,btilde_im);
deg_alpha     = length(alpha); % number of coefficients of the numerator polynomial (alpha)
deg_beta      = length(beta);  % number of coefficients of the denominator polynomial (beta)

if (deg_alpha < deg_beta)
    fprintf('Implicit method is L-stable. Hence, "L" property is satisfied.\n')
elseif (deg_alpha == deg_beta)
    fprintf('Implicit method is A-stable but not L-stable. Hence, "L" property is not satisfied.\n')
else %deg_alpha > deg_beta
    fprintf('Implicit method is not A-stable. Hence, "L" property is not satified.\n') % hence not L-stable
end


% ================= Check 'S' property using stab_function.m ==============
[alpha_ex, beta_ex] = stab_function(A_ex,b_ex); %coefficients of rational polynomial
omega_vals          = logspace(-8,0);           
s_property          = false;

for ik = 1:length(omega_vals)
    num_val   = 0.0 + 0.0i;
    denom_val = 0.0 + 0.0i;

    for k = 1:length(alpha_ex)
        num_val = num_val + alpha_ex(k)*(1i * omega_vals(ik))^(k-1); %value of the numerator polynomial at the imaginary value
    end

    for k = 1:length(beta_ex)
        denom_val = denom_val + beta_ex(k)*(1i * omega_vals(ik))^(k-1); %value of the denominator polynomial at the imaginary value
    end

    stab_val = num_val/denom_val;
    if abs(stab_val) <= 1 % is absolute value of rational polynomial at the imaginary value is less than or equal to 1
        s_property = true;
        % break;  %exit the loop once an omega value results in an interval on the imaginary axis for the explicit method
    end
end

if (s_property)
    fprintf('Stability region for the explicit part contains an interval on the imaginary axis. Hence, "S" property is satisfied.\n');
else
    fprintf('Stability region for the explicit part DOES NOT contain an interval on the imaginary axis. Hence, "S" property is not satisfied.\n');
end


% =========================== Check 'P' property ==========================
[alpha_im, beta_im] = stab_function(Atilde_im,btilde_im); %coefficients of rational polynomial
omega_vals          = logspace(-8,0);          
p_property          = false;

for ik = 1:length(omega_vals)
    num_val   = 0.0 + 0.0i;
    denom_val = 0.0 + 0.0i;

    for k = 1:length(alpha_im)
        num_val = num_val + alpha_im(k)*(1i * omega_vals(ik))^(k-1); %value of the numerator polynomial at the imaginary value
    end

    for k = 1:length(beta_im)
        denom_val = denom_val + beta_im(k)*(1i * omega_vals(ik))^(k-1); %value of the denominator polynomial at the imaginary value
    end

    stab_val = num_val/denom_val;
    if abs(stab_val) > 0 % is absolute value of rational polynomial at the imaginary value is less than or equal to 1
        p_property = true;
        % break;  %exit the loop once an omega value results in an interval on the imaginary axis for the explicit method
    end
end

if (p_property)
    fprintf('Amplification factor for the implicit method is positive. Hence, "P" property is satisfied.\n');
else
    fprintf('Amplification factor for the implicit method is NOT positive. Hence, "P" property is not satisfied.\n');
end


% ======================== Check 'U' property =============================
uc_val = btilde_im*(Atilde_im\c_ex);
if (uc_val==1)
    fprintf(['ImEx method is uniformly convergent. ' ...
        'Hence, "U" property is satisfied.\n'])
else
    fprintf(['ImEx method is NOT uniformly convergent. ' ...
        'Hence, "U" property is not satisfied.\n'])
end


% ========================= Check 'M' property ============================
r1_vals    = logspace(-8,0);
r2_vals    = logspace(-8,0);
e          = ones(stages_ex,1);
m_property = false;

% I + r1*A + r2*Atilde is nonsigular, inv(I + r1*A + r2*Atilde)*e >=0
% inv(I + r1*A + r2*Atilde)*A >=0,    inv(I + r1*A + r2*Atilde)*Atilde >=0
for r1 = r1_vals
    for r2 = r2_vals
        matrixM = eye(stages_ex) + r1 * A_ex + r2 * Atilde_im;
        if (det(matrixM)~=0)
            if (all(matrixM\e >= 0, 'all') && all(matrixM\A_ex >= 0, 'all') && all(matrixM\Atilde_im >= 0, 'all'))
                m_property = true;
                % fprintf('ImEx method has a nontrivial region of absolute monotonicity at (%.8f, %.8f). Hence, "M" property is satisfied.\n', r1,r2);
                % break;
            end
        end
    end
end

if (m_property)      
    fprintf('ImEx method has a nontrivial region of absolute monotonicity. Hence, "M" property is satisfied.\n');
else
    fprintf('ImEx method DOES NOT have a nontrivial region of absolute monotonicity. Hence, "M" property is not satisfied.\n');
end

end %of function
