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
btilde_im  = butcher_im(new_rows_im, 2:end); %vector b from implicit method
ctilde_im  = butcher_im(1:end-1, 1); % vector c from implicit method

% ====== Check stability property of implicit method with embedding =======
y_vals   = [0, logspace(-4, 6, 2000)];
if (stages_im + 2 == rows_im) %embedding exists
    b_embedding = B_implicit(end, 2:end);
    [alpha_s1, beta_s1] = stab_function(Atilde_im, b_embedding);
    alpha_s1 = double(alpha_s1);
    beta_s1  = double(beta_s1);

    deg_alpha_s1 = find(abs(alpha_s1) >= eps, 1, 'last') - 1;
    deg_beta_s1  = find(abs(beta_s1)  >= eps, 1, 'last') - 1;

    A_stable_s1 = true;
    for ik = 1:length(y_vals)
        z = 1i * y_vals(ik);
        R = polyval(flip(alpha_s1), z) / polyval(flip(beta_s1), z);
        if abs(R) > 1 + 1e-12
            A_stable_s1 = false;
        end
    end
    poles = roots(flip(beta_s1(1:deg_beta_s1+1)));
    if any(real(poles) <= 0)
        A_stable_s1 = false;
    end

    if (A_stable_s1 && deg_alpha_s1 < deg_beta_s1)
        fprintf('Implicit method with embedding is L-stable.\n')
    elseif (A_stable_s1)
        fprintf('Implicit method with embedding is A-stable but not L-stable.\n')
    else
        fprintf('Implicit method with embedding is not A-stable.\n')
    end
end

[alpha_s1, beta_s1] = stab_function(Atilde_im,b_embedding);

%keep coefficeints above machine precision
alpha_s1_amcp = abs(alpha_s1)>=eps;
beta_s1_amcp  = abs(beta_s1)>=eps;
alpha_s1      = alpha_s1(alpha_s1_amcp);
beta_s1       = beta_s1(beta_s1_amcp);

deg_alpha_s1 = length(alpha_s1); % number of coefficients of the numerator polynomial (alpha)
deg_beta_s1 = length(beta_s1);  % number of coefficients of the denominator polynomial (beta)

if (deg_alpha_s1 < deg_beta_s1)
    fprintf('Implicit method with embedding is L-stable. \n')
elseif (deg_alpha_s1 == deg_beta_s1)
    fprintf('Implicit method with embedding is A-stable but not L-stable. \n')
else %deg_alpha_s1 > deg_beta_s1
    fprintf('Implicit method with embedding is not L-stable.\n') 
end

fprintf("\n")
fprintf("=================Checking LSPUM Properties=====================\n")

% -------------------------------------------------------------------------
%                   Check the L S P U M properties
% -------------------------------------------------------------------------

% ================ Check 'L' property using stab_function.m ===============
[alpha, beta] = stab_function(Atilde_im,btilde_im);
alpha = double(alpha);
beta = double(beta);

% degree = position of the LAST nonzero coefficient
deg_alpha = find(abs(alpha) >= eps, 1, 'last') - 1; % number of coefficients of the numerator polynomial (alpha)
deg_beta  = find(abs(beta)  >= eps, 1, 'last') - 1; % number of coefficients of the denominator polynomial (beta)

%check for A-stability: |R(iy)| <= 1 on the imaginary axis
y_vals   = [0, logspace(-4, 6, 2000)];
A_stable = true;
for ik = 1:length(y_vals)
    z = 1i * y_vals(ik);
    R = polyval(flip(alpha), z) / polyval(flip(beta), z);
    if abs(R) > 1 + 1e-12
        A_stable = false;
    end
end

% A-stability, all poles must be in the right half-plane
poles = roots(flip(beta(1:deg_beta+1)));%flipping because in stab function polynomial degrees go from lowest to highest
if any(real(poles) <= 0)
    A_stable = false;
end

if (A_stable && deg_alpha < deg_beta)
    fprintf('Implicit method is L-stable. Hence, "L" property is satisfied.\n')
elseif (A_stable)
    fprintf('Implicit method is A-stable but not L-stable. Hence, "L" property is not satisfied.\n')
else 
    fprintf('Implicit method is not L-stable (also not A-stable). Hence, "L" property is not satified.\n') % hence not L-stable
end


% ================= Check 'S' property using stab_function.m ==============
[alpha_ex, beta_ex] = stab_function(A_ex,b_ex); %coefficients of rational polynomial
alpha_ex = double(alpha_ex); %real numbers for polyval
beta_ex = double(beta_ex);
omega_vals          = linspace(1e-2, 5, 5000);           
w_max = 0; %size of the interval so far, initialising

for ik = 1:length(omega_vals)
    z = 1i * omega_vals(ik);
    R = polyval(flip(alpha_ex), z) / polyval(flip(beta_ex), z);
    if abs(R) <= 1
        w_max = omega_vals(ik);   
    else
        break;                    % first point outside
    end
end

s_property = (w_max > 0);

if (s_property)
    fprintf('Explicit stability region contains [-iw, iw] with w = %.4f. Hence, "S" property is satisfied.\n', w_max);
else
    fprintf('Explicit stability region DOES NOT contain an interval on the imaginary axis. Hence, "S" property is not satisfied.\n');
end


% =========================== Check 'P' property (R(z)>=0) ==========================
[alpha_im, beta_im] = stab_function(Atilde_im,btilde_im); %coefficients of rational polynomial
alpha_im = double(alpha_im); %real numbers for polyval
beta_im = double(beta_im);
z_vals  = -logspace(-8,8,2000); 
p_property = true; 
for ik = 1:length(z_vals)
    R = polyval(flip(alpha_im), z_vals(ik)) / polyval(flip(beta_im), z_vals(ik));
    if R < 0
        p_property = false;   % found a negative value, hence P fails
        break;
    end
end

if (p_property)
    fprintf('Amplification factor for/Stability function of the implicit method is non-negative for all z <= 0. Hence, "P" property is satisfied.\n');
else
    fprintf('Amplification factor for/Stability function of the implicit method becomes negative for some z < 0.. Hence, "P" property is not satisfied.\n');
end


% ======================== Check 'U' property =============================
uc_val = btilde_im*(Atilde_im\c_ex);
if  abs(uc_val - 1) < 1e-12
    fprintf(['ImEx method is uniformly convergent. ' ...
        'Hence, "U" property is satisfied.\n'])
else
    fprintf(['ImEx method is NOT uniformly convergent. ' ...
        'Hence, "U" property is not satisfied.\n'])
end


% ========================= Check 'M' property ============================
s = stages_ex;
bigA =  double([A_ex, zeros(s,1); b_ex, 0]); %A with b added as last row
bigAt = double([Atilde_im, zeros(s,1); btilde_im, 0]);
r1_vals    = logspace(-8,0);
r2_vals    = logspace(-8,0);
e          = ones(s+1,1);
tol_m = 1e-14;
m_property = false;

% I + r1*A + r2*Atilde is nonsigular, inv(I + r1*A + r2*Atilde)*e >=0
% inv(I + r1*A + r2*Atilde)*A >=0,    inv(I + r1*A + r2*Atilde)*Atilde >=0
for r1 = r1_vals
    for r2 = r2_vals
        matrixM = eye(s+1) + r1 * bigA + r2 * bigAt;
        if (det(matrixM)~=0)
            if (all(matrixM\e >= -tol_m, 'all') && ...
                    all(matrixM\bigA >= -tol_m, 'all') && ...
                    all(matrixM\bigAt >= -tol_m, 'all'))
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
