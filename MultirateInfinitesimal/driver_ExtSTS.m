% Script to perform linear stability analyses for our proposed ExtSTS methods and embeddings
%
% Daniel R. Reynolds

% control flags
maxAlpha = 0;
embedding = true;
plotRK = false;
plotMRI = false;
plotExtSTS = true;

% ARS222 embeddings
driver_ars222(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);        % ImEx-ExtSTS
driver_ars222_erk(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);    % Expl-ExtSTS
driver_ars222_sdirk(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);  % Impl-ExtSTS

% Giraldo embeddings
driver_giraldo(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);       % ImEx-ExtSTS
driver_giraldo_erk(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);   % Expl-ExtSTS
driver_giraldo_dirk(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);  % Impl-ExtSTS

% Other embeddings
driver_ralston(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);       % Expl-ExtSTS
driver_heuneuler(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);     % Expl-ExtSTS
driver_ssp_dirk(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);      % Impl-ExtSTS
driver_erk22a(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);        # Expl-ExtSTS
driver_irk21a(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);        # Impl-ExtSTS


% control flags
maxAlpha = 0;
embedding = false;
plotRK = true;
plotMRI = false;
plotExtSTS = true;

% ARS222 methods
driver_ars222(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);        % ImEx-ExtSTS
driver_ars222_erk(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);    % Expl-ExtSTS
driver_ars222_sdirk(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);  % Impl-ExtSTS

% Giraldo methods
driver_giraldo(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);       % ImEx-ExtSTS
driver_giraldo_erk(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);   % Expl-ExtSTS
driver_giraldo_dirk(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);  % Impl-ExtSTS

% Other methods
driver_ralston(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);       % Expl-ExtSTS
driver_heuneuler(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);     % Expl-ExtSTS
driver_ssp_dirk(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);      % Impl-ExtSTS
driver_erk22a(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);        # Expl-ExtSTS
driver_irk21a(maxAlpha,embedding,plotRK,plotMRI,plotExtSTS);        # Impl-ExtSTS



% end of script