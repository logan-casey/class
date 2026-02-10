function par = params(sample)
%HW.PARAMS  Load calibrated + estimated parameters for Hennessy-Whited.
%
%   par = hw.params(sample)
%
% Inputs
%   sample : string/char, e.g. "full", "small", "large"
%
% Output
%   par : struct with model parameters + numerical settings
%
% Notes
% - Calibrated constants are set here (edit as needed to match the paper).
% - Estimated parameters must be filled from the paper's tables (left as NaN).

if nargin < 1 || isempty(sample)
    sample = "full";
end
sample = string(lower(sample));

par = struct();

% --------------------
% Calibrated constants
% --------------------
% (These are the standard numbers used in the HW calibration section.
%  Please overwrite if your PDF reports different values.)
par.r          = 0.025;  % risk-free rate
par.delta      = 0.15;   % depreciation
par.tau_i      = 0.29;   % personal tax rate on interest
par.tau_c_pos  = 0.40;   % corporate tax rate on positive income
par.tau_c_neg  = 0.20;   % corporate tax rate on negative income
par.tau_d_bar  = 0.12;   % max marginal distribution tax parameter

% --------------------
% Estimated parameters
% --------------------
% Fill these from the paper's parameter table corresponding to 'sample'.
switch sample
    case "full"
        par.alpha     = 0.627;  % production curvature
        par.phi       = 0.732;  % distribution tax curvature
        par.rho       = 0.684;  % AR(1) persistence for log z
        par.sigma_eps = 0.118;  % AR(1) innovation sd for log z

        % Equity issuance cost parameters (Assumption 4)
        par.lambda0   = 0.598;
        par.lambda1   = 0.091;
        par.lambda2   = 0.0004;

        % Bankruptcy deadweight parameter
        par.xi        = 0.104;

    case "small"
        par.alpha     = 0.693;
        par.phi       = 0.831;
        par.rho       = 0.498;
        par.sigma_eps = 0.159;
        par.lambda0   = 0.951;
        par.lambda1   = 0.120;
        par.lambda2   = 0.0004;
        par.xi        = 0.151;

    case "large"
        par.alpha     = 0.577;
        par.phi       = 0.695;
        par.rho       = 0.791;
        par.sigma_eps = 0.086;
        par.lambda0   = 0.389;
        par.lambda1   = 0.053;
        par.lambda2   = 0.0002;
        par.xi        = 0.084;

    otherwise
        error("Unknown sample='%s'.", sample);
end

% --------------------
% Numerical settings
% --------------------
par.Nz_solve = 15;
par.Nz_sim   = 60;
par.tauchen_m = 4;     % support width in std deviations (m in Tauchen)

% k-grid construction defaults (paper uses kbar*(1-delta)^[0,0.5,1,...,15])
par.k_exponents = 0:0.5:15;

% Suggested baseline grid sizes
par.Nk = numel(par.k_exponents);
par.Nb = par.Nk/2;     % paper: b grid has half as many points as k grid
% note should be integer

% Basic validation hook (optional)
par = hw_validate_params_basic(par);
end

function par = hw_validate_params_basic(par)
% Keep it lightweight: only check calibrated fields exist and are sensible.
required = ["r","delta","tau_i","tau_c_pos","tau_c_neg","tau_d_bar","Nz_solve","Nz_sim"];
for f = required
    if ~isfield(par, f)
        error("Missing par.%s", f);
    end
end
if par.delta <= 0 || par.delta >= 1
    error("delta must be in (0,1).");
end
if par.r < 0
    error("r must be >= 0.");
end
end