function [zgrid, P, logzgrid] = tauchen(rho, sigma_eps, Nz, m)
%HW.TAUCHEN  Tauchen discretization for log AR(1):
%   log z' = rho log z + eps,  eps ~ N(0, sigma_eps^2)
%
%   [zgrid, P, logzgrid] = hw.tauchen(rho, sigma_eps, Nz, m)
%
% Outputs:
%   logzgrid : Nzx1 grid for log z
%   zgrid    : Nzx1 grid for z = exp(logzgrid)
%   P        : NzxNz transition matrix

if nargin < 4 || isempty(m), m = 4; end
if Nz < 2
    error("Nz must be >= 2.");
end
if abs(rho) >= 1
    error("|rho| must be < 1 for stationary discretization.");
end
if sigma_eps <= 0
    error("sigma_eps must be > 0.");
end

sig_z = sigma_eps / sqrt(1 - rho^2);

logz_min = -m * sig_z;
logz_max =  m * sig_z;

logzgrid = linspace(logz_min, logz_max, Nz)';
step = logzgrid(2) - logzgrid(1);

P = zeros(Nz, Nz);

for i = 1:Nz
    mu = rho * logzgrid(i);

    for j = 1:Nz
        if j == 1
            % (-inf, logz_1 + step/2]
            P(i,j) = normcdf_local((logzgrid(j) - mu + step/2) / sigma_eps);
        elseif j == Nz
            % (logz_N - step/2, +inf)
            P(i,j) = 1 - normcdf_local((logzgrid(j) - mu - step/2) / sigma_eps);
        else
            upper = (logzgrid(j) - mu + step/2) / sigma_eps;
            lower = (logzgrid(j) - mu - step/2) / sigma_eps;
            P(i,j) = normcdf_local(upper) - normcdf_local(lower);
        end
    end
end

% Ensure rows sum to 1 up to numerical tolerance
row_sums = sum(P, 2);
P = P ./ row_sums;

zgrid = exp(logzgrid);
end

function y = normcdf_local(x)
% Normal CDF using erf (no toolboxes)
y = 0.5 * (1 + erf(x ./ sqrt(2)));
end