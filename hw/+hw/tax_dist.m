function Td = tax_dist(X, par)
%HW.TAX_DIST  Distribution tax T_d(X) for X >= 0:
%   tau_d(x) = tau_d_bar * (1 - exp(-phi x))
%   T_d(X) = integral_0^X tau_d(x) dx
%         = tau_d_bar * [ X - (1 - exp(-phi X))/phi ]
%
% Td = hw.tax_dist(X, par)

if ~isfield(par, "tau_d_bar")
    error("par.tau_d_bar missing");
end
if ~isfield(par, "phi") || isnan(par.phi)
    error("par.phi must be set (not NaN).");
end

X = max(X, 0); % by definition only applies for distributions/payouts

phi = par.phi;
if phi <= 0
    error("phi must be > 0.");
end

Td = par.tau_d_bar .* ( X - (1 - exp(-phi .* X)) ./ phi );
end