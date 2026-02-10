function Tc = tax_corp(kp, bp, zprime, rtilde, par)
%HW.TAX_CORP  Corporate tax function on taxable income.
%
% Tax base:
%   y = z' * kp^alpha - delta*kp - rtilde * bp
%
% Taxes:
%   Tc = tau_c_pos * max(y,0) + tau_c_neg * min(y,0)
% where min(y,0) is negative -> this captures "partial loss offset" style
% if tau_c_neg < tau_c_pos.
%
% Tc = hw.tax_corp(kp, bp, zprime, rtilde, par)

if ~isfield(par, "alpha") || isnan(par.alpha)
    error("par.alpha must be set (not NaN).");
end
req = ["delta","tau_c_pos","tau_c_neg"];
for f = req
    if ~isfield(par, f)
        error("par.%s missing", f);
    end
end

if any(kp(:) < 0)
    error("kp must be >= 0.");
end

profit = hw.profit(kp, zprime, par);
taxbase = profit - par.delta .* kp - rtilde .* bp;

Tc = par.tau_c_pos .* max(taxbase, 0) + par.tau_c_neg .* min(taxbase, 0);
end