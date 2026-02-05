function Phi = equity_cost(E, par)
%HW.EQUITY_COST  Equity issuance cost Phi(E).
% Assumption: if E > 0 (issuance), cost is lambda0 + lambda1 E + lambda2 E^2
% else 0.
%
% Phi = hw.equity_cost(E, par)

req = ["lambda0","lambda1","lambda2"];
for f = req
    if ~isfield(par, f) || isnan(par.(f))
        error("par.%s must be set (not NaN).", f);
    end
end

Epos = max(E, 0);
Phi  = zeros(size(Epos));

idx = (Epos > 0);
Phi(idx) = par.lambda0 + par.lambda1 .* Epos(idx) + par.lambda2 .* (Epos(idx).^2);
end