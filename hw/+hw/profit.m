function y = profit(k, z, par)
%HW.PROFIT  Operating profit function: z * k^alpha
%
% y = hw.profit(k, z, par)

if any(k(:) < 0)
    error("k must be >= 0.");
end
if ~isfield(par, "alpha") || isnan(par.alpha)
    error("par.alpha must be set (not NaN).");
end

y = z .* (k .^ par.alpha);
end