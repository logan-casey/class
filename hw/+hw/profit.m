function y = profit(k, z, par)
%HW.PROFIT  Operating profit function: z * k^alpha
%
% y = hw.profit(k, z, par)
y = z .* (k .^ par.alpha);
end