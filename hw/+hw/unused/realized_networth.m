function w = realized_networth(kp, bp, zprime, rtilde, par)
%HW.REALIZED_NETWORTH  Compute realized net worth w after production shock.
%
% w = (1-delta)kp + z' kp^alpha - Tc(kp,bp,z',rtilde) - (1+rtilde)*bp
%
% w = hw.realized_networth(kp, bp, zprime, rtilde, par)

Tc = hw.tax_corp(kp, bp, zprime, rtilde, par);
profit = hw.profit(kp, zprime, par);

w = (1 - par.delta).*kp + profit - Tc - (1 + rtilde).*bp;
end