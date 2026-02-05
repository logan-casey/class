function R = recovery_R(kp, zprime, wbar_zprime, par)
%HW.RECOVERY_R
%
% Recovery value for debt holders upon default.
%
% Inputs
%   kp            : next capital choice
%   zprime        : productivity realization
%   wbar_zprime   : default threshold for that z'
%   par           : parameter struct
%
% Output
%   R             : recovery payoff

% operating profit
profit = hw.profit(kp, zprime, par);

% corporate tax during default
% IMPORTANT:
% interest deductibility disappears in default,
% so taxable income excludes interest term.
taxbase = profit - par.delta * kp;

Tc = par.tau_c_pos * max(taxbase,0) + ...
     par.tau_c_neg * min(taxbase,0);

% bankruptcy deadweight loss
bankrupt_cost = par.xi * (1 - par.delta) * kp;

R = (1 - par.delta)*kp + profit - Tc ...
    - bankrupt_cost - wbar_zprime;

end
