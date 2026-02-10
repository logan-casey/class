function rtilde = bond_yield_rtilde(kp, bp, iz, zgrid, Pz, wbar, par, r_for_default)
%HW.BOND_YIELD_RTILDE  Update bond yield using equation (20) (discrete version).
%
% Implements:
%   rtilde = (1/(1-tau_i))* ( ((1 + r*(1-tau_i) - ED)/PrND) - 1 )
% where
%   ED   = sum_{default z'} Q(z,z') * (R(k,z')/b)
%   PrND = sum_{nondefault z'} Q(z,z')
%
% Default classification uses realized net worth:
%   default if w_realized(kp,bp,z') < wbar(z')
% and w_realized uses r_for_default (typically the current guess rtilde_old at (kp,bp,z)).
%
% Inputs
%   kp, bp         : choices (scalar)
%   iz             : current z index
%   zgrid          : [Nz x 1] grid for z
%   Pz             : [Nz x Nz] transition matrix
%   wbar           : [1 x Nz] default boundary wbar(z')
%   par            : parameters
%   r_for_default  : (optional) scalar yield used only to compute default set
%
% Output
%   rtilde         : updated yield

if nargin < 8 || isempty(r_for_default)
    r_for_default = par.r;
end

% If bp <= 0, yield is irrelevant for "bond pricing" in this form.
% (bp<0 corresponds to holding cash/negative debt; set to risk-free.)
if bp <= 0
    rtilde = par.r;
    return;
end

Nz = length(zgrid);

ED = 0.0;
PrND = 0.0;

for izp = 1:Nz
    prob = Pz(iz, izp);
    zprime = zgrid(izp);

    % realized net worth for default test
    w_realized = hw.realized_networth(kp, bp, zprime, r_for_default, par);

    if w_realized < wbar(izp)
        % default region contribution: (R/b)*Q
        R = hw.recovery_R(kp, zprime, wbar(izp), par);
        ED = ED + prob * (R / bp);
    else
        % non-default probability mass
        PrND = PrND + prob;
    end
end

% Safety: avoid divide-by-zero if PrND numerically collapses
if PrND <= 0
    % If default with probability ~1, the conditional formula is not well-defined.
    % Return a very high yield as a signal; could handle separately.
    rtilde = 1e6;
    return;
end

gross_rf_after_tax = 1 + par.r * (1 - par.tau_i);

rtilde = (1/(1 - par.tau_i)) * ( (gross_rf_after_tax - ED)/PrND - 1 );

end