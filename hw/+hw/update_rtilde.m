function rtilde_new = update_rtilde(rtilde_old, zgrid, Pz, wbar, par, grid, opts)
%HW.UPDATE_RTILDE  Vectorized bond-yield update using discrete eq. (20).
%
% rtilde_new = hw.update_rtilde(rtilde_old, zgrid, Pz, wbar, par, grid, opts)
%
% Implements:
%   rtilde = (1/(1-tau_i)) * ( ((1 + r*(1-tau_i) - ED)/PrND) - 1 )
% where
%   ED   = sum_{default z'} Q(z,z') * (R(k,z')/b)
%   PrND = sum_{nondefault z'} Q(z,z')
%
% Default region is determined by realized net worth:
%   default if w_realized(k,b,z') < wbar(z')
% and w_realized uses rtilde_old(k,b,z) in that test.
%
% Inputs
%   rtilde_old : [Nk x Nb x Nz] current guess
%   zgrid      : [Nz x 1]
%   Pz         : [Nz x Nz]
%   wbar       : [1 x Nz] or [Nz x 1]
%   par, grid  : structs
%   opts       : (optional) struct
%       .omega     damping in (0,1], default 1
%       .cap       max yield cap, default 1e3
%       .minPrND   threshold for PrND, default 1e-12
%
% Output
%   rtilde_new : [Nk x Nb x Nz]

if nargin < 7, opts = struct(); end
if ~isfield(opts,'omega'),   opts.omega = 1.0; end
if ~isfield(opts,'cap'),     opts.cap = 1e3; end
if ~isfield(opts,'minPrND'), opts.minPrND = 1e-12; end

omega   = opts.omega;
cap     = opts.cap;
minPrND = opts.minPrND;

wbar = wbar(:); % [Nz x 1]
Nz = length(zgrid);

rtilde_new = rtilde_old; % initialize

gross_rf_after_tax = 1 + par.r * (1 - par.tau_i);
inv_one_minus_tau_i = 1/(1 - par.tau_i);

% Precompute per-ik objects that depend on kp and z'
kp = grid.kgrid(:);                    % [Nk x 1]
kp_pow_alpha = kp .^ par.alpha;        % [Nk x 1]
one_minus_delta_kp = (1 - par.delta) .* kp;  % [Nk x 1]
delta_kp = par.delta .* kp;            % [Nk x 1]

zgrid_col = zgrid(:);                  % [Nz x 1]

% For recovery R(k,z'): interest does NOT enter, so we can precompute
% R_vec(ik, izp) for all ik,izp:
% profit = z' * kp^alpha
% Tc_def uses taxbase_def = profit - delta*kp
% R = (1-d)kp + profit - Tc_def - xi(1-d)kp - wbar(z')
%
% We'll build profit_mat as [Nz x Nk]: profit_mat(izp,ik)
profit_mat = zgrid_col * (kp_pow_alpha');    % [Nz x Nk]
taxbase_def_mat = profit_mat - (delta_kp');  % [Nz x Nk]

Tc_def_mat = par.tau_c_pos .* max(taxbase_def_mat,0) + ...
             par.tau_c_neg .* min(taxbase_def_mat,0); % [Nz x Nk]

bankrupt_cost_k = par.xi .* (one_minus_delta_kp);     % [Nk x 1]

% R_mat: [Nz x Nk]
R_mat = (one_minus_delta_kp') + profit_mat - Tc_def_mat ...
        - (bankrupt_cost_k') - (wbar);  % wbar broadcasts down rows

% Loop over current z index (outermost): cheap (Nz ~ 15)
for iz = 1:Nz
    prob_row = Pz(iz, :).'; % [Nz x 1]

    % Loop over ik and ib; vectorize over izp inside.
    for ik = 1:grid.Nk
        kp_scalar = kp(ik);

        profit_vec = profit_mat(:, ik);  % [Nz x 1]
        one_minus_delta_kp_scalar = one_minus_delta_kp(ik);
        delta_kp_scalar = delta_kp(ik);

        for ib = 1:grid.Nb
            bp = grid.bgrid(ib);

            if bp <= 0
                % Not meaningful for “bond yield” if bp<=0; set to risk-free
                r_new = par.r;
            else
                r_guess = rtilde_old(ik, ib, iz);

                % Realized net worth w(k,b,z') uses r_guess (for default test)
                taxbase = profit_vec - delta_kp_scalar - r_guess * bp;
                Tc = par.tau_c_pos .* max(taxbase,0) + par.tau_c_neg .* min(taxbase,0);

                w_real = one_minus_delta_kp_scalar + profit_vec - Tc - (1 + r_guess) * bp;

                default_mask = (w_real < wbar);      % [Nz x 1]
                nd_mask = ~default_mask;

                PrND = sum(prob_row .* nd_mask);

                if PrND <= minPrND
                    r_new = cap;
                else
                    ED = sum(prob_row .* default_mask .* (R_mat(:,ik) ./ bp));
                    % exact formula you provided:
                    r_new = inv_one_minus_tau_i * ( (gross_rf_after_tax - ED)/PrND - 1 );
                    if ~isfinite(r_new), r_new = cap; end
                    r_new = min(r_new, cap);
                end
            end

            % damping + default nudge
            rtilde_new(ik, ib, iz) = (1-omega)*rtilde_old(ik, ib, iz) + omega*r_new;
        end
    end
end

end