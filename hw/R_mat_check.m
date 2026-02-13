wbar = eq.wbar(:); % [Nz x 1]
Nz = length(zgrid);

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
        - (bankrupt_cost_k') - (wbar);






ik = 8;
ib = 40;
iz = 7;
profit_vec = profit_mat(:, ik);  % [Nz x 1]
one_minus_delta_kp_scalar = one_minus_delta_kp(ik);
delta_kp_scalar = delta_kp(ik);

bp = grid.bgrid(ib);
r_guess = eq.rtilde(ik, ib, iz);

taxbase = profit_vec - delta_kp_scalar - r_guess * bp;
Tc = par.tau_c_pos .* max(taxbase,0) + par.tau_c_neg .* min(taxbase,0);
w_real = one_minus_delta_kp_scalar + profit_vec - Tc - (1 + r_guess) * bp;
