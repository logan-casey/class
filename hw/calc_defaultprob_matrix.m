% calc_defaultprob_matrix.m
% Compute PrND(k,b,z) for all (ik,ib,iz) using current eq.rtilde and eq.wbar,
% and save as a .mat file.
%
% Outputs:
%   PrND  : [Nk x Nb x Nz]  probability of non-default next period
%   Pdef  : [Nk x Nb x Nz]  probability of default next period (=1-PrND)
%   ED    : [Nk x Nb x Nz]  expected recovery ratio term E[ 1{D} * R/b ]
%
% Assumes you have loaded: eq, par, grid, zgrid, Pz

Nk = grid.Nk;
Nb = grid.Nb;
Nz = numel(zgrid);

wbar = eq.wbar(:);           % [Nz x 1]
zvec = zgrid(:);             % [Nz x 1]

PrND = nan(Nk, Nb, Nz);
Pdef = nan(Nk, Nb, Nz);
ED   = nan(Nk, Nb, Nz);

% --- Precompute kp-dependent objects ---
kp_vec = grid.kgrid(:);                    % [Nk x 1]
kp_pow_alpha = kp_vec .^ par.alpha;        % [Nk x 1]
one_minus_delta_kp = (1 - par.delta) .* kp_vec;  % [Nk x 1]
delta_kp = par.delta .* kp_vec;            % [Nk x 1]

% profit_mat(izp,ik) = z(izp) * kp(ik)^alpha
profit_mat = zvec * (kp_pow_alpha.');      % [Nz x Nk]

% --- Precompute recovery R(ik,izp) for all kp and z' (does NOT depend on b) ---
% recovery_R signature assumed: hw.recovery_R(kp, zprime, wbar_zprime, par)
R_mat = zeros(Nz, Nk);
for ik = 1:Nk
    kp = kp_vec(ik);
    for izp = 1:Nz
        R_mat(izp, ik) = hw.recovery_R(kp, zvec(izp), wbar(izp), par);
    end
end

for iz = 1:Nz
    prob_row = Pz(iz,:).'; % [Nz x 1], probabilities over z'

    for ik = 1:Nk
        kp = kp_vec(ik);
        profit_vec = profit_mat(:, ik);         % [Nz x 1]
        one_m_d_kp = one_minus_delta_kp(ik);
        d_kp       = delta_kp(ik);

        for ib = 1:Nb
            bp = grid.bgrid(ib);

            % bp <= 0: default prob defined by the same inequality; ED uses b in denom (guarded)
            r_guess = eq.rtilde(ik, ib, iz);

            % realized net worth across z'
            taxbase = profit_vec - d_kp - r_guess * bp;
            Tc = par.tau_c_pos .* max(taxbase,0) + par.tau_c_neg .* min(taxbase,0);

            wreal = one_m_d_kp + profit_vec - Tc - (1 + r_guess) * bp;

            def = (wreal < wbar);          % [Nz x 1]
            nd  = ~def;

            prnd = sum(prob_row .* nd);
            PrND(ik, ib, iz) = prnd;
            Pdef(ik, ib, iz) = 1 - prnd;

            % expected recovery ratio term ED = E[ 1{default} * R/b ]
            denom = max(abs(bp), 1e-12);   % avoid blow-ups if bp ~ 0
            ED(ik, ib, iz) = sum(prob_row .* def .* (R_mat(:,ik) ./ denom));
        end
    end
end



% save('defaultprob_matrix.mat','PrND','Pdef','ED','kp_vec','zvec','wbar','-v7.3');

% Optional quick check for one slice:
% iz0 = ceil(Nz/2);
% imagesc(squeeze(PrND(:,:,iz0))); colorbar; title('PrND(k,b) at mid z'); xlabel('b index'); ylabel('k index');


tot = 0;
deg_pol = 0;
for iw = 1:length(wgrid)
    for iz= 1:length(zgrid)
        ed_i = ED(eq.pol_ik(iw, iz),eq.pol_ib(iw, iz),iz);
        deg_pol = deg_pol + (ed_i==0) + (ed_i==1);
        tot = tot+1;
    end
end
deg_pol/tot