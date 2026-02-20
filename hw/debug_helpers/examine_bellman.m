
for ib = 1:grid.Nb %******************

iz = 13;
iw = 18;

ik = eq.pol_ik(iw,iz); % optimal k
wbar = eq.wbar;
rtilde = eq.rtilde;
V = eq.V;

wtilde = wgrid(iw);
prob_row = Pz(iz, :).';        % [Nz x 1]
zvec     = zgrid(:);           % [Nz x 1]
wbarvec  = wbar(:);            % [Nz x 1]
Nz = length(zvec);

F = cell(Nz,1);
for izp = 1:Nz
    F{izp} = griddedInterpolant(wgrid, V(:,izp), 'linear', 'nearest');
end

Vbest   = -Inf;
ik_best = NaN;
ib_best = NaN;

kp = grid.kgrid(ik);

profit_vec = zvec .* (kp^par.alpha);                % [Nz x 1]
one_minus_delta_kp = (1 - par.delta) * kp;
delta_kp = par.delta * kp;



bp = grid.bgrid(ib);

rt = rtilde(ik, ib, iz);

% Flow to equity today via net distribution D = wtilde + bp - kp
% D = wtilde + bp - kp;
% Treat bp as face value; proceeds are discounted
proceeds = bp;
if bp > 0
    proceeds = bp / (1 + rt*(1 - par.tau_i));
end
D = wtilde + proceeds - kp;

if D >= 0
    flow = D - hw.tax_dist(D, par);
else
    E = -D;
    flow = -E - hw.equity_cost(E, par);
end

% Realized net worth vector for all z'
% taxable income uses rt*bp
taxbase_vec = profit_vec - delta_kp - rt * bp;
Tc_vec = par.tau_c_pos .* max(taxbase_vec, 0) + par.tau_c_neg .* min(taxbase_vec, 0);

w_real_vec = one_minus_delta_kp + profit_vec - Tc_vec - (1 + rt) * bp;

% Revised net worth wtilde' = max(wbar(z'), w_real)
w_next_vec = max(wbarvec, w_real_vec);

% Continuation values V(w_next, z')
Vnext_vec = zeros(Nz, 1);
for izp = 1:Nz
    Vnext_vec(izp) = F{izp}(w_next_vec(izp));
end

EV = prob_row' * Vnext_vec;

beta = 1 / (1 + par.r * (1 - par.tau_i));
Vcand = flow + beta * EV;


disp(Vcand)
end

disp(eq.pol_ib(iw,iz))