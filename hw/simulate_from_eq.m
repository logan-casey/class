function stats = simulate_from_eq(eq, par, grid, wgrid, zgrid, Pz, Nf, T, L, rngseed)
%SIMULATE_FROM_EQ  Simulate firms using equilibrium policy/value functions
% stats = simulate_from_eq(eq, par, grid, wgrid, zgrid, Pz, Nf, T, rngseed)
%
% Defaults: Nf=20000, T=200, rngseed=1

if nargin < 7 || isempty(Nf), Nf = 20000; end
if nargin < 8 || isempty(T), T = 200; end
if nargin < 9 || isempty(L), L = 14; end
if nargin < 10 || isempty(rngseed), rngseed = 1; end

rng(rngseed);

Nw = length(wgrid);
Nz = length(zgrid);
Nk = grid.Nk; Nb = grid.Nb;

% Interpolants for V over w for each z
F = cell(Nz,1);
for iz = 1:Nz
    F{iz} = griddedInterpolant(wgrid(:), eq.V(:,iz), 'linear', 'nearest');
end

% Precompute stationary distribution of z for initial draws
[V_eig, D] = eig(Pz');
[~, idx] = max(diag(D));
pi = V_eig(:, idx) / sum(V_eig(:, idx));
pi = real(pi);
pi(pi<0) = 0; pi = pi / sum(pi);

% Initialize firm states
ik_now = ones(Nf,1) * round(Nk/2);
ib_now = ones(Nf,1) * round(Nb/2);
% draw initial z using cumulative sampling (avoids randsample dependency)
cpi = cumsum(pi(:))';
u0 = rand(Nf,1);
iz_now = zeros(Nf,1);
for j=1:Nf
    iz_now(j) = find(cpi >= u0(j), 1, 'first');
end

% compute initial wtilde consistent with current k,b,z
k_now = grid.kgrid(ik_now);
b_now = grid.bgrid(ib_now);
rt_now = zeros(Nf,1);
for i=1:Nf
    rt_now(i) = eq.rtilde(ik_now(i), ib_now(i), iz_now(i));
end
% proceeds = b_now;
% pos = b_now > 0;
% proceeds(pos) = b_now(pos) ./ (1 + rt_now(pos) * (1 - par.tau_i));
% compute profit at current z
profit_now = zeros(Nf,1);
for i=1:Nf
    profit_now(i) = hw.profit(k_now(i), zgrid(iz_now(i)), par);
end
taxbase = profit_now - par.delta .* k_now - rt_now .* b_now;
Tc = par.tau_c_pos .* max(taxbase,0) + par.tau_c_neg .* min(taxbase,0);
w_real = (1 - par.delta) .* k_now + profit_now - Tc - (1 + rt_now) .* b_now;
wbar_vec = eq.wbar(:);
w_now = max(wbar_vec(iz_now), w_real);

% storage for last 14 periods
store_inv = zeros(Nf, L);
store_cf = zeros(Nf, L);
store_q = zeros(Nf, L);
store_opinc = zeros(Nf, L);
store_debt_mkt = zeros(Nf, L);
store_inv_rate = zeros(Nf, L);
store_num1 = zeros(Nf, L);
store_num2 = zeros(Nf, L);
store_zk_over_k = zeros(Nf, L);
store_leverage = zeros(Nf, L);
store_equi_flow = zeros(Nf, L);
store_kbook = zeros(Nf, L);

for t = 1:T

    % map current w to nearest wgrid index (per firm)
    % compute absolute distance matrix and take min along columns
    % result: iw_now is Nf x 1
    [~, iw_now] = min(abs(w_now - wgrid'), [], 2);

    % get policy choices (indices)
    idx_lin = iw_now + (iz_now-1) * Nw;
    ik_next = eq.pol_ik(idx_lin);
    ib_next = eq.pol_ib(idx_lin);

    kp = grid.kgrid(ik_next);
    bp = grid.bgrid(ib_next);

    % rtilde depends on ik_next, ib_next, iz_now
    rt = zeros(Nf,1);
    for i=1:Nf
        rt(i) = eq.rtilde(ik_next(i), ib_next(i), iz_now(i));
    end

    % proceeds from bonds (face value bp)
    proceeds = bp;
    pos = bp > 0;
    proceeds(pos) = bp(pos) ./ (1 + rt(pos) * (1 - par.tau_i));

    % net distribution D
    D = w_now + proceeds - kp;

    dist = zeros(Nf,1);
    E = zeros(Nf,1);
    flow = zeros(Nf,1);
    for i=1:Nf
        if D(i) >= 0
            dist(i) = D(i) - hw.tax_dist(D(i), par);
            flow(i) = dist(i);
        else
            E(i) = -D(i);
            flow(i) = -E(i) - hw.equity_cost(E(i), par);
        end
    end

    % operating income: use kp^alpha * current z (following bellman convention)
    opinc = zgrid(iz_now) .* (kp .^ par.alpha);

    % Tobin's q approx: (V(w,z) + (1+rt)*bp) / kp
    Vnow = zeros(Nf,1);
    for i=1:Nf
        Vnow(i) = F{iz_now(i)}(w_now(i));
    end
    q = (Vnow + (1 + rt) .* bp) ./ kp;

    % Debt/Market value real assets
    debt_mkt = bp ./ (Vnow + bp);

    % Investment = kp - (1-par.delta)*k_now
    inv = kp - (1 - par.delta) .* k_now;

    % store if in last L periods
    if t > T - L
        col = t - (T - L);
        store_inv(:, col) = inv ./ kp;
        store_cf(:, col) = flow ./ kp;
        store_q(:, col) = q;
        store_opinc(:, col) = opinc ./ kp;
        store_debt_mkt(:, col) = debt_mkt;
        store_inv_rate(:, col) = inv ./ kp;
        store_num1(:, col) = (zgrid(iz_now) .* (kp.^par.alpha) - hw.tax_dist(zgrid(iz_now) .* (kp.^par.alpha) - par.delta .* kp - rt .* bp, par) - (1 + rt) .* bp) ./ kp;
        store_num2(:, col) = (Vnow + (1 + rt) .* bp) ./ kp;
        store_zk_over_k(:, col) = (zgrid(iz_now) .* (kp.^par.alpha)) ./ kp;
        store_leverage(:, col) = bp ./ (Vnow + bp);
        store_equi_flow(:, col) = ( (E - dist) ) ./ kp; % issuance positive, distributions negative
        store_kbook(:, col) = kp;
    end

    % Draw z' for next period using Markov transition
    rdraw = rand(Nf,1);
    iz_next = zeros(Nf,1);
    cP = cumsum(Pz,2);
    for j=1:Nf
        iz_next(j) = find(cP(iz_now(j),:) >= rdraw(j), 1, 'first');
    end

    % compute next period wtilde: w_real' = (1-delta)*kp + profit(kp, z') - Tc - (1+rt)*bp
    profit_next = zeros(Nf,1);
    for i=1:Nf
        profit_next(i) = hw.profit(kp(i), zgrid(iz_next(i)), par);
    end
    taxbase_next = profit_next - par.delta .* kp - rt .* bp;
    Tc_next = par.tau_c_pos .* max(taxbase_next,0) + par.tau_c_neg .* min(taxbase_next,0);
    w_real_next = (1 - par.delta) .* kp + profit_next - Tc_next - (1 + rt) .* bp;
    wbar_next = eq.wbar(:);
    w_next = max(wbar_next(iz_next), w_real_next);

    % advance states
    k_now = kp;
    b_now = bp;
    w_now = w_next;
    iz_now = iz_next;
    ik_now = ik_next;
    ib_now = ib_next;
end

% Flatten last L periods and compute summary stats
X_inv = store_inv(:);
X_eiss = max(-store_equi_flow(:), 0); % issuance positive
X_dist = max(store_equi_flow(:), 0); % distributions positive
X_inv_var = var(X_inv);

stats = struct();
stats.AverageEquityIssuanceAssets = mean(X_eiss);
stats.VarianceEquityIssuanceAssets = var(X_eiss);
stats.VarianceInvestmentAssets = var(X_inv);
stats.FrequencyEquityIssuance = mean(X_eiss>0);
% payout ratio: mean(distributions / operating income) where opinc>0
opinc_all = store_opinc(:);
dist_all = store_cf(:);
mask = opinc_all > 0;
if any(mask)
    stats.PayoutRatio = mean(dist_all(mask) ./ opinc_all(mask));
else
    stats.PayoutRatio = NaN;
end
stats.FrequencyNegativeDebt = mean(store_leverage(:) < 0);
stats.VarianceDistributions = var(dist_all);
stats.AverageDebtAssetsRatio = mean(store_debt_mkt(:));
stats.CovarianceInvestmentAndEquityIssuance = cov(X_inv, X_eiss); stats.CovarianceInvestmentAndEquityIssuance = stats.CovarianceInvestmentAndEquityIssuance(1,2);
stats.CovarianceInvestmentAndLeverage = cov(X_inv, store_leverage(:)); stats.CovarianceInvestmentAndLeverage = stats.CovarianceInvestmentAndLeverage(1,2);

% serial correlation of income/assets across firms: compute per-firm then average
ser_corrs = zeros(Nf,1);
for i=1:Nf
    series = store_opinc(i, :);
    if std(series) > 0
        ser_corrs(i) = corr(series(1:end-1)', series(2:end)');
    else
        ser_corrs(i) = NaN;
    end
end
stats.SerialCorrelationIncomeAssets = nanmean(ser_corrs);
stats.StdDevIncomeAssets = nanstd(store_opinc(:));

% Print table-like output
fprintf('\nSimulation summary (last %d periods, %d firms):\n', L, Nf);
fn = fieldnames(stats);
for i=1:numel(fn)
    fprintf('  %s: %.4f\n', fn{i}, stats.(fn{i}));
end

end
