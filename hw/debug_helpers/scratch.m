[~, def_diag_ns] = hw.update_rtilde(eq.rtilde, zgrid, Pz, eq.wbar, par, grid, struct('smooth_default',false));
fprintf('crosses (no-smooth) = %d / %d\n', nnz(def_diag_ns.crosses), numel(def_diag_ns.crosses));
mg = def_diag_ns.min_gap(:); xg = def_diag_ns.max_gap(:);
fprintf('min_gap: min=%.4g, p10=%.4g, med=%.4g, p90=%.4g, max=%.4g\n', min(mg), prctile(mg,10), median(mg), prctile(mg,90), max(mg));
fprintf('max_gap: min=%.4g, p10=%.4g, med=%.4g, p90=%.4g, max=%.4g\n', min(xg), prctile(xg,10), median(xg), prctile(xg,90), max(xg));


for s = [0.01, 0.05, 0.1, 0.5, 2]
    opts = struct('smooth_default',true,'smooth_s',s,'minPrND',1e-12);
    [~, def_diag_s] = hw.update_rtilde(eq.rtilde, zgrid, Pz, eq.wbar, par, grid, opts);
    pmean = def_diag_s.pdef_mean(:);
    fprintf('s=%.3g: pdef_mean: min=%.4g, p10=%.4g, med=%.4g, p90=%.4g, max=%.4g; share (0<p<1)=%.3f\n', ...
        s, min(pmean), prctile(pmean,10), median(pmean), prctile(pmean,90), max(pmean), mean(pmean>1e-6 & pmean<1-1e-6));
end


% find an example where 0<pdef<1 for s=0.5
opts = struct('smooth_default',true,'smooth_s',0.5,'minPrND',1e-12);
[~, def_diag_example] = hw.update_rtilde(eq.rtilde, zgrid, Pz, eq.wbar, par, grid, opts);
idx = find(def_diag_example.pdef_mean>1e-6 & def_diag_example.pdef_mean<1-1e-6, 1);
if isempty(idx)
    disp('no interior pdef cells found for s=0.5');
else
    [ik, ib, iz] = ind2sub(size(def_diag_example.pdef_mean), idx);
    fprintf('example (ik=%d, ib=%d, iz=%d), pdef_mean=%.4g\n', ik, ib, iz, def_diag_example.pdef_mean(idx));
    % rebuild w_real vector used inside update_rtilde for this triple:
    kp = grid.kgrid(ik);
    bp = grid.bgrid(ib);
    r_guess = eq.rtilde(ik, ib, iz);
    profit_vec = zgrid(:) * (kp^par.alpha);
    taxbase = profit_vec - par.delta*kp - r_guess*bp;
    Tc = par.tau_c_pos .* max(taxbase,0) + par.tau_c_neg .* min(taxbase,0);
    w_real = (1-par.delta)*kp + profit_vec - Tc - (1 + r_guess) * bp;
    disp('zgrid  w_real   wbar   indicator (w_real<wbar)');
    disp([zgrid(:), w_real, eq.wbar(:), (w_real < eq.wbar(:))]);
end




fprintf('V min/max by z: ');
for iz=1:length(zgrid)
    fprintf('[z%d: %.4g/%.4g] \n', iz, min(eq.V(:,iz)), max(eq.V(:,iz)));
end
fprintf('\n');
fprintf('wbar: '); disp(eq.wbar);



% Check Bellman satisfaction at a few states
test_states = [1,1,1; 15,10,5; 30,15,10]; % [iw, iz] pairs
% test_states = [1,1,1; 100,10,5; 120,15,10]; % [iw, iz] pairs
for i=1:size(test_states,1)
    iw = test_states(i,1); iz = test_states(i,2);
    V_current = eq.V(iw, iz);
    
    % Get policy choice
    ik_opt = eq.pol_ik(iw, iz);
    ib_opt = eq.pol_ib(iw, iz);
    
    w = wgrid(iw);
    kp = grid.kgrid(ik_opt);
    bp = grid.bgrid(ib_opt);
    r = eq.rtilde(ik_opt, ib_opt, iz);
    
    % Distribution/issuance flow
    proceeds = bp;
    if bp > 0
        proceeds = bp / (1 + r*(1-par.tau_i));
    end
    D = w + proceeds - kp;
    if D >= 0
        flow = D - hw.tax_dist(D, par);
    else
        flow = -D - hw.equity_cost(-D, par);
    end
    
    % Continuation value (loop over izp)
    EV = 0;
    for izp = 1:length(zgrid)
        profit = zgrid(izp) * (kp^par.alpha);
        taxbase = profit - par.delta*kp - r*bp;
        Tc = par.tau_c_pos*max(taxbase,0) + par.tau_c_neg*min(taxbase,0);
        wreal = (1-par.delta)*kp + profit - Tc - (1+r)*bp;
        w_next = max(eq.wbar(izp), wreal);
        V_next = interp1(wgrid, eq.V(:,izp), w_next, 'linear', 'extrap');
        EV = EV + Pz(iz,izp) * V_next;  % correct: just transition weight
    end
    
    beta = 1/(1 + par.r*(1-par.tau_i));
    V_bellman = flow + beta * EV;
    
    fprintf('(iw=%d,iz=%d): V_curr=%.4g, V_bellman=%.4g, diff=%.4g\n', iw, iz, V_current, V_bellman, V_current - V_bellman);
end





% At this state, what does rtilde look like as b increases?
ik = 17; iz = 1;
fprintf('rtilde schedule for (ik=%d, iz=%d):\n', ik, iz);
fprintf('ib    b      rtilde   monotone_ok?\n');
last_r = -Inf;
for ib = 1:grid.Nb
    r = eq.rtilde(ik, ib, iz);
    ok = (r >= last_r - 1e-10);
    fprintf('%2d  %.4g  %.6g  %d\n', ib, grid.bgrid(ib), r, ok);
    % if ~ok, break; end
    last_r = r;
end



% Inspect PrND and ED across b for one (ik,iz)
ik = 17; iz = 1;
kp = grid.kgrid(ik);
Nz = length(zgrid);
Nb = grid.Nb;

% Precompute profit and R_mat as in update_rtilde
profit_vec = zgrid(:) * (kp^par.alpha);            % [Nz x 1]
taxbase_def = profit_vec - par.delta * kp;
Tc_def = par.tau_c_pos .* max(taxbase_def,0) + par.tau_c_neg .* min(taxbase_def,0);
bankrupt_cost_k = par.xi .* ((1-par.delta) * kp);
R_vec = (1-par.delta)*kp + profit_vec - Tc_def - bankrupt_cost_k - eq.wbar(:); % [Nz x 1]

PrND = nan(Nb,1);
ED   = nan(Nb,1);
rvec = squeeze(eq.rtilde(ik, :, iz));   % [1 x Nb] -> [Nb x1]
bvec = grid.bgrid(:);

for ib = 1:Nb
    bp = bvec(ib);
    r = rvec(ib);

    taxbase = profit_vec - par.delta * kp - r * bp;
    Tc = par.tau_c_pos .* max(taxbase,0) + par.tau_c_neg .* min(taxbase,0);
    w_real = (1-par.delta)*kp + profit_vec - Tc - (1 + r) * bp;

    default_mask = (w_real < eq.wbar(:));       % [Nz x 1]
    PrND(ib) = sum(Pz(iz,:)' .* (~default_mask));
    if bp == 0
        ED(ib) = 0;
    else
        ED(ib) = sum(Pz(iz,:)' .* default_mask .* (R_vec ./ max(bp,1e-12)));
    end
end

T = table((1:Nb)', bvec, rvec', PrND, ED, 'VariableNames', {'ib','b','rtilde','PrND','ED'});
disp(T);

% quick plots
figure;
subplot(2,1,1); plot(bvec, rvec, '-o'); ylabel('rtilde'); grid on;
subplot(2,1,2); yyaxis left; plot(bvec, PrND, '-o'); ylabel('PrND'); yyaxis right; plot(bvec, ED, '-x'); ylabel('ED'); grid on;
saveas(gcf, 'prnd.png', 'png')