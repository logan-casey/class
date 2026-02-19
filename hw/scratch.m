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